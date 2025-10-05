
module write_utils

    use parameters
    use geometry_utils
    use simulation_state
    use, intrinsic :: iso_fortran_env, only: real64

    implicit none

contains

    subroutine WriteLAMMPSTRJ(box, lammpstrj_filename, append_mode)

        implicit none

        ! Input parameters
        logical, OPTIONAL :: append_mode    ! Flag to append to existing files
        type(type_box), intent(inout) :: box
        character(len=*), intent(in), optional :: lammpstrj_filename ! Output trajectory in lammpstrj format

        ! Local variables
        integer :: i, j, k                      ! Loop indices over residue types, residues, and atoms
        integer :: atom_id, atom_type
        integer :: UNIT_LMP = 18
        logical :: do_append                    ! Internal flag controlling append/overwrite
        integer :: dim                          ! Integer for looping over dimensions
        character(len=20) :: file_position      ! 'APPEND' or 'ASIS' for file positioning
        real(real64), dimension(3) :: pos, com  ! Center of mass and atomic position vectors

        ! Determine whether to append or overwrite
        do_append = .false.
        if (PRESENT(append_mode)) do_append = append_mode
        if (do_append) then
            file_position = 'APPEND'
        else
            file_position = 'ASIS'
        end if

        ! Open file for writing
        open(UNIT=UNIT_LMP, FILE=TRIM(output_path) // lammpstrj_filename, &
            STATUS='UNKNOWN', ACTION='write', POSITION=file_position)

        ! Write LAMMPS-style header
        write(UNIT_LMP, '(A)') "ITEM: TIMESTEP"
        write(UNIT_LMP, '(I10)') input%nb_block
        write(UNIT_LMP, '(A)') "ITEM: NUMBER OF ATOMS"
        write(UNIT_LMP, '(I10)') box%num_atoms
        write(UNIT_LMP, '(A)') "ITEM: BOX BOUNDS pp pp pp"
        write(UNIT_LMP, '(F15.8,1X,F15.8)') -box%matrix(1, 1)/2, box%matrix(1, 1)/2
        write(UNIT_LMP, '(F15.8,1X,F15.8)') -box%matrix(2, 2)/2, box%matrix(2, 2)/2
        write(UNIT_LMP, '(F15.8,1X,F15.8)') -box%matrix(3, 3)/2, box%matrix(3, 3)/2

        ! Atom data header
        write(UNIT_LMP, '(A)') "ITEM: ATOMS id type x y z"

        atom_id = 0
        do i = 1, nb%type_residue
            do j = 1, box%num_residues(i)

                ! Extract CoM
                do dim = 1, 3
                    com(dim) = box%mol_com(dim, i, j)
                end do

                ! Wrap CoM into box for active molecules
                if (input%is_active(i) == 1) then
                    call WrapIntoBox(com, box)
                end if

                do k = 1, nb%atom_in_residue(i)

                    atom_id = atom_id + 1
                    atom_type = box%atom_types(i, k)

                    do dim = 1, 3
                        pos(dim) = com(dim) + box%site_offset(dim, i, j, k)
                    end do

                    ! Wrap position for inactive structure
                    if (input%is_active(i) == 0) then
                        call WrapIntoBox(pos, box)
                    end if

                    write(UNIT_LMP, '(I6,1X,I4,3(1X,F12.7))') atom_id, atom_type, pos(1), pos(2), pos(3)
                end do
            end do
        end do

        close(UNIT_LMP)

    end subroutine WriteLAMMPSTRJ

    subroutine WriteEnergyAndCount()

        implicit none

        integer :: UNIT_ENERGY = 18
        integer :: UNIT_COUNT  = 19
        integer :: UNIT_MOVES  = 20
        character(len=8) :: file_status
        character(len=:), allocatable :: filename
        integer :: resi

        real(real64) :: e_recip_coulomb_kcalmol
        real(real64) :: e_non_coulomb_kcalmol
        real(real64) :: e_coulomb_kcalmol
        real(real64) :: e_ewald_self_kcalmol
        real(real64) :: e_intra_coulomb__kcalmol
        real(real64) :: e_total_kcalmol

        ! Decide whether to create new files or append
        ! For current_block == 0, we want to recreate the file from scratch
        ! For later blocks, we append to existing files
        if (current_block == 0) then
            file_status = 'REPLACE'
        else
            file_status = 'OLD'
        end if

        ! Convert energy components to kcal/mol
        e_recip_coulomb_kcalmol = energy%recip_coulomb * KB_kcalmol
        e_non_coulomb_kcalmol = energy%non_coulomb * KB_kcalmol
        e_coulomb_kcalmol = energy%coulomb * KB_kcalmol
        e_ewald_self_kcalmol = energy%ewald_self * KB_kcalmol
        e_intra_coulomb__kcalmol = energy%intra_coulomb * KB_kcalmol
        e_total_kcalmol = energy%total * KB_kcalmol

        ! ------------------------------------------------------
        ! Write to energy.dat
        ! ------------------------------------------------------
        open(UNIT=UNIT_ENERGY, FILE=TRIM(output_path) // 'energy.dat', &
            STATUS=file_status, ACTION='write', POSITION='APPEND')
        if (current_block == 0) then
            write(UNIT_ENERGY, '(A)') '#    block        total        recipCoulomb' // &
                                    '     non-coulomb      coulomb     ewald_self    intramolecular-coulomb'
        end if
        write(UNIT_ENERGY, '(I10, 1X, F16.6, 1X, F16.6, 1X, F16.6, 1X, F16.6, 1X, F16.6, 1X, F16.6)') &
            current_block, e_total_kcalmol, e_recip_coulomb_kcalmol, e_non_coulomb_kcalmol, &
            e_coulomb_kcalmol, e_ewald_self_kcalmol, e_intra_coulomb__kcalmol
        close(UNIT_ENERGY)

        ! ------------------------------------------------------
        ! Loop over residues and write to number_RESNAME.dat
        ! ------------------------------------------------------
        do resi = 1, nb%type_residue

            if (primary%num_residues(resi) /= 0 .and. input%is_active(resi) == 1) then
                ! Construct the filename for this residue
                filename = TRIM(output_path) // 'number_' // TRIM(res%names_1d(resi)) // '.dat'

                ! Open file for append (create if not exists)
                open(UNIT=UNIT_COUNT, FILE=filename, STATUS=file_status, ACTION='write', POSITION='APPEND')

                ! Write header only once (when current_block == 0)
                if (current_block == 0) then
                    write(UNIT_COUNT, '(A)') '# Block   Active_Molecules'
                end if

                ! Write the block number and count for this residue
                write(UNIT_COUNT, '(I10, 1X, I10)') current_block, primary%num_residues(resi)

                close(UNIT_COUNT)
            end if

        end do

        ! ------------------------------------------------------
        ! Write to moves.dat
        ! ------------------------------------------------------
        open(UNIT=UNIT_MOVES, FILE=TRIM(output_path) // 'moves.dat', &
            STATUS=file_status, ACTION='write', POSITION='APPEND')
        if (current_block == 0) then
            write(UNIT_MOVES, '(A)') &
                '# Block   Trial_Trans   Trans_Moves   Trial_Create   Create_Moves   ' // &
                'Trial_Delete   Delete_Moves   Trial_Rotate   Rotate_Moves   ' // &
                'Trial_BigMove   Big_Moves'
        end if
        write(UNIT_MOVES, '(I12, 1X, I12, 1X, I12, 1X, I12, 1X, I12, 1X, I12, 1X,' // &
                          'I12, 1X, I12, 1X, I12)') &
            current_block, &
            counter%trial_translations, counter%translations, &
            counter%trial_creations, counter%creations, &
            counter%trial_deletions, counter%deletions, &
            counter%trial_rotations, counter%deletions
        close(UNIT_MOVES)

    end subroutine WriteEnergyAndCount

    subroutine WriteLAMMPSData(box)

        implicit none

        ! Input parameters
        type(type_box), intent(inout) :: box

        ! Local variables
        integer :: i, j, k, atom_id, mol_id, atom_type
        integer :: cpt_bond, cpt_atom, cpt_angle, cpt_dihedral, cpt_improper
        integer :: UNIT_DATA = 19
        real(real64) :: charge
        real(real64), dimension(3) :: pos
        integer :: dim ! Integer for looping over dimensions

        ! Update bond number count
        cpt_bond = 0
        do i = 1, nb%type_residue
            do j = 1, box%num_residues(i)
                do k = 1, nb%bonds_per_residue(i)
                    cpt_bond = cpt_bond + 1
                end do
            end do
        end do
        box%num_bonds = cpt_bond

        ! Update angle number count
        cpt_angle = 0
        do i = 1, nb%type_residue
            do j = 1, box%num_residues(i)
                do k = 1, nb%angles_per_residue(i)
                    cpt_angle = cpt_angle + 1
                end do
            end do
        end do
        box%num_angles = cpt_angle

        ! Update dihedral number count
        cpt_dihedral = 0
        do i = 1, nb%type_residue
            do j = 1, box%num_residues(i)
                do k = 1, nb%dihedrals_per_residue(i)
                    cpt_dihedral = cpt_dihedral + 1
                end do
            end do
        end do
        box%num_dihedrals = cpt_dihedral

        ! Update improper number count
        cpt_improper = 0
        do i = 1, nb%type_residue
            do j = 1, box%num_residues(i)
                do k = 1, nb%impropers_per_residue(i)
                    cpt_improper = cpt_improper + 1
                end do
            end do
        end do
        box%num_impropers = cpt_improper

        ! Open file
        open(UNIT=UNIT_DATA, FILE=TRIM(output_path) // data_filename, STATUS='REPLACE', ACTION='write')

        write(UNIT_DATA, *) "! LAMMPS data file (atom_style full)"
        write(UNIT_DATA, *) box%num_atoms, " atoms"
        write(UNIT_DATA, *) box%num_atomtypes, " atom types"
        write(UNIT_DATA, *) box%num_bonds, " bonds"
        write(UNIT_DATA, *) box%num_bondtypes, " bond types"
        write(UNIT_DATA, *) box%num_angles, " angles"
        write(UNIT_DATA, *) box%num_angletypes, " angle types"
        write(UNIT_DATA, *) box%num_dihedrals, " dihedrals"
        write(UNIT_DATA, *) box%num_dihedraltypes, " dihedral types"
        write(UNIT_DATA, *) box%num_impropers, " impropers"
        write(UNIT_DATA, *) box%num_impropertypes, " improper types"
        write(UNIT_DATA, *)

        ! X bounds

        write(UNIT_DATA, '(2(F15.8,1X))', ADVANCE='NO') box%bounds(1,1), box%bounds(1,2)
        write(UNIT_DATA, '(A)') "xlo xhi"

        ! Y bounds
        write(UNIT_DATA, '(2(F15.8,1X))', ADVANCE='NO') box%bounds(2,1), box%bounds(2,2)
        write(UNIT_DATA, '(A)') "ylo yhi"

        ! Z bounds
        write(UNIT_DATA, '(2(F15.8,1X))', ADVANCE='NO') box%bounds(3,1), box%bounds(3,2)
        write(UNIT_DATA, '(A)') "zlo zhi"

        if (box%is_triclinic) then
            write(UNIT_DATA, '(3(F15.8,1X))') box%tilt(1), box%tilt(2), box%tilt(3)
            write(UNIT_DATA, '(A)') "xy xz yz"
        end if

        write(UNIT_DATA, *)

        ! Masses section (assumes atomic mass array `atom_masses`)
        write(UNIT_DATA, *) "Masses"
        write(UNIT_DATA, *)
        do atom_id = 1, primary%num_atomtypes
            write(UNIT_DATA, '(I5, 1X, F12.6)') atom_id, box%site_masses_vector(atom_id)
        end do

        write(UNIT_DATA, *)

        ! Atoms section
        write(UNIT_DATA, *) "Atoms"
        write(UNIT_DATA, *)

        atom_id = 0
        mol_id = 0
        do i = 1, nb%type_residue
            do j = 1, box%num_residues(i)
                mol_id = mol_id + 1
                do k = 1, nb%atom_in_residue(i)

                    atom_id = atom_id + 1
                    atom_type = box%atom_types(i,k)
                    charge = box%atom_charges(i,k)

                    do dim = 1, 3
                        pos(dim) = box%mol_com(dim, i, j) + box%site_offset(dim, i, j, k)
                    end do

                    ! Only wrap atom of inative species (i.e. leave active molecules continuous at the pbc)
                    if (input%is_active(i) == 0) then
                        call WrapIntoBox(pos, box)
                    end if

                    write(UNIT_DATA, '(I6,1X,I6,1X,I4,1X,F12.8,3(1X,F12.7))') &
                        atom_id, mol_id, atom_type, charge, pos(1), pos(2), pos(3)

                end do
            end do
        end do

        if (primary%num_bonds > 0 .or. reservoir%num_bonds > 0) then
            ! Atoms section
            write(UNIT_DATA, *)
            write(UNIT_DATA, *) "Bonds"
            write(UNIT_DATA, *)
            cpt_bond = 1
            cpt_atom = 0
            do i = 1, nb%type_residue
                do j = 1, box%num_residues(i)
                    do k = 1, nb%bonds_per_residue(i)
                        write(UNIT_DATA, *) cpt_bond, res%bond_type_2d(i, k, 1), &
                            cpt_atom + res%bond_type_2d(i, k, 2), &
                            cpt_atom + res%bond_type_2d(i, k, 3)
                        cpt_bond = cpt_bond + 1
                    end do
                    cpt_atom = cpt_atom + nb%atom_in_residue(i)
                end do
            end do
        end if

        if (primary%num_angles > 0 .or. reservoir%num_angles > 0) then
            ! Atoms section
            write(UNIT_DATA, *)
            write(UNIT_DATA, *) "Angles"
            write(UNIT_DATA, *)
            cpt_angle = 1
            cpt_atom = 0
            do i = 1, nb%type_residue
                do j = 1, box%num_residues(i)
                    do k = 1, nb%angles_per_residue(i)
                        write(UNIT_DATA, *) cpt_angle, res%angle_type_2d(i, k, 1), &
                            cpt_atom + res%angle_type_2d(i, k, 2), &
                            cpt_atom + res%angle_type_2d(i, k, 3), &
                            cpt_atom + res%angle_type_2d(i, k, 4)
                        cpt_angle = cpt_angle + 1
                    end do
                    cpt_atom = cpt_atom + nb%atom_in_residue(i)
                end do
            end do
        end if

        if (primary%num_dihedrals > 0 .or. reservoir%num_dihedrals > 0) then
            ! Atoms section
            write(UNIT_DATA, *)
            write(UNIT_DATA, *) "Dihedrals"
            write(UNIT_DATA, *)
            cpt_dihedral = 1
            cpt_atom = 0
            do i = 1, nb%type_residue
                do j = 1, box%num_residues(i)
                    do k = 1, nb%dihedrals_per_residue(i)
                        write(UNIT_DATA, *) cpt_dihedral, res%dihedral_type_2d(i, k, 1), &
                            cpt_atom + res%dihedral_type_2d(i, k, 2), &
                            cpt_atom + res%dihedral_type_2d(i, k, 3), &
                            cpt_atom + res%dihedral_type_2d(i, k, 4), &
                            cpt_atom + res%dihedral_type_2d(i, k, 5)
                        cpt_dihedral= cpt_dihedral + 1
                    end do
                    cpt_atom = cpt_atom + nb%atom_in_residue(i)
                end do
            end do
        end if

        if (primary%num_impropers > 0 .or. reservoir%num_impropers > 0) then
            ! Atoms section
            write(UNIT_DATA, *)
            write(UNIT_DATA, *) "Impropers"
            write(UNIT_DATA, *)
            cpt_improper = 1
            cpt_atom = 0
            do i = 1, nb%type_residue
                do j = 1, box%num_residues(i)
                    do k = 1, nb%impropers_per_residue(i)
                        write(UNIT_DATA, *) cpt_improper, res%improper_type_2d(i, k, 1), &
                            cpt_atom + res%improper_type_2d(i, k, 2), &
                            cpt_atom + res%improper_type_2d(i, k, 3), &
                            cpt_atom + res%improper_type_2d(i, k, 4), &
                            cpt_atom + res%improper_type_2d(i, k, 5)
                        cpt_improper= cpt_improper + 1
                    end do
                    cpt_atom = cpt_atom + nb%atom_in_residue(i)
                end do
            end do
        end if

        close(UNIT_DATA)

    end subroutine WriteLAMMPSData

    !------------------------------------------------------------------------------
    ! subroutine UpdateFiles
    ! Updates all relevant output files for the current Monte Carlo step.
    !------------------------------------------------------------------------------
    subroutine UpdateFiles(later_step)

        implicit none

        logical, intent(in) :: later_step  ! Flag to distinguish first vs. later steps

        ! Write LAMMPS trajectory
        call WriteLAMMPSTRJ(primary, "trajectory.lammpstrj", later_step)
        if (has_reservoir) call WriteLAMMPSTRJ(reservoir, "reservoir.lammpstrj", later_step)

        ! Write energies and move counts to data file
        call WriteEnergyAndCount()

        ! Write LAMMPS data file
        call WriteLAMMPSData(primary)

    end subroutine UpdateFiles

end module write_utils
