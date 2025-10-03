module data_parser

    use parameters
    use simulation_state
    use output_utils
    use geometry_utils
    use check_utils
    use readers_utils
    use, intrinsic :: iso_fortran_env, only: real64

    implicit none

    ! 1D atom information
    real(real64), allocatable :: atom_x_1d(:), atom_y_1d(:), atom_z_1d(:) ! Array of X Y Z coordinate of atoms
    real(real64), allocatable :: atom_charges_1d(:) ! Partial charges on sites
    real(real64), allocatable :: atom_masses_1d(:) ! Masses of atoms
    integer, allocatable :: atom_types_1d(:) ! Array of atom types
    integer, allocatable :: atom_ids_1d(:) ! Atom_id
    integer, allocatable :: atom_original_1d(:) ! Atom_id (from data file)
    character(10), allocatable :: atom_names_1d(:) ! Array of atom names

    integer, allocatable :: bond_ids_1d(:)
    integer, allocatable :: bond_types_1d(:)
    integer, allocatable :: bond_atoms_1d(:, :)

    integer, allocatable :: angle_ids_1d(:)
    integer, allocatable :: angle_types_1d(:)
    integer, allocatable :: angle_atoms_1d(:, :)

    integer, allocatable :: dihedral_ids_1d(:)
    integer, allocatable :: dihedral_types_1d(:)
    integer, allocatable :: dihedral_atoms_1d(:, :)

    integer, allocatable :: improper_ids_1d(:)
    integer, allocatable :: improper_types_1d(:)
    integer, allocatable :: improper_atoms_1d(:, :)

contains

    !-----------------------------------------------------------------------------
    ! Subroutine: ReadSystemData
    !
    ! Purpose:
    !   Reads the main system (primary) data file and optionally a reservoir file.
    !   Ensures consistency between primary and reservoir if applicable.
    !-----------------------------------------------------------------------------
    subroutine ReadSystemData()

        implicit none

        !-----------------------------
        ! Read primary system
        !-----------------------------
        call ReadLMPData(data_file, primary, .true.)

        !-----------------------------
        ! Read reservoir if provided
        !-----------------------------
        if (trim(res_file) == '') then
            has_reservoir = .false.
        else
            has_reservoir = .true.
            call ReadLMPData(res_file, reservoir, .false.)
            call CheckConsistency()   ! Ensure reservoir is consistent with primary
        end if

    end subroutine ReadSystemData

    subroutine ReadLMPData(data_file_name, box, is_primary)

        ! Input parameters
        character(len=*), intent(in) :: data_file_name   ! Path to the LAMMPS data file to read
        type(type_box), intent(inout) :: box
        logical, intent(in) :: is_primary

        ! Local variables
        integer :: INFILE = 100   ! Fortran unit number used to open the input file
        integer :: ios            ! I/O status flag (nonzero if file open/read fails)

        ! === Step 1: Open file ===
        open(UNIT=INFILE, FILE=data_file_name, STATUS='OLD', ACTION='read', IOSTAT=ios)

        if (ios /= 0) then
            call AbortRun("Error opening file: " // trim(data_file_name))
        end if

        ! === Step 2: Read header ===
        call ReadLMPHeaderInfo(INFILE, box)

        ! === Step 3: Allocate arrays ===
        allocate(atom_ids_1d(box%num_atoms))
        allocate(atom_original_1d(box%num_atoms))
        allocate(atom_types_1d(box%num_atoms))
        allocate(atom_charges_1d(box%num_atoms))
        allocate(atom_x_1d(box%num_atoms))
        allocate(atom_y_1d(box%num_atoms))
        allocate(atom_z_1d(box%num_atoms))
        allocate(atom_names_1d(box%num_atoms))
        allocate(bond_ids_1d(box%num_bonds))
        allocate(bond_types_1d(box%num_bonds))
        allocate(bond_atoms_1d(box%num_bonds,2))
        allocate(angle_ids_1d(box%num_angles))
        allocate(angle_types_1d(box%num_angles))
        allocate(angle_atoms_1d(box%num_angles,3))
        allocate(dihedral_ids_1d(box%num_dihedrals))
        allocate(dihedral_types_1d(box%num_dihedrals))
        allocate(dihedral_atoms_1d(box%num_dihedrals,4))
        allocate(improper_ids_1d(box%num_impropers))
        allocate(improper_types_1d(box%num_impropers))
        allocate(improper_atoms_1d(box%num_impropers,4))

        ! === Step 4: Read box dimensions ===
        rewind(INFILE)
        call ParseLAMMPSBox(INFILE, box)

        ! Detect box type and other informations
        call PrepareSimulationBox(box)

        ! === Step 5: Read masses ===
        rewind(INFILE)
        call ReadLAMMPSMasses(INFILE, box)

        ! === Step 6: Read atoms ===
        rewind(INFILE)
        call ReadLAMMPSAtoms(INFILE, data_file_name, box)
        call AssignAtomNames(box)

        ! === Step 6: Read bonds ===
        rewind(INFILE)
        call ReadLAMMPSBonds(INFILE, data_file_name, box)

        ! === Step 6: Read angles ===
        rewind(INFILE)
        call ReadLAMMPSAngles(INFILE, data_file_name, box)

        ! === Step 7: Read dihedrals ===
        rewind(INFILE)
        call ReadLAMMPSDihedrals(INFILE, data_file_name, box)

        ! === Step 7: Read improper ===
        rewind(INFILE)
        call ReadLAMMPSImpropers(INFILE, data_file_name, box)

        close(INFILE)

        ! === Step 7: Post-processing ===
        call SortAtomsByOriginalID(box)
        call DetectMolecules(box)
        call RepairActiveMolecules(box)
        call TransformCoordinate(box)
        call CheckMolecule(box)

        ! === Step 7: Detect bond, angle, dihedral, improper per residue ===
        call DetectBondPerResidue(box)
        call DetectAnglePerResidue(box)
        call DetectDihedralPerResidue(box)
        call DetectImproperPerResidue(box)

        ! === Step 8: Deallocate arrays ===
        if (allocated(atom_ids_1d)) deallocate(atom_ids_1d)
        if (allocated(atom_original_1d)) deallocate(atom_original_1d)
        if (allocated(atom_types_1d)) deallocate(atom_types_1d)
        if (allocated(atom_charges_1d)) deallocate(atom_charges_1d)
        if (allocated(atom_x_1d)) deallocate(atom_x_1d)
        if (allocated(atom_y_1d)) deallocate(atom_y_1d)
        if (allocated(atom_z_1d)) deallocate(atom_z_1d)
        if (allocated(atom_names_1d)) deallocate(atom_names_1d)
        if (allocated(bond_ids_1d)) deallocate(bond_ids_1d)
        if (allocated(bond_types_1d)) deallocate(bond_types_1d)
        if (allocated(bond_atoms_1d)) deallocate(bond_atoms_1d)
        if (allocated(angle_ids_1d)) deallocate(angle_ids_1d)
        if (allocated(angle_types_1d)) deallocate(angle_types_1d)
        if (allocated(angle_atoms_1d)) deallocate(angle_atoms_1d)
        if (allocated(dihedral_ids_1d)) deallocate(dihedral_ids_1d)
        if (allocated(dihedral_types_1d)) deallocate(dihedral_types_1d)
        if (allocated(dihedral_atoms_1d)) deallocate(dihedral_atoms_1d)
        if (allocated(improper_ids_1d)) deallocate(improper_ids_1d)
        if (allocated(improper_types_1d)) deallocate(improper_types_1d)
        if (allocated(improper_atoms_1d)) deallocate(improper_atoms_1d)

        ! === Setep 9: Log data ===
        call LogData(data_file_name, box, is_primary)
        call LogConnectivity(box)

    end subroutine ReadLMPData

    !-----------------------------------------------------------------------------
    ! subroutine ReadLAMMPSMasses(INFILE, box)
    !
    ! Reads atomic masses from the "Masses" section of a LAMMPS data file and maps them
    ! into a 2D mass array for each atom type and residue type.
    !-----------------------------------------------------------------------------
    subroutine ReadLAMMPSMasses(INFILE, box)
        implicit none

        ! Input parameters
        type(type_box), intent(inout) :: box
        integer, intent(in) :: INFILE

        ! Local variables
        character(len=256) :: line, trimmed_line
        integer :: ios, i, j, k, tmp_int, mass_found
        real(real64) :: tmp_flt

        ! Initialize
        mass_found = 0

        allocate(box%site_masses_vector(box%num_atomtypes))
        box%site_masses_vector = 0.0d0

        ! Rewind file to start (optional, remove if file pointer is already positioned)
        rewind(INFILE)

        ! Look for "Masses" section
        do
            read(INFILE, '(A)', IOSTAT=ios) line
            if (ios /= 0) then
                exit
            end if

            ! Trim leading and trailing spaces
            trimmed_line = adjustl(trim(line))

            ! Skip empty lines or comments
            if (len_trim(trimmed_line) == 0 .or. trimmed_line(1:1) == '!') then
                cycle
            end if

            ! Check for "Masses" keyword
            if (trim(trimmed_line) == "Masses") then
                ! Skip next line (usually blank or a comment)
                read(INFILE, '(A)', IOSTAT=ios) line
                if (ios /= 0) then
                    exit
                end if
                trimmed_line = adjustl(trim(line))
                if (len_trim(trimmed_line) > 0 .and. trimmed_line(1:1) /= '!') then
                end if

                ! Read num_atoms(SYSTEM) lines for masses
                do i = 1, box%num_atoms
                    read(INFILE, '(A)', IOSTAT=ios) line
                    if (ios /= 0) then
                        exit
                    end if

                    trimmed_line = adjustl(trim(line))
                    ! Skip empty or comment lines
                    if (len_trim(trimmed_line) == 0 .or. trimmed_line(1:1) == '!') then
                        cycle
                    end if

                    ! Parse species ID and mass
                    read(trimmed_line, *, IOSTAT=ios) tmp_int, tmp_flt
                    if (ios /= 0) then
                        exit
                    end if

                    if (tmp_int < 1 .or. tmp_int > box%num_atoms) then
                        exit
                    end if

                    box%site_masses_vector(tmp_int) = tmp_flt
                    mass_found = mass_found + 1
                end do

                ! Map site masses into 2D mass array (assuming nb%type_residue and nb%types_per_residue are defined)
                i = 1
                do j = 1, nb%type_residue
                    do k = 1, nb%types_per_residue(j)
                        if ((i <= box%num_atoms) .and. (i <= box%num_atomtypes)) then
                            box%atom_masses(j, k) = box%site_masses_vector(i)
                            i = i + 1
                        else
                            exit
                        end if
                    end do
                end do

                exit  ! Exit after processing Masses section
            end if
        end do

        ! Trigger error if no masses found or mismatch with declared atom types
        if (mass_found == 0) then
            call AbortRun("No masses found in data file", 12)
        else if (mass_found /= box%num_atomtypes) then
            call AbortRun("Number of masses found in data file differs from declared atom types", 13)
        end if

    end subroutine ReadLAMMPSMasses

    subroutine AssignAtomNames(box)

        implicit none

        type(type_box), intent(inout) :: box

        integer :: i, j, k

        ! Loop over all atoms
        do k = 1, box%num_atoms
            atom_names_1d(k) = '' ! Initialize to empty string
            do i = 1, nb%type_residue
                do j = 1, nb%types_per_residue(i)
                    if (res%types_2d(i, j) == atom_types_1d(k)) then
                        atom_names_1d(k) = res%names_2d(i, j)
                        exit  ! Exit inner loop once match is found
                    end if
                end do
                if (len_trim(atom_names_1d(k)) > 0) exit  ! Exit outer loop if found
            end do
            if (len_trim(atom_names_1d(k)) == 0) then
                atom_names_1d(k) = 'Unknown'  ! Optional default value
            end if
        end do

    end subroutine AssignAtomNames

    subroutine DetectBondPerResidue(box)

        implicit none

        type(type_box), intent(inout) :: box

        integer :: i, j
        integer :: id1, id2
        integer :: k, l
        integer :: nb_bond
        character(len=100) :: formatted_msg

        do i = 1, nb%type_residue
            nb_bond = 0

            ! Loop directly over the bonds
            do j = 1, box%num_bonds

                ! Ids in bonds
                id1 = bond_atoms_1d(j,1)
                id2 = bond_atoms_1d(j,2)

                ! Check if both atoms belong to residue i
                if (isInResidue(box, i, id1) .and. isInResidue(box, i, id2)) then

                    ! Get atom local index in the residue
                    k = atomIndexInResidue(box, i, id1)
                    l = atomIndexInResidue(box, i, id2)
                    ! Store bond info
                    nb_bond = nb_bond + 1

                    if (nb_bond > NB_MAX_BOND) then
                        write(formatted_msg, '(A, I0)') &
                            "The number of bonds exceeds the maximum allowed = ", NB_MAX_BOND
                        call AbortRun(trim(formatted_msg) // &
                            " Please increase 'NB_MAX_BOND' in src/parameters.f90 and recompile.", 11)
                    end if

                    res%bond_type_2d(i, nb_bond, 1) = bond_types_1d(j)

                    if (l < k) then
                        res%bond_type_2d(i, nb_bond, 2) = l
                        res%bond_type_2d(i, nb_bond, 3) = k
                    else
                        res%bond_type_2d(i, nb_bond, 2) = k
                        res%bond_type_2d(i, nb_bond, 3) = l
                    end if
                end if

            end do

            nb%bonds_per_residue(i) = nb_bond

        end do

    end subroutine DetectBondPerResidue

    subroutine DetectAnglePerResidue(box)

        implicit none

        type(type_box), intent(inout) :: box

        integer :: i, j, k, l, m
        integer :: id1, id2, id3
        integer :: nb_angle
        character(len=100) :: formatted_msg

        do i = 1, nb%type_residue
            nb_angle = 0

            ! Loop directly over the angles
            do j = 1, box%num_angles
                id1 = angle_atoms_1d(j,1)
                id2 = angle_atoms_1d(j,2)
                id3 = angle_atoms_1d(j,3)

                ! Check if all 3 atoms belong to residue i
                if (isInResidue(box, i, id1) .and. isInResidue(box, i, id2) .and. isInResidue(box, i, id3)) then

                    ! Get local indices inside the residue
                    k = atomIndexInResidue(box, i, id1)
                    l = atomIndexInResidue(box, i, id2)
                    m = atomIndexInResidue(box, i, id3)

                    ! Store angle info
                    nb_angle = nb_angle + 1

                    if (nb_angle > NB_MAX_ANGLE) then
                        write(formatted_msg, '(A, I0)') &
                            "The number of angles exceeds the maximum allowed = ", NB_MAX_ANGLE
                        call AbortRun(trim(formatted_msg) // &
                            " Please increase 'NB_MAX_ANGLE' in src/parameters.f90 and recompile.", 11)
                    end if

                    res%angle_type_2d(i, nb_angle, 1) = angle_types_1d(j)
                    if (k < m) then
                        res%angle_type_2d(i, nb_angle, 2) = k
                        res%angle_type_2d(i, nb_angle, 3) = l
                        res%angle_type_2d(i, nb_angle, 4) = m
                    else
                        res%angle_type_2d(i, nb_angle, 2) = m
                        res%angle_type_2d(i, nb_angle, 3) = l
                        res%angle_type_2d(i, nb_angle, 4) = k
                    end if
                end if
            end do
            nb%angles_per_residue(i) = nb_angle
        end do

    end subroutine DetectAnglePerResidue

    subroutine DetectDihedralPerResidue(box)

        implicit none

        type(type_box), intent(inout) :: box

        integer :: i, j, k, l, m, n
        integer :: id1, id2, id3, id4
        integer :: nb_dihedral
        character(len=100) :: formatted_msg

        do i = 1, nb%type_residue
            nb_dihedral = 0

            ! Loop directly over the dihedrals
            do j = 1, box%num_dihedrals
                id1 = dihedral_atoms_1d(j,1)
                id2 = dihedral_atoms_1d(j,2)
                id3 = dihedral_atoms_1d(j,3)
                id4 = dihedral_atoms_1d(j,4)

                ! Check if all 4 atoms belong to residue i
                if (isInResidue(box, i, id1) .and. isInResidue(box, i, id2) .and. &
                    isInResidue(box, i, id3) .and. isInResidue(box, i, id4)) then

                    ! Get local indices inside the residue
                    k = atomIndexInResidue(box, i, id1)
                    l = atomIndexInResidue(box, i, id2)
                    m = atomIndexInResidue(box, i, id3)
                    n = atomIndexInResidue(box, i, id4)

                    ! Store dihedral info
                    nb_dihedral = nb_dihedral + 1

                    if (nb_dihedral > NB_MAX_DIHEDRAL) then
                        write(formatted_msg, '(A, I0)') &
                            "The number of dihedrals exceeds the maximum allowed = ", NB_MAX_DIHEDRAL
                        call AbortRun(trim(formatted_msg) // &
                            " Please increase 'NB_MAX_DIHEDRAL' in src/parameters.f90 and recompile.", 11)
                    end if

                    res%dihedral_type_2d(i, nb_dihedral, 1) = dihedral_types_1d(j)
                    if (k < n) then
                        res%dihedral_type_2d(i, nb_dihedral, 2) = k
                        res%dihedral_type_2d(i, nb_dihedral, 3) = l
                        res%dihedral_type_2d(i, nb_dihedral, 4) = m
                        res%dihedral_type_2d(i, nb_dihedral, 5) = n
                    else
                        res%dihedral_type_2d(i, nb_dihedral, 2) = n
                        res%dihedral_type_2d(i, nb_dihedral, 3) = m
                        res%dihedral_type_2d(i, nb_dihedral, 4) = l
                        res%dihedral_type_2d(i, nb_dihedral, 5) = k
                    end if
                end if
            end do
            nb%dihedrals_per_residue(i) = nb_dihedral
        end do

    end subroutine DetectDihedralPerResidue

    subroutine DetectImproperPerResidue(box)

        implicit none

        type(type_box), intent(inout) :: box

        integer :: i, j, k, l, m, n
        integer :: id1, id2, id3, id4
        integer :: nb_improper
        character(len=100) :: formatted_msg

        do i = 1, nb%type_residue
            nb_improper = 0

            ! Loop directly over the impropers
            do j = 1, box%num_impropers
                id1 = improper_atoms_1d(j,1)
                id2 = improper_atoms_1d(j,2)
                id3 = improper_atoms_1d(j,3)
                id4 = improper_atoms_1d(j,4)

                ! Check if all 4 atoms belong to residue i
                if (isInResidue(box, i, id1) .and. isInResidue(box, i, id2) .and. &
                    isInResidue(box, i, id3) .and. isInResidue(box, i, id4)) then

                    ! Get local indices inside the residue
                    k = atomIndexInResidue(box, i, id1)
                    l = atomIndexInResidue(box, i, id2)
                    m = atomIndexInResidue(box, i, id3)
                    n = atomIndexInResidue(box, i, id4)

                    ! Store improper info
                    nb_improper = nb_improper + 1

                    if (nb_improper > NB_MAX_IMPROPER) then
                        write(formatted_msg, '(A, I0)') &
                            "The number of impropers exceeds the maximum allowed = ", NB_MAX_IMPROPER
                        call AbortRun(trim(formatted_msg) // &
                            " Please increase 'NB_MAX_IMPROPER' in src/parameters.f90 and recompile.", 11)
                    end if

                    res%improper_type_2d(i, nb_improper, 1) = improper_types_1d(j)
                    if (k < n) then
                        res%improper_type_2d(i, nb_improper, 2) = k
                        res%improper_type_2d(i, nb_improper, 3) = l
                        res%improper_type_2d(i, nb_improper, 4) = m
                        res%improper_type_2d(i, nb_improper, 5) = n
                    else
                        res%improper_type_2d(i, nb_improper, 2) = n
                        res%improper_type_2d(i, nb_improper, 3) = m
                        res%improper_type_2d(i, nb_improper, 4) = l
                        res%improper_type_2d(i, nb_improper, 5) = k
                    end if
                end if
            end do
            nb%impropers_per_residue(i) = nb_improper
        end do

    end subroutine DetectImproperPerResidue

    subroutine ReadLAMMPSAtoms(INFILE, data_file_name, box)
        implicit none

        ! Input
        integer, intent(in) :: INFILE
        character(len=*), intent(in) :: data_file_name
        type(type_box), intent(inout) :: box

        ! Local
        character(len=200) :: line, formatted_msg
        integer :: ios, k
        integer :: atom_found
        integer :: tmp1_int, tmp2_int, tmp3_int
        real(real64) :: tmp1_flt, tmp2_flt, tmp3_flt, tmp4_flt
        integer, parameter :: ERROR_NO_ATOMS = 13
        integer, parameter :: ERROR_READ = 14
        integer, parameter :: ERROR_PARSE = 15
        integer, parameter :: ERROR_ATOM_TYPE = 16
        integer, parameter :: ERROR_BOUNDS = 17

        atom_found = 0

        ! Read until "Atoms" section or end of file
        do
            read(INFILE, '(A)', IOSTAT=ios) line

            if (ios < 0) then
                write(formatted_msg, '(A, A)') "No atoms found in data file: ", trim(data_file_name)
                call AbortRun(trim(formatted_msg), ERROR_NO_ATOMS)
            end if

            if (ios > 0) then
                write(formatted_msg, '(A, A)') "I/O error reading data file: ", trim(data_file_name)
                call AbortRun(trim(formatted_msg), ERROR_READ)
            end if

            ! Check for "Atoms" or "Atoms # full" (case-insensitive)
            if (index(ADJUSTL(line), 'Atoms') == 1) then
                exit
            end if
        end do

        ! Read and discard a blank line, if present
        read(INFILE, '(A)', IOSTAT=ios) line
        if (ios < 0) then
            write(formatted_msg, '(A, A)') "No atoms found in data file: ", trim(data_file_name)
            call AbortRun(trim(formatted_msg), ERROR_READ)
        end if

        if (TRIM(ADJUSTL(line)) /= "") then
            ! If not blank, rewind to re-read as atom data
            backspace(INFILE, IOSTAT=ios)
            if (ios /= 0) then
                write(formatted_msg, '(A, A)') "Failed to rewind file after non-blank line in: ", trim(data_file_name)
                call AbortRun(trim(formatted_msg), ERROR_READ)
            end if
        end if

        ! Loop over all atoms
        do k = 1, box%num_atoms
            read(INFILE, '(A)', IOSTAT=ios) line

            if (ios < 0) then
                write(formatted_msg, '(A, I0, A)') &
                    "Unexpected end of file at atom line ", k, " in: " // trim(data_file_name)
                call AbortRun(trim(formatted_msg), ERROR_READ)
            end if

            if (ios > 0) then
                write(formatted_msg, '(A, I0, A)') &
                    "I/O error reading atom line ", k, " in: " // trim(data_file_name)
                call AbortRun(trim(formatted_msg), ERROR_READ)
            end if

            ! Parse atom info: id, molecule-ID, type, q, x, y, z
            read(line, *, IOSTAT=ios) tmp1_int, tmp2_int, tmp3_int, &
                                    tmp1_flt, tmp2_flt, tmp3_flt, tmp4_flt

            if (ios /= 0) then
                write(formatted_msg, '(A, A)') "Failed to parse atom line: '" &
                    // TRIM(line) // "' in: " // trim(data_file_name)
                call AbortRun(trim(formatted_msg), ERROR_PARSE)
            end if

            ! Validate atom type
            if (tmp3_int < 1 .or. tmp3_int > box%num_atomtypes) then
                write(formatted_msg, '(A, I0, A, I0, A)') "Invalid atom type ", tmp3_int, &
                                                        " (max allowed: ", box%num_atomtypes, &
                                                        ") in: " // trim(data_file_name)
                call AbortRun(trim(formatted_msg), ERROR_ATOM_TYPE)
            end if

            ! Validate array bounds
            if (k < 1 .or. k > size(atom_ids_1d)) then
                write(formatted_msg, '(A, I0, A, I0, A)') "Invalid array index k=", k, &
                                                        " (array size: ", size(atom_ids_1d), &
                                                        ") in: " // trim(data_file_name)
                call AbortRun(trim(formatted_msg), ERROR_BOUNDS)
            end if

            ! Store atom properties
            atom_ids_1d(k) = k
            atom_original_1d(k) = tmp1_int
            atom_types_1d(k) = tmp3_int
            atom_charges_1d(k) = tmp1_flt
            atom_x_1d(k) = tmp2_flt
            atom_y_1d(k) = tmp3_flt
            atom_z_1d(k) = tmp4_flt
            atom_found = atom_found + 1
        end do

        ! Check if the number of atoms found matches the expected number
        if (atom_found /= box%num_atoms) then
            write(formatted_msg, '(A, I0, A, I0, A)') "Found ", atom_found, &
                                                    " atoms, expected ", box%num_atoms, &
                                                    " in: " // trim(data_file_name)
            call AbortRun(trim(formatted_msg), ERROR_NO_ATOMS)
        end if

    end subroutine ReadLAMMPSAtoms

    subroutine ReadLAMMPSBonds(INFILE, data_file_name, box)
        implicit none

        ! Input
        integer, intent(in) :: INFILE
        character(len=*), intent(in) :: data_file_name
        type(type_box), intent(inout) :: box

        ! Local
        character(len=200) :: line, formatted_msg
        integer :: ios, k
        integer :: bond_found
        integer :: tmp1_int, tmp2_int, tmp3_int, tmp4_int
        integer, parameter :: ERROR_READ = 24
        integer, parameter :: ERROR_PARSE = 25

        bond_found = 0

        ! If no bonds are expected, skip the whole routine

        if (box%num_bonds == 0) then
            call InfoMessage("No bonds expected in data file: " // trim(data_file_name))
            return
        end if

        ! Read until "Bonds" section or end of file
        do
            read(INFILE, '(A)', IOSTAT=ios) line

            if (ios < 0) then
                write(formatted_msg, '(A, A)') "No bonds found in data file: ", trim(data_file_name)
                call InfoMessage(trim(formatted_msg))
            end if

            if (ios > 0) then
                write(formatted_msg, '(A, A)') "I/O error reading data file: ", trim(data_file_name)
                call AbortRun(trim(formatted_msg), ERROR_READ)
            end if

            ! Check for "Bonds" or "Bonds # full" (case-insensitive)
            if (index(ADJUSTL(line), 'Bonds') == 1) then
                exit
            end if

        end do

        ! Read and discard a blank line, if present
        read(INFILE, '(A)', IOSTAT=ios) line

        if (ios < 0) then
            write(formatted_msg, '(A, A)') "No bonds found in data file: ", trim(data_file_name)
            call AbortRun(trim(formatted_msg), ERROR_READ)
        end if

        if (TRIM(ADJUSTL(line)) /= "") then
            ! If not blank, rewind to re-read as bond data
            backspace(INFILE, IOSTAT=ios)
            if (ios /= 0) then
                write(formatted_msg, '(A, A)') "Failed to rewind file after non-blank line in: " // trim(data_file_name)
                call AbortRun(trim(formatted_msg), ERROR_READ)
            end if
        end if

        ! Loop over all bonds
        do k = 1, box%num_bonds
            read(INFILE, '(A)', IOSTAT=ios) line

            if (ios < 0) then
                write(formatted_msg, '(A, I0, A)') "Unexpected end of file at bond line ", k, " in: " // trim(data_file_name)
                call AbortRun(trim(formatted_msg), ERROR_READ)
            end if

            if (ios > 0) then
                write(formatted_msg, '(A, I0, A)') "I/O error reading bond line ", k, " in: " // trim(data_file_name)
                call AbortRun(trim(formatted_msg), ERROR_READ)
            end if

            ! Parse bond info: id, type, atom1, atom2
            read(line, *, IOSTAT=ios) tmp1_int, tmp2_int, tmp3_int, tmp4_int
            if (ios /= 0) then
                write(formatted_msg, '(A, A)') "Failed to parse bond line: '" // TRIM(line) // "' in: " // trim(data_file_name)
                call AbortRun(trim(formatted_msg), ERROR_PARSE)
            end if

            ! Store bond properties
            bond_ids_1d(k) = tmp1_int
            bond_types_1d(k) = tmp2_int
            bond_atoms_1d(k,1) = tmp3_int
            bond_atoms_1d(k,2) = tmp4_int
            bond_found = bond_found + 1
        end do

    end subroutine ReadLAMMPSBonds

    subroutine ReadLAMMPSAngles(INFILE, data_file_name, box)
        implicit none

        ! Input
        integer, intent(in) :: INFILE
        character(len=*), intent(in) :: data_file_name
        type(type_box), intent(inout) :: box

        ! Local
        character(len=200) :: line, formatted_msg
        integer :: ios, k
        integer :: angle_found
        integer :: tmp1_int, tmp2_int, tmp3_int, tmp4_int, tmp5_int
        integer, parameter :: ERROR_READ = 24
        integer, parameter :: ERROR_PARSE = 25

        angle_found = 0

        ! If no angles are expected, skip the whole routine
        if (box%num_angles == 0) then
            write(formatted_msg, '(A, A)') "No angles expected in data file: ", trim(data_file_name)
            call InfoMessage(trim(formatted_msg))
            return
        end if

        ! Read until "Angles" section or end of file
        do
            read(INFILE, '(A)', IOSTAT=ios) line

            if (ios < 0) then
                write(formatted_msg, '(A, A)') "No angles found in data file: ", trim(data_file_name)
                call InfoMessage(trim(formatted_msg))
            end if

            if (ios > 0) then
                write(formatted_msg, '(A, A)') "I/O error reading data file: ", trim(data_file_name)
                call AbortRun(trim(formatted_msg), ERROR_READ)
            end if

            ! Check for "Angles" or "Angles # full" (case-insensitive)
            if (index(ADJUSTL(line), 'Angles') == 1) then
                exit
            end if
        end do

        ! Read and discard a blank line, if present
        read(INFILE, '(A)', IOSTAT=ios) line

        if (ios < 0) then
            write(formatted_msg, '(A, A)') "No angles found in data file: ", trim(data_file_name)
            call AbortRun(trim(formatted_msg), ERROR_READ)
        end if

        if (TRIM(ADJUSTL(line)) /= "") then
            ! If not blank, rewind to re-read as angle data
            backspace(INFILE, IOSTAT=ios)
            if (ios /= 0) then
                write(formatted_msg, '(A, A)') "Failed to rewind file after non-blank line in: ", trim(data_file_name)
                call AbortRun(trim(formatted_msg), ERROR_READ)
            end if
        end if

        ! Loop over all angles
        do k = 1, box%num_angles
            read(INFILE, '(A)', IOSTAT=ios) line

            if (ios < 0) then
                write(formatted_msg, '(A, I0, A)') "Unexpected end of file at angle line ", k, " in: " // trim(data_file_name)
                call AbortRun(trim(formatted_msg), ERROR_READ)
            end if

            if (ios > 0) then
                write(formatted_msg, '(A, I0, A)') "I/O error reading angle line ", k, " in: " // trim(data_file_name)
                call AbortRun(trim(formatted_msg), ERROR_READ)
            end if

            ! Parse angle info: id, type, atom1, atom2
            read(line, *, IOSTAT=ios) tmp1_int, tmp2_int, tmp3_int, tmp4_int, tmp5_int
            if (ios /= 0) then
                write(formatted_msg, '(A, A)') "Failed to parse angle line: '" // TRIM(line) // "' in: " // trim(data_file_name)
                call AbortRun(trim(formatted_msg), ERROR_PARSE)
            end if

            ! Store bond properties
            angle_ids_1d(k) = tmp1_int
            angle_types_1d(k) = tmp2_int
            angle_atoms_1d(k,1) = tmp3_int
            angle_atoms_1d(k,2) = tmp4_int
            angle_atoms_1d(k,3) = tmp5_int
            angle_found = angle_found + 1
        end do

    end subroutine ReadLAMMPSAngles

    subroutine ReadLAMMPSDihedrals(INFILE, data_file_name, box)
        implicit none

        ! Input
        integer, intent(in) :: INFILE
        character(len=*), intent(in) :: data_file_name
        type(type_box), intent(inout) :: box

        ! Local
        character(len=200) :: line, formatted_msg
        integer :: ios, k
        integer :: dihedral_found
        integer :: tmp1_int, tmp2_int, tmp3_int, tmp4_int, tmp5_int, tmp6_int
        integer, parameter :: ERROR_READ = 24
        integer, parameter :: ERROR_PARSE = 25

        dihedral_found = 0

        ! If no dihedral are expected, skip the whole routine
        if (box%num_dihedrals == 0) then
            write(formatted_msg, '(A, A)') "No dihedrals expected in data file: ", trim(data_file_name)
            call InfoMessage(trim(formatted_msg))
            return
        end if

        ! Read until "Dihedrals" section or end of file
        do
            read(INFILE, '(A)', IOSTAT=ios) line

            if (ios < 0) then
                write(formatted_msg, '(A, A)') "No dihedrals found in data file: ", trim(data_file_name)
                call InfoMessage(trim(formatted_msg))
            end if

            if (ios > 0) then
                write(formatted_msg, '(A, A)') "I/O error reading data file: ", trim(data_file_name)
                call AbortRun(trim(formatted_msg), ERROR_READ)
            end if

            ! Check for "Dihedrals" or "Dihedrals # full" (case-insensitive)
            if (index(ADJUSTL(line), 'Dihedrals') == 1) then
                exit
            end if
        end do

        ! Read and discard a blank line, if present
        read(INFILE, '(A)', IOSTAT=ios) line

        if (ios < 0) then
            write(formatted_msg, '(A, A)') "No dihedrals found in data file: ", trim(data_file_name)
            call AbortRun(trim(formatted_msg), ERROR_READ)
        end if

        if (TRIM(ADJUSTL(line)) /= "") then
            ! If not blank, rewind to re-read as bond data
            backspace(INFILE, IOSTAT=ios)
            if (ios /= 0) then
                write(formatted_msg, '(A, A)') "Failed to rewind file after non-blank line in: ", trim(data_file_name)
                call AbortRun(trim(formatted_msg), ERROR_READ)
            end if
        end if

        ! Loop over all dihedrals
        do k = 1, box%num_dihedrals
            read(INFILE, '(A)', IOSTAT=ios) line

            if (ios < 0) then
                write(formatted_msg, '(A, I0, A)') "Unexpected end of file at dihedrals line ", &
                    k, " in: " // trim(data_file_name)
                call AbortRun(trim(formatted_msg), ERROR_READ)
            end if

            if (ios > 0) then
                write(formatted_msg, '(A, I0, A)') "I/O error reading dihedrals line ", &
                    k, " in: " // trim(data_file_name)
                call AbortRun(trim(formatted_msg), ERROR_READ)
            end if

            ! Parse dihedrals info: id, type, atom1, atom2
            read(line, *, IOSTAT=ios) tmp1_int, tmp2_int, tmp3_int, tmp4_int, tmp5_int, tmp6_int
            if (ios /= 0) then
                write(formatted_msg, '(A, A)') "Failed to parse dihedrals line: '" // &
                    TRIM(line) // "' in: " // trim(data_file_name)
                call AbortRun(trim(formatted_msg), ERROR_PARSE)
            end if

            ! Store bond properties
            dihedral_ids_1d(k) = tmp1_int
            dihedral_types_1d(k) = tmp2_int
            dihedral_atoms_1d(k,1) = tmp3_int
            dihedral_atoms_1d(k,2) = tmp4_int
            dihedral_atoms_1d(k,3) = tmp5_int
            dihedral_atoms_1d(k,4) = tmp6_int
            dihedral_found = dihedral_found + 1
        end do

    end subroutine ReadLAMMPSDihedrals

    subroutine ReadLAMMPSImpropers(INFILE, data_file_name, box)
        implicit none

        ! Input
        integer, intent(in) :: INFILE
        character(len=*), intent(in) :: data_file_name
        type(type_box), intent(inout) :: box

        ! Local
        character(len=200) :: line, formatted_msg
        integer :: ios, k
        integer :: improper_found
        integer :: tmp1_int, tmp2_int, tmp3_int, tmp4_int, tmp5_int, tmp6_int
        integer, parameter :: ERROR_READ = 24
        integer, parameter :: ERROR_PARSE = 25

        improper_found = 0

        ! If no improper are expected, skip the whole routine
        if (box%num_impropers == 0) then
            write(formatted_msg, '(A, A)') "No impropers expected in data file: ", trim(data_file_name)
            call InfoMessage(trim(formatted_msg))
            return
        end if

        ! Read until "Impropers" section or end of file
        do
            read(INFILE, '(A)', IOSTAT=ios) line
            if (ios < 0) then
                write(formatted_msg, '(A, A)') "No impropers found in data file: ", trim(data_file_name)
                call InfoMessage(trim(formatted_msg))
            end if
            if (ios > 0) then
                write(formatted_msg, '(A, A)') "I/O error reading data file: ", trim(data_file_name)
                call AbortRun(trim(formatted_msg), ERROR_READ)
            end if

            ! Check for "Impropers" or "Impropers # full" (case-insensitive)
            if (index(ADJUSTL(line), 'Impropers') == 1) then
                exit
            end if
        end do

        ! Read and discard a blank line, if present
        read(INFILE, '(A)', IOSTAT=ios) line

        if (ios < 0) then
            write(formatted_msg, '(A, A)') "No impropers found in data file: ", trim(data_file_name)
            call AbortRun(trim(formatted_msg), ERROR_READ)
        end if

        if (TRIM(ADJUSTL(line)) /= "") then
            ! If not blank, rewind to re-read as bond data
            backspace(INFILE, IOSTAT=ios)
            if (ios /= 0) then
                write(formatted_msg, '(A, A)') "Failed to rewind file after non-blank line in: ", trim(data_file_name)
                call AbortRun(trim(formatted_msg), ERROR_READ)
            end if
        end if

        ! Loop over all impropers
        do k = 1, box%num_impropers
            read(INFILE, '(A)', IOSTAT=ios) line

            if (ios < 0) then
                write(formatted_msg, '(A, I0, A)') &
                    "Unexpected end of file at impropers line ", k, " in: " // trim(data_file_name)
                call AbortRun(trim(formatted_msg), ERROR_READ)
            end if

            if (ios > 0) then
                write(formatted_msg, '(A, I0, A)') "I/O error reading impropers line ", &
                    k, " in: " // trim(data_file_name)
                call AbortRun(trim(formatted_msg), ERROR_READ)
            end if

            ! Parse impropers info: id, type, atom1, atom2
            read(line, *, IOSTAT=ios) tmp1_int, tmp2_int, tmp3_int, tmp4_int, tmp5_int, tmp6_int
            if (ios /= 0) then
                write(formatted_msg, '(A, A)') "Failed to parse impropers line: '" &
                    // TRIM(line) // "' in: " // trim(data_file_name)
                call AbortRun(trim(formatted_msg), ERROR_PARSE)
            end if

            ! Store bond properties
            improper_ids_1d(k) = tmp1_int
            improper_types_1d(k) = tmp2_int
            improper_atoms_1d(k,1) = tmp3_int
            improper_atoms_1d(k,2) = tmp4_int
            improper_atoms_1d(k,3) = tmp5_int
            improper_atoms_1d(k,4) = tmp6_int
            improper_found = improper_found + 1
        end do

    end subroutine ReadLAMMPSImpropers

    subroutine SortAtomsByOriginalID(box)

        implicit none

        ! Input parameters
        type(type_box), intent(inout) :: box

        ! Local variables
        integer :: i, j, temp_index, tmp_orig
        integer, allocatable :: sort_index(:)
        character(10), allocatable :: tmp_atom_names_1d(:)
        real(real64), allocatable :: tmp_atom_charges_1d(:)
        real(real64), allocatable :: tmp_atom_x_1d(:)
        real(real64), allocatable :: tmp_atom_y_1d(:)
        real(real64), allocatable :: tmp_atom_z_1d(:)
        integer, allocatable :: tmp_atom_types_1d(:)
        integer, allocatable :: tmp_atom_original_1d(:)

        ! Allocate sorting index array
        allocate(sort_index(box%num_atoms))
        do i = 1, box%num_atoms
            sort_index(i) = i
        end do

        ! Insertion sort based on atom_original_1d
        do i = 2, box%num_atoms
            temp_index = sort_index(i)
            tmp_orig = atom_original_1d(temp_index)
            j = i - 1

            do while (j >= 1)
                if (atom_original_1d(sort_index(j)) > tmp_orig) then
                    sort_index(j + 1) = sort_index(j)
                    j = j - 1
                else
                    exit
                end if
            end do

            sort_index(j + 1) = temp_index
        end do

        ! Allocate temporary arrays
        allocate(tmp_atom_names_1d(box%num_atoms))
        allocate(tmp_atom_charges_1d(box%num_atoms))
        allocate(tmp_atom_x_1d(box%num_atoms))
        allocate(tmp_atom_y_1d(box%num_atoms))
        allocate(tmp_atom_z_1d(box%num_atoms))
        allocate(tmp_atom_types_1d(box%num_atoms))
        allocate(tmp_atom_original_1d(box%num_atoms))

        ! Apply sorted indices
        do i = 1, box%num_atoms
            tmp_atom_names_1d(i)   = atom_names_1d(sort_index(i))
            tmp_atom_charges_1d(i) = atom_charges_1d(sort_index(i))
            tmp_atom_x_1d(i)       = atom_x_1d(sort_index(i))
            tmp_atom_y_1d(i)       = atom_y_1d(sort_index(i))
            tmp_atom_z_1d(i)       = atom_z_1d(sort_index(i))
            tmp_atom_types_1d(i)   = atom_types_1d(sort_index(i))
            tmp_atom_original_1d(i) = atom_original_1d(sort_index(i))
        end do

        ! Update global arrays with sorted data
        do i = 1, box%num_atoms
            atom_names_1d(i) = tmp_atom_names_1d(i)
            atom_charges_1d(i) = tmp_atom_charges_1d(i)
            atom_x_1d(i) = tmp_atom_x_1d(i)
            atom_y_1d(i) = tmp_atom_y_1d(i)
            atom_z_1d(i) = tmp_atom_z_1d(i)
            atom_types_1d(i) = tmp_atom_types_1d(i)
            atom_ids_1d(i) = i
            atom_original_1d(i) = tmp_atom_original_1d(i)
        end do

        call DetectResiduePattern(box, tmp_atom_types_1d)

        if (allocated(tmp_atom_names_1d)) deallocate(tmp_atom_names_1d)
        if (allocated(tmp_atom_charges_1d)) deallocate(tmp_atom_charges_1d)
        if (allocated(tmp_atom_x_1d)) deallocate(tmp_atom_x_1d)
        if (allocated(tmp_atom_y_1d)) deallocate(tmp_atom_y_1d)
        if (allocated(tmp_atom_z_1d)) deallocate(tmp_atom_z_1d)
        if (allocated(tmp_atom_types_1d)) deallocate(tmp_atom_types_1d)
        if (allocated(tmp_atom_original_1d)) deallocate(tmp_atom_original_1d)

    end subroutine SortAtomsByOriginalID

    subroutine DetectResiduePattern(box, tmp_atom_types_1d)

        implicit none

        type(type_box), intent(inout) :: box
        integer, allocatable, intent(in) :: tmp_atom_types_1d(:)

        integer :: idx, i, j, residue
        logical :: found
        integer, allocatable :: cpt_per_residue(:)

        ! one counter per residue
        allocate(cpt_per_residue(nb%type_residue))
        cpt_per_residue = 1

        do idx = 1, box%num_atoms
            found = .false.
            residue = -1

            ! find which residue this atom belongs to
            do i = 1, nb%type_residue
                do j = 1, nb%types_per_residue(i)
                    if (res%types_2d(i,j) == tmp_atom_types_1d(idx)) then
                        found = .true.
                        residue = i
                        exit  ! break from j-loop
                    end if
                end do
                if (found) exit  ! break from i-loop
            end do

            if (found) then
                nb%types_pattern(residue, cpt_per_residue(residue)) = tmp_atom_types_1d(idx)
                cpt_per_residue(residue) = cpt_per_residue(residue) + 1
                if (cpt_per_residue(residue) > nb%atom_in_residue(residue)) then
                    cpt_per_residue(residue) = 1
                end if
            end if
        end do

    end subroutine DetectResiduePattern

    logical function name_in_expected(atom_name, n_atoms, names2d, row)
        character(len=*), intent(in) :: atom_name
        integer, intent(in) :: n_atoms, row
        character(len=*), dimension(:,:), intent(in) :: names2d
        integer :: k

        name_in_expected = .false.
        do k = 1, n_atoms
            if (trim(atom_name) == trim(names2d(row,k))) then
                name_in_expected = .true.
                return
            end if
        end do
    end function name_in_expected

    !-----------------------------------------------------------------------
    !  DetectMolecules
    !
    !  Purpose:
    !    Scan through the list of atoms in the simulation box and group them
    !    into molecules (residues) according to predefined residue templates.
    !-----------------------------------------------------------------------
    subroutine DetectMolecules(box)

        implicit none

        ! Input parameters
        type(type_box), intent(inout) :: box

        ! Local variables
        integer :: i, j, k, cpt                 ! loop counters
        integer :: detected_nb_max_molecule     ! track max molecules found
        character(len=200) :: formatted_msg     ! buffer for error/info messages

        detected_nb_max_molecule = 0

        ! Initialize arrays inside box
        box%atom_types = 0
        box%atom_ids = 0
        box%atom_names = ""
        box%atom_charges = 0.0D0

        ! Loop over residue types
        do i = 1, nb%type_residue

            box%num_residues(i) = 0
            k = 1

            ! Scan atoms until we reach the end of the box
            do
                if (k > box%num_atoms) exit ! stop when no atoms left

                ! Reset atom counter for matching residue pattern
                cpt = 1

                ! If atom at position k matches the first atom of residue i
                if (atom_types_1d(k) == nb%types_pattern(i, cpt)) then

                    ! Ensure there are enough atoms left to build a full residue
                    if (k + nb%atom_in_residue(i) - 1 > box%num_atoms) then
                        call AbortRun("Not enough atoms left in box to complete residue type ")
                    end if

                    ! Copy all atoms belonging to this residue
                    do j = 1, nb%atom_in_residue(i)
                        box%atom_types(i, j) = atom_types_1d(k)
                        box%atom_ids(i, j) = atom_original_1d(k)
                        box%atom_names(i, j) = atom_names_1d(k)
                        box%atom_charges(i, j) = atom_charges_1d(k)

                        ! Check atom order against residue pattern
                        if (input%is_active(i) == 1) then
                            if (atom_types_1d(k) /= nb%types_pattern(i, cpt)) then
                                call AbortRun("Issue with atom order in data file")
                            end if
                        end if

                        ! Move to next atom and increment residue pattern index
                        k = k + 1
                        cpt = cpt + 1

                    end do

                    ! Increment residue counter for type i
                    box%num_residues(i) = box%num_residues(i) + 1

                else
                    k = k + 1
                end if
            end do

            ! Track maximum number of molecules detected among all types
            if (box%num_residues(i) > detected_nb_max_molecule) detected_nb_max_molecule = box%num_residues(i)

        end do

        ! Global check: maximum allowed molecules
        if (detected_nb_max_molecule > NB_MAX_MOLECULE) then
            write(formatted_msg, '(A, I0)') "The number of molecules exceeds the maximum allowed = ", NB_MAX_MOLECULE
            call InfoMessage(trim(formatted_msg))
            write(formatted_msg, '(A, I0)') "Number of molecules found = ", detected_nb_max_molecule
            call AbortRun(trim(formatted_msg) // " Please increase 'NB_MAX_MOLECULE' in src/parameters.f90 and recompile.", 11)
        end if

    end subroutine DetectMolecules

    !-----------------------------------------------------------------------
    !  RepairActiveMolecules
    !
    !  Purpose:
    !    Loop through all molecules in the box and, for those marked as active,
    !    extract their coordinates, apply a repair procedure, and write the
    !    updated coordinates back into the global arrays.
    !-----------------------------------------------------------------------
    subroutine RepairActiveMolecules(box)

        implicit none

        ! Input parameters
        type(type_box), intent(inout) :: box

        ! Local variables
        integer :: natoms, nmolecules             ! atoms per residue, molecules per type
        integer :: i, j, k, l, m, k_temp          ! loop counters and atom index tracking
        real(real64), allocatable :: tmp_x(:)     ! temporary x coordinates
        real(real64), allocatable :: tmp_y(:)     ! temporary y coordinates
        real(real64), allocatable :: tmp_z(:)     ! temporary z coordinates
        real(real64) :: dist                      ! Intramoleclar distance for sanity check

        ! Allocate temporary arrays for one molecule's atom positions
        allocate(tmp_x(maxval(nb%atom_in_residue)))
        allocate(tmp_y(maxval(nb%atom_in_residue)))
        allocate(tmp_z(maxval(nb%atom_in_residue)))

        k = 1  ! Global atom index

        ! Loop over all residue types
        do i = 1, nb%type_residue

            natoms = nb%atom_in_residue(i)      ! number of atoms per molecule of type i
            nmolecules = box%num_residues(i)    ! number of molecules of this type

            ! Loop over each molecule of this residue type
            do j = 1, nmolecules

                if (input%is_active(i) == 1) then ! ACTIVE molecule → repair coordinates

                    k_temp = k ! Save starting index for this molecule

                    ! Copy coordinates into temporary arrays
                    do l = 1, natoms
                        tmp_x(l) = atom_x_1d(k)
                        tmp_y(l) = atom_y_1d(k)
                        tmp_z(l) = atom_z_1d(k)
                        k = k + 1
                    end do

                    ! Repair the molecule if its marked as active
                    call RepairMolecule(tmp_x, tmp_y, tmp_z, natoms, box)

                    k = k_temp  ! Reset k to overwrite the same atoms

                    ! Sanity check : make sure distance between atoms of a same
                    ! molecules are not too large, and are not overlapping
                    do l = 1, natoms-1
                        do m = l+1, natoms
                            dist = sqrt((tmp_x(l)-tmp_x(m))**2 + &
                                        (tmp_y(l)-tmp_y(m))**2 + &
                                        (tmp_z(l)-tmp_z(m))**2)
                            if (dist > 10.0_real64) then
                                call WarnUser("Unusually large distance (> 1 nm) detected in active residue")
                            else if (dist < 1.0e-5_real64) then
                                call WarnUser("Overlapping atoms detected in molecule")
                            end if
                        end do
                    end do

                    ! Update coordinates with repaired values
                    do l = 1, natoms
                        atom_x_1d(k) = tmp_x(l)
                        atom_y_1d(k) = tmp_y(l)
                        atom_z_1d(k) = tmp_z(l)
                        k = k + 1
                    end do

                ! For inactve residue, update counter
                else
                    k = k + natoms ! INACTIVE molecule → skip without modification
                end if
            end do
        end do

        ! Free temporary arrays
        deallocate(tmp_x, tmp_y, tmp_z)

    end subroutine RepairActiveMolecules

    !===============================================================
    ! Subroutine: TransformCoordinate
    ! Purpose:
    !   Transform absolute atom coordinates into relative coordinates
    !   with respect to the molecule's center of mass (CoM).
    !===============================================================
    subroutine TransformCoordinate(box)

        use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

        implicit none

        ! Input parameters
        type(type_box), intent(inout) :: box

        ! Local variables
        integer :: i, j, k, l, m, cpt
        real(real64), dimension(3) :: com, original_com
        integer :: dim                      ! Dimension index (x, y, z)
        integer :: nb_res                   ! Number of molecules of this type

        ! Temporary arrays (assumed declared module-wide or can be moved here)
        real(real64), allocatable :: tmp_atom_masses_1d(:)
        real(real64), allocatable :: tmp_atom_x_1d(:)
        real(real64), allocatable :: tmp_atom_y_1d(:)
        real(real64), allocatable :: tmp_atom_z_1d(:)

        ! Allocate temporaries (size depends on max atoms per residue)
        allocate(tmp_atom_masses_1d(box%num_atoms))

        if (.not. allocated(res%masses_1d)) then
            allocate(res%masses_1d(primary%num_atomtypes))
        end if

        allocate(tmp_atom_x_1d(box%num_atoms))
        allocate(tmp_atom_y_1d(box%num_atoms))
        allocate(tmp_atom_z_1d(box%num_atoms))

        cpt = 1
        do i = 1, nb%type_residue

            nb_res = 0

            ! Build mass vector for the residue type i
            do j = 1, nb%atom_in_residue(i)
                do l = 1, nb%types_per_residue(i)
                    if (res%types_2d(i, l) == primary%atom_types(i, j)) then
                        res%masses_1d(primary%atom_types(i, j)) = primary%atom_masses(i, l)
                        tmp_atom_masses_1d = primary%atom_masses(i, l)
                        cpt = cpt + 1
                    end if
                end do
            end do

            ! Sanity check
            if (any(tmp_atom_masses_1d(1:nb%atom_in_residue(i)) <= 0.0_real64)) then
                call WarnUser("Zero or negative atomic mass detected in residue type")
            else if (sum(tmp_atom_masses_1d(1:nb%atom_in_residue(i))) < 1.0e-6_real64) then
                call WarnUser("Total molecular mass nearly zero in residue")
            end if

            ! Loop over all atoms and group them into molecules
            k = 1
            l = 1
            do while (k <= box%num_atoms)
                if (atom_types_1d(k) == box%atom_types(i, 1)) then

                    nb_res = nb_res + 1

                    ! Extract atom coordinates for this molecule
                    do j = 1, nb%atom_in_residue(i)
                        tmp_atom_x_1d(j) = atom_x_1d(k)
                        tmp_atom_y_1d(j) = atom_y_1d(k)
                        tmp_atom_z_1d(j) = atom_z_1d(k)
                        k = k + 1
                    end do

                    ! Compute Center of Mass
                    call ComputeCOM(tmp_atom_x_1d, tmp_atom_y_1d, tmp_atom_z_1d, &
                                    nb%atom_in_residue(i), tmp_atom_masses_1d, &
                                    com(1), com(2), com(3))

                    ! Save original CoM before applying periodic boundary conditions
                    original_com = com ! Save the CoM

                    ! Ensure CoM is inside simulation box
                    call ApplyPBC(com, box)

                    ! Sanity check (division by zero or corrupted input)
                    if (.not. all(ieee_is_finite(com))) then
                        call AbortRun("Invalid (NaN/Inf) CoM detected in residue")
                    end if

                    ! Sanity check, CoM outside simulation box
                    do dim = 1, 3
                        if (com(dim) < box%bounds(dim,1) .or. com(dim) > box%bounds(dim,2)) then
                            call WarnUser("Molecule COM outside simulation box")
                        end if
                    end do

                    ! Sanity check, CoM far from atoms in active molecule
                    if (input%is_active(i) == 1) then
                        if (minval(sqrt((tmp_atom_x_1d(1:nb%atom_in_residue(i)) - original_com(1))**2 + &
                                        (tmp_atom_y_1d(1:nb%atom_in_residue(i)) - original_com(2))**2 + &
                                        (tmp_atom_z_1d(1:nb%atom_in_residue(i)) - original_com(3))**2)) > 10.0_real64) then
                            call WarnUser("CoM unusually far from all atoms in residue type")
                        end if
                    end if

                    ! Store CoM position
                    do dim = 1, 3
                        box%mol_com(dim, i, l) = com(dim)
                    end do

                    ! Compute site offsets from CoM
                    do m = 1, nb%atom_in_residue(i)

                        box%site_offset(1, i, l, m) = tmp_atom_x_1d(m) - original_com(1)
                        box%site_offset(2, i, l, m) = tmp_atom_y_1d(m) - original_com(2)
                        box%site_offset(3, i, l, m) = tmp_atom_z_1d(m) - original_com(3)

                    end do
                    l = l + 1 ! Move to next molecule slot
                else
                    k = k + 1 ! Not start of this residue, keep scanning
                end if
            end do
        end do

        if (allocated(tmp_atom_masses_1d)) deallocate(tmp_atom_masses_1d)

    end subroutine TransformCoordinate

    !===============================================================
    !> Check if a given atom ID belongs to a specific residue.
    !---------------------------------------------------------------
    logical function isInResidue(box, res, atom_id)

        integer, intent(in) :: res, atom_id
        type(type_box), intent(inout) :: box

        integer :: n

        isInResidue = .false.

        do n = 1, nb%atom_in_residue(res)

            if (box%atom_ids(res, n) == atom_id) then
                isInResidue = .true.

                return
            end if
        end do

    end function isInResidue

    !===============================================================
    !> Get the local index of an atom within a residue.
    !---------------------------------------------------------------
    integer function atomIndexInResidue(box, res, atom_id)

        integer, intent(in) :: res, atom_id
        type(type_box), intent(inout) :: box

        integer :: n

        atomIndexInResidue = -1
        do n = 1, nb%atom_in_residue(res)
            if (box%atom_ids(res, n) == atom_id) then
                atomIndexInResidue = n
                return
            end if
        end do
    end function atomIndexInResidue

end module data_parser
