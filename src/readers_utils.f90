module readers_utils

    use simulation_state
    use output_utils
    use, intrinsic :: iso_fortran_env, only: real64

    implicit none

contains

    subroutine ComputeCOM(x, y, z, nb_atoms, mass, com_x, com_y, com_z)

        implicit none

        integer, intent(in) :: nb_atoms
        real(real64), intent(in) :: x(nb_atoms), y(nb_atoms), z(nb_atoms), mass(nb_atoms)
        real(real64), intent(out) :: com_x, com_y, com_z
        real(real64) :: total_mass
        character(len=32) :: mass_str
        integer :: i

        com_x = 0.d0
        com_y = 0.d0
        com_z = 0.d0
        total_mass = 0.d0

        do i = 1, nb_atoms
            com_x = com_x + mass(i) * x(i)
            com_y = com_y + mass(i) * y(i)
            com_z = com_z + mass(i) * z(i)
            total_mass = total_mass + mass(i)
        end do

        if (total_mass <= 0.0d0) then
            write(mass_str, '(F12.5)') total_mass
            call AbortRun("Total mass is zero or negative: " // trim(mass_str), 1)
        end if

        if (total_mass > 0.d0) then
            com_x = com_x / total_mass
            com_y = com_y / total_mass
            com_z = com_z / total_mass
        end if

    end subroutine ComputeCOM

    ! Read header info (e.g., number of atoms, atom types) from LAMMPS data file
    subroutine ReadLMPHeaderInfo(INFILE, box)

        implicit none

        ! Input parameters
        integer, intent(in) :: INFILE
        type(type_box), intent(inout) :: box

        ! Local variables
        character(len=256) :: line, trimmed_line
        integer :: ios, tmp_int
        logical :: found_atoms, found_atom_types

        ! Initialize outputs and flags
        found_atoms = .false.
        found_atom_types = .false.

        ! Rewind to the start of the file (optional, if you want to ensure starting from the top)
        rewind(INFILE)

        do
            ! Read a line from the file
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

            ! Try to read the line for "atoms" or "atom types"
            ! Use a simple string search for robustness
            if (index(trimmed_line, "atoms") > 0 .and. .not. found_atoms) then
                read(trimmed_line, *, IOSTAT=ios) tmp_int
                if (ios == 0) then
                    box%num_atoms = tmp_int
                    found_atoms = .true.
                end if
            else if (index(trimmed_line, "atom types") > 0 .and. .not. found_atom_types) then
                read(trimmed_line, *, IOSTAT=ios) tmp_int
                if (ios == 0) then
                    box%num_atomtypes = tmp_int
                    found_atom_types = .true.
                end if
            else if (index(trimmed_line, "bonds") > 0) then
                read(trimmed_line, *, IOSTAT=ios) tmp_int
                if (ios == 0) then
                    box%num_bonds = tmp_int
                end if
            else if (index(trimmed_line, "bond types") > 0) then
                read(trimmed_line, *, IOSTAT=ios) tmp_int
                if (ios == 0) then
                    box%num_bondtypes = tmp_int
                end if
            else if (index(trimmed_line, "angles") > 0) then
                read(trimmed_line, *, IOSTAT=ios) tmp_int
                if (ios == 0) then
                    box%num_angles = tmp_int
                end if
            else if (index(trimmed_line, "angle types") > 0) then
                read(trimmed_line, *, IOSTAT=ios) tmp_int
                if (ios == 0) then
                    box%num_angletypes = tmp_int
                end if
            else if (index(trimmed_line, "dihedrals") > 0) then
                read(trimmed_line, *, IOSTAT=ios) tmp_int
                if (ios == 0) then
                    box%num_dihedrals = tmp_int
                end if
            else if (index(trimmed_line, "dihedral types") > 0) then
                read(trimmed_line, *, IOSTAT=ios) tmp_int
                if (ios == 0) then
                    box%num_dihedraltypes = tmp_int
                end if
            else if (index(trimmed_line, "impropers") > 0) then
                read(trimmed_line, *, IOSTAT=ios) tmp_int
                if (ios == 0) then
                    box%num_impropers = tmp_int
                end if
            else if (index(trimmed_line, "improper types") > 0) then
                read(trimmed_line, *, IOSTAT=ios) tmp_int
                if (ios == 0) then
                    box%num_impropertypes = tmp_int
                end if
            end if
        end do

        ! Check if we found the required values
        if (.not. found_atoms)      call WarnUser("Number of atoms not found in header")
        if (.not. found_atom_types) call WarnUser("Number of atom types not found in header")

    end subroutine ReadLMPHeaderInfo

    !-----------------------------------------------------------------------------
    ! subroutine ParseLAMMPSBox(INFILE, box)
    !
    ! Parses the simulation box geometry from a LAMMPS data file, extracting both
    ! orthogonal and triclinic box parameters.
    ! Calculate 3x3 box matrix in lower-triangular convention (MDAnalysis convention)
    !-----------------------------------------------------------------------------
    subroutine ParseLAMMPSBox(INFILE, box)

        implicit none

        ! Input parameters
        type(type_box), intent(inout) :: box
        integer, intent(in) :: INFILE

        ! Local variables
        character(len=256) :: line
        character(len=32) :: keyword1, keyword2
        character(len=8) :: s1, s2, s3
        real(real64) :: tmp1, tmp2, tmp3
        real(real64) :: lx, ly, lz
        real(real64) :: xy, xz, yz
        real(real64) :: zero = 0.0_real64
        integer :: ios

        box%is_triclinic = .false.
        box%tilt(:) = zero

        ! Initialize bounds to a tiny value to indicate "not set"
        box%bounds(:,:) = 0.0_real64

        do
            read(INFILE, '(A)', IOSTAT=ios) line
            if (ios /= 0) exit

            ! Read xlo xhi
            read(line, *, IOSTAT=ios) tmp1, tmp2, keyword1, keyword2
            if (ios == 0) then
                if (TRIM(keyword1)//' '//TRIM(keyword2) == 'xlo xhi') then
                    box%bounds(1, 1) = tmp1
                    box%bounds(1, 2) = tmp2
                end if
            end if

            ! Read ylo yhi
            read(line, *, IOSTAT=ios) tmp1, tmp2, keyword1, keyword2
            if (ios == 0) then
                if (TRIM(keyword1)//' '//TRIM(keyword2) == 'ylo yhi') then
                    box%bounds(2, 1) = tmp1
                    box%bounds(2, 2) = tmp2
                end if
            end if

            ! Read zlo zhi
            read(line, *, IOSTAT=ios) tmp1, tmp2, keyword1, keyword2
            if (ios == 0) then
                if (TRIM(keyword1)//' '//TRIM(keyword2) == 'zlo zhi') then
                    box%bounds(3, 1) = tmp1
                    box%bounds(3, 2) = tmp2
                end if
            end if

            ! xy xz yz (triclinic tilt factors)
            read(line, *, IOSTAT=ios) tmp1, tmp2, tmp3, s1, s2, s3
            if (ios == 0 .and. TRIM(s1)//' '//TRIM(s2)//' '//TRIM(s3) == 'xy xz yz') then
                box%tilt(1) = tmp1
                box%tilt(2) = tmp2
                box%tilt(3) = tmp3
                box%is_triclinic = .true.
            end if

        end do

        ! Make sure that box size were provided
        if (abs(box%bounds(1,1)) < 1.0e-11_real64 .or. abs(box%bounds(1,2)) < 1.0e-11_real64) then
            call AbortRun("ParseLAMMPSBox: xlo xhi not found in input file!")
        end if

        if (abs(box%bounds(2,1)) < 1.0e-11_real64 .or. abs(box%bounds(2,2)) < 1.0e-11_real64) then
            call AbortRun("ParseLAMMPSBox: ylo yhi not found in input file!")
        end if

        if (abs(box%bounds(3,1)) < 1.0e-11_real64 .or. abs(box%bounds(3,2)) < 1.0e-11_real64) then
            call AbortRun("ParseLAMMPSBox: zlo zhi not found in input file!")
        end if

        ! Box lengths
        lx = box%bounds(1,2) - box%bounds(1,1)
        ly = box%bounds(2,2) - box%bounds(2,1)
        lz = box%bounds(3,2) - box%bounds(3,1)

        ! Triclinic tilt
        xy = box%tilt(1)
        xz = box%tilt(2)
        yz = box%tilt(3)

        ! Assign row-major matrix following
        box%matrix(1,:) = [lx, zero, zero]  ! row 1 (a vector)
        box%matrix(2,:) = [xy, ly, zero]    ! row 2 (b vector)
        box%matrix(3,:) = [xz, yz, lz]      ! row 3 (c vector)

    end subroutine ParseLAMMPSBox

    !-----------------------------------------------------------
    ! Subroutine: RepairMolecule
    ! Ensures that bonded atoms in a molecule remain contiguous across
    ! periodic boundaries. Works for cubic, orthorhombic, and triclinic boxes.
    !
    ! Note: This routine does NOT re-center the molecule in the simulation box.
    !       It only adjusts atomic positions so that the distance between
    !       bonded atoms is minimized under PBC.
    !-----------------------------------------------------------
    subroutine RepairMolecule(atom_x, atom_y, atom_z, nb_atoms, box)

        use, intrinsic :: iso_fortran_env, only: real64
        implicit none

        ! Input parameters
        integer, intent(in) :: nb_atoms             ! Number of atoms in the molecule
        type(type_box), intent(inout) :: box        ! Simulation box: contains box type, matrix, and reciprocal matrix for PBC handling
        real(real64), intent(inout) :: atom_x(nb_atoms) ! x-coordinates of atoms, updated to ensure molecular contiguity under PBC
        real(real64), intent(inout) :: atom_y(nb_atoms) ! y-coordinates of atoms, updated to ensure molecular contiguity under PBC
        real(real64), intent(inout) :: atom_z(nb_atoms) ! z-coordinates of atoms, updated to ensure molecular contiguity under PBC

        ! Local variables
        real(real64), dimension(3) :: delta_r_cart  ! Cartesian displacement vector: (dx, dy, dz) between current and previous atom
        real(real64), dimension(3) :: delta_r_frac  ! Fractional displacement vector: delta_r_cart transformed into fractional coordinates for PBC wrapping
        real(real64), dimension(3) :: ref_coords    ! Cartesian coordinates of the reference atom (i-1) in the molecule
        real(real64), dimension(3) :: curr_coords   ! Cartesian coordinates of the current atom (i), updated after PBC correction
        integer :: idim, iatom                      ! Loop index for iteration

        do iatom = 2, nb_atoms
            ! Take reference atom (iatom-1) and current atom (iatom)
            ref_coords  = [atom_x(iatom-1), atom_y(iatom-1), atom_z(iatom-1)]
            curr_coords = [atom_x(iatom), atom_y(iatom), atom_z(iatom)]

            ! Raw displacement between the two atoms (may cross PBC boundary)
            delta_r_cart = curr_coords - ref_coords

            if (box%type == 1 .or. box%type == 2) then
                
                ! Cubic/orthorhombic: simple wrapping along box edges
                do idim = 1, 3
                    delta_r_cart(idim) = wrap_nearest(delta_r_cart(idim), box%matrix(idim, idim))
                end do

            else if (box%type == 3) then
                ! Triclinic: use fractional coords
                delta_r_frac = matmul(box%reciprocal, delta_r_cart)

                ! Wrap into [-0.5, 0.5)
                do idim = 1, 3
                    delta_r_frac(idim) = wrap_nearest(delta_r_frac(idim), 1.0_real64)
                end do

                ! Back to Cartesian
                delta_r_cart = matmul(box%matrix, delta_r_frac)
            else
                call AbortRun("ERROR in RepairMolecule: box%type is invalid.", 1)
            end if

            ! Repair atom position
            curr_coords = ref_coords + delta_r_cart
            atom_x(iatom) = curr_coords(1)
            atom_y(iatom) = curr_coords(2)
            atom_z(iatom) = curr_coords(3)

        end do

    end subroutine RepairMolecule

    ! wrap_nearest:
    !   Given a displacement x and a box length boxlen, returns the
    !   equivalent displacement wrapped into the nearest image of the
    !   simulation box.
    pure real(real64) function wrap_nearest(x, boxlen) result(xwrap)

        implicit none
        real(real64), intent(in) :: x, boxlen
        ! Returns value in [-boxlen/2, boxlen/2)
        xwrap = modulo(x + 0.5_real64*boxlen, boxlen) - 0.5_real64*boxlen

    end function wrap_nearest


    subroutine sort_pairs(n, p1, p2, eps, sig)
        integer, intent(in) :: n
        integer, intent(inout) :: p1(n), p2(n)
        real(real64), intent(inout) :: eps(n), sig(n)
        integer :: i, j
        integer :: temp1, temp2
        real(real64) :: tempe, temps
        do i = 1, n-1
            do j = i+1, n
                if ((p1(i) > p1(j)) .or. ((p1(i) == p1(j)) .and. (p2(i) > p2(j)))) then
                    temp1 = p1(i)
                    temp2 = p2(i)
                    tempe = eps(i)
                    temps = sig(i)

                    p1(i) = p1(j)
                    p2(i) = p2(j)
                    eps(i) = eps(j)
                    sig(i) = sig(j)

                    p1(j) = temp1
                    p2(j) = temp2
                    eps(j) = tempe
                    sig(j) = temps
                end if
            end do
        end do
    end subroutine sort_pairs

end module readers_utils

