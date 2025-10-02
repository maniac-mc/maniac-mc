module geometry_utils

    use simulation_state
    use output_utils
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none

contains

    !-----------------------------------------------------------
    ! Subroutine: PrepareSimulationBox
    ! Purpose: Initialize a simulation box by determining its 
    !          symmetry, computing geometric properties, and 
    !          calculating the inverse matrix for reciprocal-space use.
    !-----------------------------------------------------------
    subroutine PrepareSimulationBox(box)

        implicit none

        type(type_box), intent(inout) :: box

        character(200) :: formatted_msg ! Buffer for formatted output messages

        call DetermineBoxSymmetry(box)

        ! Print determined box type
        call LogMessage("====== Simulation preparation ======")
        call LogMessage("")

        select case (box%type)
            case (1)
                call LogMessage("Box symmetry type: Cubic")
            case (2)
                call LogMessage("Box symmetry type: Orthorhombic")
            case (3)
                call LogMessage("Box symmetry type: Triclinic")
            case default
                write(formatted_msg, '(A, I0)') 'Box symmetry type determined: ', box%type
                call LogMessage(formatted_msg)
        end select

        call ComputeCellProperties(box)

        write(formatted_msg, '(A, F20.4)') 'Cell volume (Å^3): ', box%volume
        call LogMessage(formatted_msg)

        call ComputeInverse(box)

    end subroutine PrepareSimulationBox

    !-----------------------------------------------------------
    ! Subroutine: DetermineBoxSymmetry
    ! Determine the symmetry type of a simulation box.
    ! box type = 1 for cubic box,
    ! box type = 2 for orthorhombic box,
    ! box type = 3 for triclinic box.
    !-----------------------------------------------------------
    subroutine DetermineBoxSymmetry(box)

        implicit none

        ! Input parameters
        type(type_box), intent(inout) :: box

        ! Local variables
        real(real64), dimension(6) :: off_diag

        ! Collect off-diagonal elements
        off_diag = [box%matrix(1,2), box%matrix(1,3), &
            box%matrix(2,1), box%matrix(2,3), &
            box%matrix(3,1), box%matrix(3,2)]

        if (maxval(abs(off_diag)) > 1.0D-12) then
            box%type = 3
        else if (abs(box%matrix(1,1) - box%matrix(2,2)) > 1.0D-12 .or. &
                 abs(box%matrix(1,1) - box%matrix(3,3)) > 1.0D-12) then
            box%type = 2
        else
            box%type = 1
        end if
        
    end subroutine DetermineBoxSymmetry

    !-----------------------------------------------------------
    ! Subroutine: ComputeCellProperties
    ! Purpose: Calculate lengths, cosines of angles,
    !          and perpendicular widths of a simulation cell.
    !-----------------------------------------------------------
    subroutine ComputeCellProperties(box)

        implicit none

        ! Input parameters
        type(type_box), intent(inout) :: box

        ! Local variables
        real(real64), dimension(3) :: vec_lengths ! Lengths of the three cell vectors
        real(real64), dimension(3) :: axb ! Cross product of vector a and b
        real(real64), dimension(3) :: bxc ! Cross product of vector b and c
        real(real64), dimension(3) :: cxa ! Cross product of vector c and a
        real(real64) :: a(3), b(3), c(3) ! Local vectors for cross products

        ! Calculate lengths of cell vectors
        box%metrics(1:3) = sqrt(sum(box%matrix(:, 1:3)**2, dim=1))

        ! Copy cell vectors from the box matrix
        a = box%matrix(:, 1)
        b = box%matrix(:, 2)
        c = box%matrix(:, 3)

        ! Calculate vector lengths
        vec_lengths(1) = sqrt(sum(a**2))
        vec_lengths(2) = sqrt(sum(b**2))
        vec_lengths(3) = sqrt(sum(c**2))

        ! Calculate cosines of angles using dot products divided by product of lengths
        box%metrics(4) = dot_product(a, b) / (vec_lengths(1) * vec_lengths(2))
        box%metrics(5) = dot_product(a, c) / (vec_lengths(1) * vec_lengths(3))
        box%metrics(6) = dot_product(b, c) / (vec_lengths(2) * vec_lengths(3))

        ! Calculate cross products
        axb = CrossProduct(a, b)
        bxc = CrossProduct(b, c)
        cxa = CrossProduct(c, a)

        ! Calculate volume of cell
        box%volume = abs(dot_product(box%matrix(:,1), bxc))

        ! Calculate cell perpendicular widths
        box%metrics(7) = box%volume / norm2(bxc)
        box%metrics(8) = box%volume / norm2(cxa)
        box%metrics(9) = box%volume / norm2(axb)

    end subroutine ComputeCellProperties

    !-----------------------------------------------------------
    ! Subroutine: ApplyPBC
    ! Wraps a 3D position into the simulation box
    ! Works for orthogonal or triclinic boxes with any lower bounds
    !-----------------------------------------------------------
    subroutine ApplyPBC(pos, box)

        implicit none

        ! Input parameters
        type(type_box), intent(in) :: box
        real(real64), dimension(3), intent(inout) :: pos

        ! Local variables
        real(real64), dimension(3) :: frac_coords    ! Fractional coordinates inside the box
        real(real64), dimension(3) :: box_lengths    ! Box lengths along x, y, z (diagonal)
        integer :: dim                               ! Loop index

        ! Orthogonal / Rectangular Box
        if (.not. box%is_triclinic) then
            do dim = 1,3
                box_lengths(dim) = box%matrix(dim,dim)
                pos(dim) = box%bounds(dim,1) + modulo(pos(dim) - box%bounds(dim,1), box_lengths(dim))
            end do

        ! Triclinic Box
        else
            ! Convert Cartesian -> fractional coordinates
            frac_coords = matmul(box%reciprocal, pos - box%bounds(:,1))

            ! Wrap fractional coords into [0,1)
            do dim = 1,3
                frac_coords(dim) = modulo(frac_coords(dim), 1.0_real64)
            end do

            ! Convert back to Cartesian
            pos = box%bounds(:,1) + matmul(box%matrix, frac_coords)
        end if

    end subroutine ApplyPBC

    !-----------------------------------------------------------
    ! Subroutine: WrapIntoBox
    ! Ensures pos (x, y, z) are inside the simulation box
    ! Convention: [-L/2, L/2] for cubic/orthorhombic,
    !             [-0.5, 0.5] in fractional space for triclinic
    !-----------------------------------------------------------
    subroutine WrapIntoBox(pos, box)

        implicit none

        ! Input parameters
        real(real64), dimension(3), intent(inout) :: pos
        type(type_box), intent(inout) :: box

        ! Local variables
        real(real64), dimension(3) :: fractional_pos
        integer :: dim

        ! Cubic or Orthorhombic box
        if (box%type == 1 .or. box%type == 2) then

            do dim = 1, 3
                ! Shift position into [-L/2, L/2]
                pos(dim) = pos(dim) - box%matrix(dim, dim) * &
                        nint(pos(dim) / box%matrix(dim, dim))
            end do

        ! Triclinic box
        else if (box%type == 3) then

            ! Convert Cartesian to fractional coordinates
            fractional_pos = matmul(box%reciprocal, pos)

            ! Wrap fractional coordinates into [-0.5, 0.5)
            do dim = 1, 3
                fractional_pos(dim) = fractional_pos(dim) - nint(fractional_pos(dim))
            end do

            ! Convert back to Cartesian
            pos = matmul(box%matrix, fractional_pos)

        end if

    end subroutine WrapIntoBox

    !----------------------------------------------------------------------------
    ! subroutine ComputeInverse
    ! Utility for computing the inverse and determinant of a 3×3 box matrix
    ! used in reciprocal space calculations
    !----------------------------------------------------------------------------
    subroutine ComputeInverse(box)

        implicit none

        ! Input parameters
        type(type_box), intent(inout) :: box

        ! Local variables
        real(real64), dimension(3,3) :: adjugate ! Adjugate (cofactor transpose) matrix of box_mat
        real(real64) :: reciprocal ! Reciprocal of the determinant of box_mat
        real(real64), dimension(3) :: a, b
        integer :: i, j ! Loop indices for matrix element iteration
        character(200) :: formatted_msg ! buffer for logging

        ! Numerical thresholds
        real(real64), parameter :: tol_det   = 1.0D-15  ! determinant tolerance
        real(real64), parameter :: tiny_val  = tiny(1.0_real64)  ! smallest safe value

        ! Compute the adjugate (cofactor transpose) using cross products
        ! First column of adjugate
        a = box%matrix(:, 2)
        b = box%matrix(:, 3)
        adjugate(:,1) = CrossProduct(a, b)

        ! Second column
        a = box%matrix(:, 3)
        b = box%matrix(:, 1)
        adjugate(:,2) = CrossProduct(a, b)

        ! Third column
        a = box%matrix(:, 1)
        b = box%matrix(:, 2)
        adjugate(:,3) = CrossProduct(a, b)

        ! Determinant from first row
        box%determinant = dot_product(box%matrix(:,1), adjugate(:,1))

        ! Check for denormal/underflow values in determinant
        if (abs(box%determinant) < tiny_val) then
            write(formatted_msg,'(A,1PE12.4)') "Error: Determinant fell into denormal/underflow range. det = ", box%determinant
            call LogMessage(formatted_msg)
            stop 1
        end if

        ! Calculate reciprocal of determinant if non-zero
        reciprocal = 0.0_real64
        if (abs(box%determinant) > tol_det) then  ! small tolerance to avoid near-zero
            reciprocal = 1.0_real64 / box%determinant
        else
            ! Handle singular matrix case if needed, e.g. print error or return
            call WarnUser("Determinant is zero or near zero, inverse not computed.")
            return
        end if
        
        ! Store inverse matrix in reciprocal
        do i = 1, 3
            do j = 1, 3
                box%reciprocal(i, j) = reciprocal * adjugate(i, j)
            end do
        end do

        return
    end subroutine ComputeInverse 

    ! Compute the cross product of two 3D vectors a and b
    function CrossProduct(a, b) result(c)
    
        real(real64), intent(in) :: a(3), b(3) ! Input vectors
        real(real64) :: c(3) ! Resulting cross product vector

        ! Cross product components
        c(1) = a(2)*b(3) - a(3)*b(2) ! x-component
        c(2) = a(3)*b(1) - a(1)*b(3) ! y-component
        c(3) = a(1)*b(2) - a(2)*b(1) ! z-component
    end function CrossProduct


    !---------------------------------------------------------------------
    !> Compute the minimum-image distance between two atoms in a periodic box
    !>
    !> This function calculates the Cartesian distance between two atoms,
    !> taking into account periodic boundary conditions (PBC). It handles
    !> three types of boxes:
    !>   1. Cubic
    !>   2. Orthorhombic
    !>   3. Triclinic (tilted boxes)
    !>
    !> For cubic and orthorhombic boxes, the minimum-image convention is
    !> applied per Cartesian component. For triclinic boxes, all 27
    !> neighboring periodic images are checked to ensure the true
    !> shortest distance is returned.
    !>
    !> @param box                Box structure containing box vectors,
    !>                           molecule centers, atom offsets, and box type
    !> @param residue_type_1     Residue type of the first atom
    !> @param molecule_index_1   Molecule index of the first atom
    !> @param atom_index_1       Atom index of the first atom
    !> @param residue_type_2     Residue type of the second atom
    !> @param molecule_index_2   Molecule index of the second atom
    !> @param atom_index_2       Atom index of the second atom
    !>
    !> @return distance          Minimum-image Cartesian distance between the two atoms
    !---------------------------------------------------------------------
    function ComputeDistance(box, residue_type_1, molecule_index_1, atom_index_1, &
                            residue_type_2, molecule_index_2, atom_index_2) result(distance)

        use, intrinsic :: iso_fortran_env, only: real64

        implicit none

        type(type_box), intent(in) :: box
        integer, intent(in) :: residue_type_1, molecule_index_1, atom_index_1
        integer, intent(in) :: residue_type_2, molecule_index_2, atom_index_2

        integer :: shift_x, shift_y, shift_z
        integer :: dim
        real(real64), dimension(3) :: delta
        real(real64), dimension(3) :: trial_delta
        real(real64) :: distance
        real(real64) :: min_dist2
        real(real64) :: trial_dist2

        ! Compute Cartesian difference between the two atoms
        delta = (box%mol_com(:,residue_type_2,molecule_index_2) + &
                box%site_offset(:,residue_type_2,molecule_index_2,atom_index_2)) - &
                (box%mol_com(:,residue_type_1,molecule_index_1) + &
                box%site_offset(:,residue_type_1,molecule_index_1,atom_index_1))

        ! Cubic or Orthorhombic box
        if (box%type == 1 .or. box%type == 2) then

            ! Apply PBC to delta vector directly
            do dim = 1,3
                delta(dim) = modulo(delta(dim) + 0.5_real64*box%matrix(dim,dim), &
                    box%matrix(dim,dim)) - 0.5_real64*box%matrix(dim,dim)
            end do

            ! Compute distance
            distance = sqrt(dot_product(delta, delta))

        ! Triclinic box
        else if (box%type == 3) then

            min_dist2 = huge(1.0_real64)
            do shift_x = -1,1
                do shift_y = -1,1
                    do shift_z = -1,1
                        trial_delta = delta + shift_x*box%matrix(:,1) + &
                                            shift_y*box%matrix(:,2) + &
                                            shift_z*box%matrix(:,3)
                        trial_dist2 = dot_product(trial_delta, trial_delta)
                        if (trial_dist2 < min_dist2) min_dist2 = trial_dist2
                    end do
                end do
            end do
            distance = sqrt(min_dist2)

        end if

    end function ComputeDistance

end module geometry_utils
