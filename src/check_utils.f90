module check_utils

    use output_utils
    use simulation_state
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none

contains

    subroutine CheckMolecule(box)

        implicit none

        ! Input parameters
        type(type_box), intent(inout) :: box

        ! Local variables
        integer :: i, j, k, dim
        real(real64), dimension(3) :: offset, com
        logical :: warned_large_offset
        warned_large_offset = .false.

        do i = 1, nb%type_residue
            do j = 1, box%num_residues(i)
                com = box%mol_com(:, i, j)

                ! Assert COM is inside the box
                do dim = 1,3
                    if (com(dim) < box%bounds(dim,1) .or. com(dim) > box%bounds(dim,2)) then
                        call WarnUser("Error: Molecule COM outside simulation box!")
                    end if
                end do

                ! Check atom offsets
                do k = 1, nb%atom_in_residue(i)
                    offset = box%site_offset(:, i, j, k)
                    if (any(abs(offset) > 10.0_real64) .and. (input%is_active(i)) == 1) then
                        if (.not. warned_large_offset) then
                            call WarnUser("One of the active molecules has an offset larger than 1 nanometer.")
                            warned_large_offset = .true.
                            write(*,*) box%site_offset(1, i, j, k), box%site_offset(2, i, j, k), box%site_offset(3, i, j, k)
                            stop 3
                        end if
                    end if
                end do
            end do
        end do

    end subroutine CheckMolecule

    !===============================================================
    ! Subroutine: CheckConsistency
    ! Purpose:
    !   Verify that atom masses in the primary simulation system
    !   and the external reservoir are consistent
    ! !===============================================================
    subroutine CheckConsistency()

        implicit none

        real(real64), parameter :: tol = 1.0d-6
        integer :: j, k
        character(len=128) :: msg

        if (has_reservoir) then
            ! Check that system and reservoir atom masses are consistent
            do j = 1, nb%type_residue
                do k = 1, nb%types_per_residue(j)
                    if (abs(primary%atom_masses(j, k) - reservoir%atom_masses(j, k)) > tol) then

                        ! Generic header warning
                        call WarnUser("Reservoir and system mass don't match.")

                        ! Detailed mismatch info
                        write(msg, '(A,I3,A,I3,A,F10.5,A,F10.5)') &
                            "Mismatch at residue ", j, ", site ", k, &
                            ": system mass = ", primary%atom_masses(j, k), &
                            ", reservoir mass = ", reservoir%atom_masses(j, k)
                        call WarnUser(trim(msg))

                        ! Extra guidance
                        call WarnUser("Make sure atom types are consistent between the two files.")
                    end if
                end do
            end do
        end if

    end subroutine CheckConsistency

end module check_utils
