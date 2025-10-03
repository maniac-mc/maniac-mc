module tabulated_utils

    use simulation_state
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none

contains

    ! Initialize all precomputed tables for the simulation
    subroutine PrecomputeTable()
        ! Example: only erfc/r table for now
        call InitializeTabulatedErfcR(erfc_r_table, ewald%alpha, input%real_space_cutoff)
    end subroutine PrecomputeTable

    ! Precomputes a tabulated version of the function erfc(alpha*r)/r for faster
    ! evaluation during simulations. Useful for Ewald direct-space Coulomb energy.
    subroutine InitializeTabulatedErfcR(table, alpha, r_cut)

        ! Input variables
        type(tabulated), intent(inout) :: table  ! Tabulated function structure to store x and f arrays
        real(real64), intent(in) :: alpha        ! Ewald screening parameter for erfc(alpha*r)/r
        real(real64), intent(in) :: r_cut        ! Maximum distance (real-space cutoff) for the table

        ! Local variables
        integer :: i                              ! Loop index over table points
        real(real64) :: r                         ! Current distance corresponding to table grid point

        ! Allocate arrays for grid points (x) and function values (f)
        ! Using global parameter TABULATED_POINTS from parameters module
        allocate(table%x(0:TABULATED_POINTS))
        allocate(table%f(0:TABULATED_POINTS))
        table%n = TABULATED_POINTS
        table%dx = r_cut / real(TABULATED_POINTS, real64)

        ! Fill the table
        do i = 0, TABULATED_POINTS
            r = i * table%dx
            table%x(i) = r
            if (r < error) then
                table%f(i) = two * alpha / sqrt(pi)
            else
                table%f(i) = erfc(alpha*r) / r
            end if
        end do

        ! Mark table as initialized
        table%initialized = .true.
    end subroutine InitializeTabulatedErfcR

    ! Linear interpolation in a tabulated function
    pure function LookupTabulated(table, r) result(f_r)

        ! Input variables
        type(tabulated), intent(in) :: table    ! Tabulated function structure to store x and f arrays
        real(real64), intent(in) :: r           ! Query distance

        ! Local variables
        real(real64) :: f_r                     ! Interpolated function value
        real(real64) :: t, f1, f2               ! Interpolation weight and function values
        integer :: i                            ! Index of lower grid point

        if (r <= zero) then
            f_r = table%f(0)
            return
        end if

        if (r >= table%x(table%n)) then
            f_r = zero
            return
        end if

        i = int(r / table%dx)
        f1 = table%f(i)
        f2 = table%f(i+1)
        t = (r - table%x(i)) / table%dx
        f_r = (one - t) * f1 + t * f2
    end function LookupTabulated

end module tabulated_utils
