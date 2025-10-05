module helper_utils

    use simulation_state
    use constants
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none

contains

    !-------------------------------------------------------
    ! Function: repeat_char
    ! Creates a string of repeated characters
    !-------------------------------------------------------
    pure function repeat_char(ch, n) result(res)

        character(len=*), intent(in) :: ch
        integer, intent(in) :: n
        character(len=n) :: res
        integer :: i

        do i = 1, n
            res(i:i) = ch
        end do
    end function repeat_char

    !========================================================
    ! Function: RotationMatrix
    !
    ! Returns a 3x3 rotation matrix for a given axis (X=1, Y=2, Z=3)
    ! and rotation angle theta (radians).
    !
    ! Inputs:
    !   axis  - integer, rotation axis (1=X, 2=Y, 3=Z)
    !   theta - real(real64), rotation angle in radians
    !
    ! Output:
    !   rotation_matrix - real(real64) 3x3 rotation matrix
    !========================================================
    function RotationMatrix(axis, theta) result(rotation_matrix)

        implicit none
            
        integer, intent(in) :: axis               ! Rotation axis (1=X, 2=Y, 3=Z)
        real(real64), intent(in) :: theta         ! Rotation angle in radians
        real(real64) :: rotation_matrix(3,3)      ! 3x3 rotation matrix to be returned
        real(real64) :: cos_theta, sin_theta      ! Cosine and sine of theta

        ! Compute trigonometric values
        cos_theta = cos(theta)
        sin_theta = sin(theta)

        ! Initialize as identity
        rotation_matrix = zero
        rotation_matrix(1,1) = one
        rotation_matrix(2,2) = one
        rotation_matrix(3,3) = one

        ! Fill rotation matrix based on axis
        select case(axis)
        case(1) ! X-axis
            rotation_matrix(2,2) = cos_theta
            rotation_matrix(2,3) = -sin_theta
            rotation_matrix(3,2) = sin_theta
            rotation_matrix(3,3) = cos_theta
        case(2) ! Y-axis
            rotation_matrix(1,1) = cos_theta
            rotation_matrix(1,3) = sin_theta
            rotation_matrix(3,1) = -sin_theta
            rotation_matrix(3,3) = cos_theta
        case(3) ! Z-axis
            rotation_matrix(1,1) = cos_theta
            rotation_matrix(1,2) = -sin_theta
            rotation_matrix(2,1) = sin_theta
            rotation_matrix(2,2) = cos_theta
        end select

    end function RotationMatrix

    !--------------------------------------------------------------------
    ! present_or_false
    !
    ! Utility function to safely handle optional logical arguments.
    !
    ! If the optional argument `opt_flag` is present, its value is returned.
    ! If it is not present, the function returns `.false.` by default.
    !--------------------------------------------------------------------
    pure logical function present_or_false(opt_flag)

        logical, intent(in), optional :: opt_flag
    
        if (present(opt_flag)) then
            present_or_false = opt_flag
        else
            present_or_false = .false.
        end if
    
    end function present_or_false

    !--------------------------------------------------------------------
    ! amplitude_squared
    !
    ! Purpose:
    !   Computes the **squared modulus** (magnitude squared) of a complex number.
    !
    ! Description:
    !   For a complex number z = x + i*y, the squared modulus is defined as:
    !
    !       |z|^2 = x^2 + y^2 = z * conjg(z)
    !
    !   This is commonly used in Fourier/Ewald calculations to compute
    !   |A(k)|^2 for structure factor amplitudes.
    !
    ! Input:
    !   z : complex(real64)
    !       Complex number whose squared modulus is to be computed.
    !
    ! Output:
    !   val : real(real64)
    !       Squared modulus of z, i.e., |z|^2
    !
    ! Notes:
    !   - Pure function with no side effects.
    !   - Uses `z*conjg(z)` to compute |z|^2 efficiently.
    !--------------------------------------------------------------------
    pure function amplitude_squared(z) result(val)

        ! Input argument
        complex(real64), intent(in) :: z
        ! Output rgument
        real(real64) :: val

        val = real(z*conjg(z), kind=real64)

    end function amplitude_squared

end module helper_utils
