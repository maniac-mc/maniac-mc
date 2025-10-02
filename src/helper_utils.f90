module helper_utils

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

end module helper_utils
