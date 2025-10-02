module random_utils

    use, intrinsic :: iso_fortran_env, only: real64

    implicit none

    contains

    ! Returns double precision random number uniformly distributed in [0,1)
    function rand_uniform() result(val)
        real(real64) :: val
        call random_number(val)
    end function rand_uniform

    ! Seeds the RNG: if seed=0, generate a random seed from system_clock
    subroutine seed_rng(seed)
        integer, intent(in) :: seed
        integer :: n, i, useed
        integer, allocatable :: put(:)
        integer :: count, count_rate, count_max

        if (seed == 0) then
            ! Get system clock count as a seed
            call system_clock(count, count_rate, count_max)
            useed = abs(mod(count, 2147483647)) ! keep in 32-bit integer range
        else
            useed = seed
        end if

        call random_seed(size=n)
        allocate(put(n))

        put = useed + 37 * [(i-1, i=1,n)]
        call random_seed(put=put)

        deallocate(put)
    end subroutine seed_rng

end module random_utils
