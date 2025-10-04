module ewald_kvectors

    use simulation_state
    use output_utils
    use constants
    use geometry_utils
    use, intrinsic :: iso_fortran_env, only: real64

    implicit none

contains

    !-------------------------------------------------------------------
    ! This subroutine precomputes all valid reciprocal lattice vectors
    ! (kx, ky, kz) for the Ewald summation. It also precomputes the
    ! normalized squared magnitude k^2 for each vector.
    !
    ! This eliminates the need for triple nested loops during the
    ! reciprocal-space energy computation.
    !-------------------------------------------------------------------
    subroutine PrecomputeKVectors()

        implicit none

        integer :: kx_idx, ky_idx, kz_idx   ! Loop indices for each k-component
        integer :: count                    ! Counter for valid k-vectors
        real(real64) :: k_squared           ! Normalized squared magnitude of k-vector
        real(real64) :: k_squared_mag       ! Squared magnitude of the k-vector
        real(real64), dimension(3,3) :: kvec_matrix ! Columns are reciprocal lattice vectors b1, b2, b3

        ! Store reciprocal lattice vectors as columns of a 3x3 matrix
        kvec_matrix = TWOPI * reshape(primary%reciprocal, shape(kvec_matrix))

        ! Fill the array with the actual k-vectors
        count = 0
        do kx_idx = 0, ewald%kmax(1)
            do ky_idx = -ewald%kmax(2), ewald%kmax(2)
                do kz_idx = -ewald%kmax(3), ewald%kmax(3)

                    if (kx_idx == 0 .and. ky_idx == 0 .and. kz_idx == 0) cycle

                    ! Compute normalized k^2 again
                    k_squared = NormalizedK2(kx_idx, ky_idx, kz_idx, ewald%kmax)

                    ! Skip invalid k-vectors
                    if (.not. IsValidKVector(k_squared)) cycle

                    k_squared_mag = KVectorSquaredMag(kx_idx, ky_idx, kz_idx, kvec_matrix)

                    ! Increment counter and store k-vector components
                    count = count + 1
                    ewald%kvectors(count)%kx = kx_idx
                    ewald%kvectors(count)%ky = ky_idx
                    ewald%kvectors(count)%kz = kz_idx
                    ewald%kvectors(count)%k_squared = k_squared
                    ewald%kvectors(count)%k_squared_mag = k_squared_mag

                    ! Precompute form factor for current kx index: 1 for zero, 2 otherwise
                    ewald%form_factor(count) = FormFactor(kx_idx)
                end do
            end do
        end do

    end subroutine PrecomputeKVectors

    ! Computes the normalized squared magnitude of a reciprocal lattice 
    ! vector in index space.
    pure function NormalizedK2(kx, ky, kz, kmax) result(k_squared)

        ! Input arguments
        integer, intent(in) :: kx, ky, kz
        integer, intent(in) :: kmax(3)
        ! Output argument
        real(real64) :: k_squared

        ! Compute normalized squared magnitude in index space
        k_squared = (dble(kx)/dble(kmax(1)))**2 + (dble(ky)/dble(kmax(2)))**2 + (dble(kz)/dble(kmax(3)))**2
    
    end function NormalizedK2

    ! Computes the squared magnitude of a reciprocal lattice vector
    ! in Cartesian space given the lattice indices and the reciprocal
    ! lattice vectors. This is used to compute the Ewald weighting factor.
    pure function KVectorSquaredMag(kx, ky, kz, kvec_matrix) result(k2_mag)

        ! Input arguments
        integer, intent(in) :: kx, ky, kz
        real(real64), intent(in) :: kvec_matrix(3,3)
        ! Output argument
        real(real64) :: k2_mag
        ! Local argument
        real(real64) :: kvec(3)

        ! Build the 3D k-vector
        kvec = dble(kx) * kvec_matrix(:,1) + &
            dble(ky) * kvec_matrix(:,2) + &
            dble(kz) * kvec_matrix(:,3)

        ! Return squared magnitude
        k2_mag = dot_product(kvec, kvec)

    end function KVectorSquaredMag

    ! Computes the multiplicative form factor used in reciprocal-space
    ! Ewald summation to account for positive/negative k-vector symmetry.
    ! For a given k-index along one axis:
    !   - If the index is zero, the factor is 1
    !   - Otherwise (non-zero index), the factor is 2
    pure function FormFactor(idx) result(factor)

        ! Input arguments
        integer, intent(in) :: idx
        ! Output argument
        real(real64) :: factor

        if (idx == 0) then
            factor = one
        else
            factor = two
        end if

    end function FormFactor

    !---------------------------------------------------------------------
    ! Tests whether a k-vector (expressed as its normalized squared length)
    ! is valid for inclusion in the reciprocal-space Ewald summation.
    !
    ! A k-vector is considered valid if:
    !   - Its squared magnitude is larger than a small tolerance (not ~0),
    !   - Its squared magnitude does not exceed 1.0 in normalized space.
    !
    ! This prevents division by zero (k ≈ 0) and excludes k-vectors lying
    ! outside the unit sphere defined by the cutoff in normalized space.
    !---------------------------------------------------------------------
    pure function IsValidKVector(k_squared) result(valid)

        ! Input arguments
        real(real64), intent(in) :: k_squared
        ! Output argument
        logical :: valid

        ! Reject near-zero k-vectors (avoid singularity at k=0)
        ! and any vectors outside the normalized unit sphere.
        valid = (abs(k_squared) >= error) .and. (k_squared <= one)
        
    end function IsValidKVector

    !------------------------------------------------------------------------------
    ! Precomputes the reciprocal-space weighting factors for the Ewald summation.
    ! Loops over the precomputed reciprocal lattice vectors stored in ewald%kvectors
    ! and calculates the corresponding weighting factor for each vector:
    !    W(k) = exp(-|k|^2 / (4 * alpha^2)) / |k|^2
    ! This factor is used in the reciprocal-space sum of the Ewald method.
    ! Only done once at the start of the simulation.
    !------------------------------------------------------------------------------
    subroutine InitializeReciprocalWeights()

        implicit none

        ! Local variables
        integer :: idx                  ! Loop index over precomputed k-vectors
        real(real64) :: k_squared_mag   ! Squared magnitude of the k-vector
        real(real64) :: alpha_squared   ! Precompute alpha^2 for efficiency (alpha = screening parameter)

        ! Calculate the square of the screening parameter
        alpha_squared = ewald%alpha**2

        ! Loop over all precomputed reciprocal lattice vectors
        do idx = 1, ewald%num_kvectors

            ! Compute the reciprocal-space weighting factor for this k-vector
            k_squared_mag = ewald%kvectors(idx)%k_squared_mag
            ewald%recip_constants(idx) = exp(-k_squared_mag / (four * alpha_squared)) / k_squared_mag

        end do

    end subroutine InitializeReciprocalWeights

end module ewald_kvectors
