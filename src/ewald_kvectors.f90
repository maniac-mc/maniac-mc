module ewald_kvectors

    !----------------------------------------------------------------------
    ! Module: ewald_kvectors
    !
    ! Purpose:
    !   Provides routines to precompute all reciprocal lattice vectors
    !   and associated quantities for the reciprocal-space part of the
    !   Ewald summation method.
    !
    ! Features:
    !   - Precompute valid k-vectors and their squared magnitudes.
    !   - Precompute multiplicative form factors for symmetry.
    !   - Precompute reciprocal-space weighting constants W(k).
    !   - Avoid repeated loops over k-space during energy computations.
    !----------------------------------------------------------------------

    use simulation_state
    use output_utils
    use constants
    use geometry_utils
    use, intrinsic :: iso_fortran_env, only: real64

    implicit none

contains

    !-------------------------------------------------------------------
    ! PrecomputeValidReciprocalVectors
    !
    ! Purpose:
    !   Generates all valid reciprocal lattice vectors for the system
    !   and stores them in ewald%kvectors. Also precomputes the normalized
    !   squared magnitude (k^2) and Cartesian squared magnitude.
    !
    ! Description:
    !   Loops over all possible kx, ky, kz indices within the specified
    !   kmax bounds and stores only valid vectors:
    !     - k ≠ 0 (to avoid singularity)
    !     - k within normalized cutoff sphere
    !   Also precomputes the symmetry form factor for each vector:
    !     1 for zero index, 2 for non-zero.
    !-------------------------------------------------------------------
    subroutine PrecomputeValidReciprocalVectors()

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
                    k_squared = NormalizedKSquared(kx_idx, ky_idx, kz_idx, ewald%kmax)

                    ! Skip invalid k-vectors
                    if (.not. CheckValidReciprocalVector(k_squared)) cycle

                    k_squared_mag = ComputeCartesianKSquared(kx_idx, ky_idx, kz_idx, kvec_matrix)

                    ! Increment counter and store k-vector components
                    count = count + 1
                    ewald%kvectors(count)%kx = kx_idx
                    ewald%kvectors(count)%ky = ky_idx
                    ewald%kvectors(count)%kz = kz_idx
                    ewald%kvectors(count)%k_squared = k_squared
                    ewald%kvectors(count)%k_squared_mag = k_squared_mag

                    ! Precompute form factor for current kx index: 1 for zero, 2 otherwise
                    ewald%form_factor(count) = ComputeSymmetryFormFactor(kx_idx)
                end do
            end do
        end do

    end subroutine PrecomputeValidReciprocalVectors

    !-------------------------------------------------------------------
    ! NormalizedKSquared
    !
    ! Purpose:
    !   Computes the squared magnitude of a k-vector in normalized index space.
    !
    ! Formula:
    !   k^2_normalized = (kx/kmax_x)^2 + (ky/kmax_y)^2 + (kz/kmax_z)^2
    !
    ! Input:
    !   kx, ky, kz : integer indices of the reciprocal lattice vector
    !   kmax(3)    : maximum k-indices along each axis
    !
    ! Output:
    !   k_squared  : normalized squared magnitude
    !-------------------------------------------------------------------
    pure function NormalizedKSquared(kx, ky, kz, kmax) result(k_squared)

        ! Input arguments
        integer, intent(in) :: kx, ky, kz
        integer, intent(in) :: kmax(3)
        ! Output argument
        real(real64) :: k_squared

        ! Compute normalized squared magnitude in index space
        k_squared = (dble(kx)/dble(kmax(1)))**2 + &
                    (dble(ky)/dble(kmax(2)))**2 + &
                    (dble(kz)/dble(kmax(3)))**2

    end function NormalizedKSquared

    !-------------------------------------------------------------------
    ! ComputeCartesianKSquared
    !
    ! Purpose:
    !   Computes the squared magnitude of a reciprocal lattice vector
    !   in Cartesian space using the reciprocal lattice vectors.
    !
    ! Input:
    !   kx, ky, kz   : integer lattice indices
    !   kvec_matrix  : 3x3 matrix of reciprocal lattice vectors (columns)
    !
    ! Output:
    !   k2_mag       : squared magnitude |k|^2 in Cartesian units
    !-------------------------------------------------------------------
    pure function ComputeCartesianKSquared(kx, ky, kz, kvec_matrix) result(k2_mag)

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

    end function ComputeCartesianKSquared

    !-------------------------------------------------------------------
    ! ComputeSymmetryFormFactor
    !
    ! Purpose:
    !   Returns the multiplicative symmetry factor for a k-index along one axis.
    !   1 for zero index, 2 for non-zero.
    !
    ! Input:
    !   idx : integer index along an axis
    !
    ! Output:
    !   factor : multiplicative form factor
    !-------------------------------------------------------------------
    pure function ComputeSymmetryFormFactor(idx) result(factor)

        ! Input arguments
        integer, intent(in) :: idx
        ! Output argument
        real(real64) :: factor

        if (idx == 0) then
            factor = one
        else
            factor = two
        end if

    end function ComputeSymmetryFormFactor

    !-------------------------------------------------------------------
    ! CheckValidReciprocalVector
    !
    ! Purpose:
    !   Determines whether a normalized k-vector is valid for inclusion
    !   in the reciprocal-space Ewald summation.
    !
    ! Criteria:
    !   - Reject near-zero vectors (to avoid division by zero)
    !   - Exclude vectors outside the normalized unit sphere (k^2 > 1)
    !
    ! Input:
    !   k_squared : normalized squared magnitude
    !
    ! Output:
    !   valid     : .true. if vector is valid
    !-------------------------------------------------------------------
    pure function CheckValidReciprocalVector(k_squared) result(valid)

        ! Input arguments
        real(real64), intent(in) :: k_squared
        ! Output argument
        logical :: valid

        ! Reject near-zero k-vectors (avoid singularity at k=0)
        ! and any vectors outside the normalized unit sphere.
        valid = (abs(k_squared) >= error) .and. (k_squared <= one)

    end function CheckValidReciprocalVector

    !-------------------------------------------------------------------
    ! ComputeReciprocalWeights
    !
    ! Purpose:
    !   Precomputes reciprocal-space weighting factors W(k) for the
    !   Ewald sum, using the standard formula:
    !
    !       W(k) = exp(-|k|^2 / (4 * alpha^2)) / |k|^2
    !
    ! Description:
    !   Loops over all precomputed k-vectors and stores W(k) in
    !   ewald%recip_constants. This avoids recomputation during energy sums.
    !-------------------------------------------------------------------
    subroutine ComputeReciprocalWeights()

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

    end subroutine ComputeReciprocalWeights

end module ewald_kvectors
