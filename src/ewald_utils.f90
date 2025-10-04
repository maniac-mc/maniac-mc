module ewald_utils

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

    !--------------------------------------------------------------------
    ! Computes the phase vector for a single atom in reciprocal space.
    !
    !   For Ewald reciprocal-space calculations, the phase for each axis is
    !   computed as the dot product of the atom position with the corresponding
    !   reciprocal lattice vector, multiplied by 2π. This gives the argument
    !   for the complex exponential exp(i k · r) along each axis.
    !--------------------------------------------------------------------
    pure subroutine ComputeAtomPhase(atom_pos, reciprocal_box, phase)

        real(real64), intent(in) :: atom_pos(3)
        real(real64), intent(in) :: reciprocal_box(3,3)
        real(real64), intent(out) :: phase(3)

        integer :: i, j

        do i = 1, 3
            phase(i) = zero
            do j = 1, 3
                phase(i) = phase(i) + reciprocal_box(j,i) * atom_pos(j)
            end do
            phase(i) = TWOPI * phase(i)
        end do

    end subroutine ComputeAtomPhase

    !--------------------------------------------------------------------
    ! Computes 1D complex exponential phase factors along a single axis
    ! for a given atom.
    !
    !   The array phase_factor_axis stores exp(i * k * phase_component)
    !   for all k-indices, including negative k-values. These factors are
    !   reused in the reciprocal-space energy computation to avoid recalculating
    !   trigonometric functions repeatedly.
    !
    ! For each k index (including negative k), we store the complex
    ! exponential. Exploiting the property:
    !    exp(-i * k * θ) = conjg(exp(i * k * θ))
    ! avoids redundant computations for negative k values.
    !--------------------------------------------------------------------
    pure subroutine ComputePhaseFactors1D(phase_factor_axis, phase_component, kmax)

        ! Input arguments
        complex(real64), intent(inout) :: phase_factor_axis(-kmax:kmax)  ! Array of complex phase factors for all k indices along one axis
        real(real64), intent(in) :: phase_component                  ! Phase angle for this atom along the current axis
        integer, intent(in) :: kmax                                  ! Maximum k-index in the positive direction
        ! Local arguments
        integer :: k                                                 ! Loop index over k-values
        complex(real64) :: phase_exp                                        ! Temporary complex exponential for current k

        do k = 0, kmax
            ! Compute complex exponential for positive k
            phase_exp = cmplx(dcos(k*phase_component), dsin(k*phase_component), kind=real64)
            phase_factor_axis(k) = phase_exp

            ! Exploit conjugate symmetry for negative k values
            if (k /= 0) then
                phase_factor_axis(-k) = conjg(phase_exp)
            end if
        end do

    end subroutine ComputePhaseFactors1D

    !------------------------------------------------------------------------------
    ! ComputeAllFourierTerms
    ! Computes the Fourier structure factors e^(i·k·r) for all atoms in all molecules.
    ! These are used in the reciprocal-space part of the Ewald summation to avoid
    ! redundant calculations of complex exponentials during the k-sum.
    !------------------------------------------------------------------------------
    subroutine ComputeAllFourierTerms()

        implicit none

        integer :: residue_type_1
        integer :: molecule_index_1

        ! Loop over all residue types
        do residue_type_1 = 1, nb%type_residue

            ! Loop over all molecules of this residue type
            do molecule_index_1 = 1, primary%num_residues(residue_type_1)

                ! Compute Fourier terms for a single molecule
                call SingleMolFourierTerms(residue_type_1, molecule_index_1)

            end do
        end do

    end subroutine ComputeAllFourierTerms

    !--------------------------------------------------------------------
    ! Computes the reciprocal-space contribution to the Coulomb energy
    !
    ! This subroutine calculates the Ewald reciprocal-space energy using
    ! precomputed k-vectors stored in ewald%kvectors. Each k-vector's
    ! structure factor amplitude is computed, squared, and weighted by
    ! the precomputed reciprocal constant and the form factor.
    !
    ! Advantages of using precomputed k-vectors:
    !   - Eliminates triple nested loops over kx, ky, kz indices.
    !   - Avoids repeated computation of k-vector magnitudes.
    !   - k=0 and vectors outside the unit-sphere cutoff are already excluded.
    !
    ! Energy formula for each k-vector:
    !   E_k = form_factor * recip_constant * |recip_amplitude|^2
    !
    ! The total reciprocal energy is scaled to eV per molecule.
    !--------------------------------------------------------------------
    subroutine ComputeReciprocalEnergy(u_recipCoulomb)

        implicit none

        ! Input arguments
        real(real64), intent(out) :: u_recipCoulomb
        ! Internal variables
        integer :: idx                   ! Index over precomputed reciprocal vectors
        real(real64) :: form_factor      ! Factor to account for symmetry (k vs -k)
        complex(real64) :: recip_amplitude    ! Structure factor for the current k-vector
        real(real64) :: recip_constant   ! Precomputed Ewald reciprocal-space weight
        real(real64) :: amplitude_sq     ! Squared modulus of the structure factor amplitude

        ! Initialize
        u_recipCoulomb = zero

        ! Loop over all precomputed reciprocal lattice vectors
        do idx = 1, ewald%num_kvectors

            ! Use pre-computed form factor
            form_factor = ewald%form_factor(idx)

            ! Compute the structure factor (complex amplitude) for this k-vector
            recip_amplitude = computeRecipAmplitude(ewald%kvectors(idx)%kx, &
                ewald%kvectors(idx)%ky, ewald%kvectors(idx)%kz)

            ! Retrieve the precomputed reciprocal-space weight
            recip_constant = ewald%recip_constants(idx)

            ! Compute squared modulus of the structure factor amplitude
            amplitude_sq = amplitude_squared(recip_amplitude)

            ! Accumulate reciprocal-space energy:
            ! form_factor * reciprocal constant * |amplitude|^2
            ! E_k = form_factor * recip_constant * |recip_amplitude|^2
            u_recipCoulomb = u_recipCoulomb + form_factor * recip_constant * amplitude_sq

        end do

        ! Convert accumulated energy to correct units (eV per molecule)
        u_recipCoulomb = u_recipCoulomb * EPS0_INV_eVA / KB_eVK * TWOPI / primary%volume

    end subroutine ComputeReciprocalEnergy

    !--------------------------------------------------------------------
    ! Computes the reciprocal-space amplitude for a given k-vector
    !
    ! The function sums over all atoms in all molecules of all residue types,
    ! multiplying the atom charge by the corresponding phase factors in x, y, z.
    !--------------------------------------------------------------------
    pure function ComputeRecipAmplitude(kx_idx, ky_idx, kz_idx) result(amplitude)

        implicit none

        ! Input arguments
        integer, intent(in) :: kx_idx, ky_idx, kz_idx
        ! Internal variables
        complex(real64) :: amplitude
        integer :: residue_type, molecule_index, atom_index

        ! Initialize amplitude to zero (complex)
        amplitude = (zero, zero)

        ! Loop over all residue types
        do residue_type = 1, nb%type_residue
            ! Loop over all molecules of this residue type
            do molecule_index = 1, primary%num_residues(residue_type)
                ! Loop over sites in molecule
                do atom_index = 1, nb%atom_in_residue(residue_type)
                    ! Accumulate contribution from this atom:
                    ! charge * exp(i k · r) factor in x, y, z directions
                    amplitude = amplitude + primary%atom_charges(residue_type, atom_index) * &
                                ewald%phase_factor_x(residue_type, molecule_index, atom_index, kx_idx) * &
                                ewald%phase_factor_y(residue_type, molecule_index, atom_index, ky_idx) * &
                                ewald%phase_factor_z(residue_type, molecule_index, atom_index, kz_idx)
                end do
            end do
        end do
    end function ComputeRecipAmplitude

    !--------------------------------------------------------------------
    ! SaveSingleMolFourierTerms
    !
    ! Saves the current phase factors and reciprocal amplitudes for a
    ! given molecule (or residue). This allows rollback of the Fourier
    ! state if a Monte Carlo or MD move is later rejected.
    !
    ! Steps:
    !   1. Save phase factors (IKX, IKY, IKZ) for each atom in the residue.
    !      - Handles both positive and negative indices for ky and kz.
    !   2. Save reciprocal amplitudes A(k) for all valid k-vectors
    !      into the backup array.
    !--------------------------------------------------------------------
    subroutine SaveSingleMolFourierTerms(residue_type, molecule_index)

        implicit none

        ! Input arguments
        integer, intent(in) :: residue_type    ! Residue type identifier
        integer, intent(in) :: molecule_index  ! Index of the molecule to save
        ! Local variables
        integer :: atom_index_1                ! Atom index within the residue
        integer :: kx_idx, ky_idx, kz_idx      ! Reciprocal vector indices
        integer :: idx                         ! Index over precomputed k-vectors

        do atom_index_1 = 1, nb%atom_in_residue(residue_type)

            ! Save IKX terms (kx_idx from 0 to kmax(1))
            do kx_idx = 0, ewald%kmax(1)
                ewald%phase_factor_x_old(atom_index_1, kx_idx) = &
                    ewald%phase_factor_x(residue_type, molecule_index, atom_index_1, kx_idx)
            end do

            ! Save IKY terms (KY from 0 to kmax(2)), include negative KY only if non-zero
            do ky_idx = 0, ewald%kmax(2)
                ewald%phase_factor_y_old(atom_index_1, ky_idx) = &
                    ewald%phase_factor_y(residue_type, molecule_index, atom_index_1, ky_idx)
                if (ky_idx /= 0) then
                    ewald%phase_factor_y_old(atom_index_1, -ky_idx) = &
                        ewald%phase_factor_y(residue_type, molecule_index, atom_index_1, -ky_idx)
                end if
            end do

            ! Save IKZ terms (KZ from 0 to kmax(3)), include negative KZ only if non-zero
            do kz_idx = 0, ewald%kmax(3)
                ewald%phase_factor_z_old(atom_index_1, kz_idx) = &
                    ewald%phase_factor_z(residue_type, molecule_index, atom_index_1, kz_idx)
                if (kz_idx /= 0) then
                    ewald%phase_factor_z_old(atom_index_1, -kz_idx) = &
                        ewald%phase_factor_z(residue_type, molecule_index, atom_index_1, -kz_idx)
                end if
            end do
        end do

        !------------------------------------------------------
        ! Step 2: Save reciprocal amplitudes A(k)
        !         Loop directly over precomputed valid k-vectors
        !------------------------------------------------------
        do idx = 1, ewald%num_kvectors
            ewald%recip_amplitude_old(idx) = ewald%recip_amplitude(idx)
        end do

    end subroutine SaveSingleMolFourierTerms

    !--------------------------------------------------------------------
    ! ComputeRecipEnergySingleMol
    !
    ! Purpose:
    !   Computes the reciprocal-space Coulomb energy contribution from a
    !   single molecule or residue using the Ewald summation method.
    !
    ! Description:
    !   For each precomputed reciprocal lattice vector (k-vector):
    !     1. Update the structure factor amplitude A(k) depending on the
    !        operation:
    !          - Creation:   A(k) ← A(k) + Σ q_i e^(i k·r_i,new)
    !          - Deletion:   A(k) ← A(k) - Σ q_i e^(i k·r_i,old)
    !          - Displacement/Update:
    !                        A(k) ← A(k) + Σ q_i [ e^(i k·r_i,new) - e^(i k·r_i,old) ]
    !
    !     2. Accumulate the reciprocal-space energy contribution:
    !            E_k = form_factor * W(k) * |A(k)|^2
    !        where:
    !          - form_factor accounts for symmetry (1 for kx=0, 2 otherwise),
    !          - W(k) is the precomputed reciprocal constant,
    !          - A(k) is the structure factor amplitude.
    !
    !   After processing all k-vectors, the total energy is scaled into
    !   physical units (eV per molecule).
    !
    ! Input arguments:
    !   residue_type   - index of the residue type
    !   molecule_index - index of the molecule
    !   is_creation    - optional flag, true if molecule is being created
    !   is_deletion    - optional flag, true if molecule is being deleted
    !
    ! Output arguments:
    !   u_recipCoulomb_new - reciprocal-space Coulomb energy (in eV)
    !
    ! Notes:
    !   - If neither is_creation nor is_deletion is present, the routine
    !     assumes a standard displacement/MC move and applies the
    !     difference (new - old) update.
    !   - Updates both the reciprocal amplitudes and the total reciprocal
    !     energy consistently.
    !--------------------------------------------------------------------
    subroutine ComputeRecipEnergySingleMol(residue_type, molecule_index, u_recipCoulomb_new, is_creation, is_deletion)

        implicit none

        ! Input arguments
        integer, intent(in) :: residue_type    ! Index of the residue type
        integer, intent(in) :: molecule_index  ! Index of the molecule in the system
        real(real64), intent(out) :: u_recipCoulomb_new  ! Output: reciprocal-space Coulomb energy
        logical, intent(in), optional :: is_creation
        logical, intent(in), optional :: is_deletion
        ! Local variables
        integer :: kx_idx, ky_idx, kz_idx      ! Components of current reciprocal lattice vector
        integer :: idx                         ! Loop index over precomputed k-vectors
        integer :: natoms                      ! Number of atoms in this residue type
        real(real64) :: amplitude_sq           ! Squared modulus of the structure factor amplitude
        real(real64) :: form_factor            ! Symmetry factor: 1 if kx=0, 2 otherwise
        real(real64), dimension(:), allocatable :: charges   ! Partial charges of atoms
        complex(real64), dimension(:), allocatable :: phase_new   ! Updated phase factor product
        complex(real64), dimension(:), allocatable :: phase_old   ! Previous phase factor product
        logical :: creation_flag
        logical :: deletion_flag

        creation_flag = present_or_false(is_creation)
        deletion_flag = present_or_false(is_deletion)

        ! Initialize energy accumulator
        u_recipCoulomb_new = zero

        ! Atom charges in this residue
        natoms = nb%atom_in_residue(residue_type)
        allocate(charges(natoms))
        charges = primary%atom_charges(residue_type, 1:natoms)

        ! Allocate phase_new and phase_old
        allocate(phase_new(natoms), phase_old(natoms))
        
        ! Loop over all precomputed reciprocal lattice vectors
        do idx = 1, ewald%num_kvectors

            ! Use pre-computed form factor
            form_factor = ewald%form_factor(idx)

            ! Current k-vector indices
            kx_idx = ewald%kvectors(idx)%kx
            ky_idx = ewald%kvectors(idx)%ky
            kz_idx = ewald%kvectors(idx)%kz

            ! Compute total phase factors (new vs. old configuration)
            phase_new = ewald%phase_factor_x(residue_type, molecule_index, 1:natoms, kx_idx) * &
                        ewald%phase_factor_y(residue_type, molecule_index, 1:natoms, ky_idx) * &
                        ewald%phase_factor_z(residue_type, molecule_index, 1:natoms, kz_idx)

            phase_old = ewald%phase_factor_x_old(1:natoms, kx_idx) * &
                        ewald%phase_factor_y_old(1:natoms, ky_idx) * &
                        ewald%phase_factor_z_old(1:natoms, kz_idx)

            ! Update Fourier coefficient A(k)
            if (creation_flag) then
                ! Molecule creation
                ! A(k) ← A(k) + Σ q_i [ e^(i k·r_i,new) ]
                ewald%recip_amplitude(idx) = ewald%recip_amplitude(idx) + sum(charges * phase_new)
            else if (deletion_flag) then
                ! Molecule deletion
                ! A(k) ← A(k) + Σ q_i [ - e^(i k·r_i,old) ]
                ewald%recip_amplitude(idx) = ewald%recip_amplitude(idx) - sum(charges * phase_old)
            else
                ! Standard move (translation, rotation)
                ! A(k) ← A(k) + Σ q_i [ e^(i k·r_i,new) - e^(i k·r_i,old) ]
                ewald%recip_amplitude(idx) = ewald%recip_amplitude(idx) + sum(charges * (phase_new - phase_old))
            end if

            ! Compute squared modulus of the structure factor amplitude
            amplitude_sq = amplitude_squared(ewald%recip_amplitude(idx))

            !----------------------------------------------
            ! Accumulate reciprocal-space energy:
            ! E_k = form_factor * W(k) * |A(k)|^2
            ! where W(k) is the precomputed reciprocal constant
            !----------------------------------------------
            u_recipCoulomb_new = u_recipCoulomb_new + form_factor * ewald%recip_constants(idx) * amplitude_sq
        end do

        !----------------------------------------------
        ! Convert accumulated energy to physical units:
        !----------------------------------------------
        u_recipCoulomb_new = u_recipCoulomb_new * EPS0_INV_eVA / KB_eVK * TWOPI / primary%volume

    end subroutine ComputeRecipEnergySingleMol

    !--------------------------------------------------------------------
    ! RestoreSingleMolFourier
    !
    ! Restores the phase factors and reciprocal amplitudes of a single
    ! molecule (or residue) from their saved "old" values. This is used
    ! to roll back Fourier-space state after a rejected MC move.
    !
    ! Steps:
    !   1. Restore phase factors (IKX, IKY, IKZ) for each atom in the residue.
    !      - Handles both positive and negative indices for ky and kz.
    !   2. Restore reciprocal amplitudes A(k) for all valid k-vectors
    !      from their backup array.
    !--------------------------------------------------------------------
    subroutine RestoreSingleMolFourier(residue_type, molecule_index)

        implicit none

        ! Input arguments
        integer, intent(in) :: residue_type    ! Residue type identifier
        integer, intent(in) :: molecule_index  ! Index of the molecule to restore
        ! Local variables
        integer :: atom_index_1                ! Atom index within the residue
        integer :: kx_idx, ky_idx, kz_idx      ! Reciprocal vector indices
        integer :: idx                         ! Index over precomputed k-vectors

        ! Restore per-atom phase factors (IKX, IKY, IKZ)
        do atom_index_1 = 1, nb%atom_in_residue(residue_type)

            ! Restore IKX
            do kx_idx = 0, ewald%kmax(1)
                ewald%phase_factor_x(residue_type, molecule_index, atom_index_1, kx_idx) = &
                    ewald%phase_factor_x_old(atom_index_1, kx_idx)
            end do

            ! Restore IKY (include negative KY only if non-zero)
            do ky_idx = 0, ewald%kmax(2)
                ewald%phase_factor_y(residue_type, molecule_index, atom_index_1, ky_idx) = &
                    ewald%phase_factor_y_old(atom_index_1, ky_idx)
                if (ky_idx /= 0) then
                    ewald%phase_factor_y(residue_type, molecule_index, atom_index_1, -ky_idx) = &
                        ewald%phase_factor_y_old(atom_index_1, -ky_idx)
                end if
            end do

            ! Restore IKZ (include negative KZ only if non-zero)
            do kz_idx = 0, ewald%kmax(3)
                ewald%phase_factor_z(residue_type, molecule_index, atom_index_1, kz_idx) = &
                    ewald%phase_factor_z_old(atom_index_1, kz_idx)
                if (kz_idx /= 0) then
                    ewald%phase_factor_z(residue_type, molecule_index, atom_index_1, -kz_idx) = &
                        ewald%phase_factor_z_old(atom_index_1, -kz_idx)
                end if
            end do
        end do

        !------------------------------------------------------
        ! Restore reciprocal amplitudes A(k) from backup
        ! Loop directly over precomputed valid k-vectors
        !------------------------------------------------------
        do idx = 1, ewald%num_kvectors
            ewald%recip_amplitude(idx) = ewald%recip_amplitude_old(idx)
        end do

    end subroutine RestoreSingleMolFourier

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
        
    !--------------------------------------------------------------------
    ! Returns squared modulus of a complex number:
    !   |z|^2 = Re(z)^2 + Im(z)^2
    !--------------------------------------------------------------------
    pure function amplitude_squared(z) result(val)

        ! Input argument
        complex(real64), intent(in) :: z
        ! Output rgument
        real(real64) :: val

        val = real(z*conjg(z), kind=real64)

    end function amplitude_squared

    !------------------------------------------------------------------------------
    ! ComputeEwaldSelfInteractionSingleMol
    ! Computes the Ewald self-energy correction for a single molecule.
    !------------------------------------------------------------------------------
    subroutine ComputeEwaldSelfInteractionSingleMol(residue_type, self_energy_1)
        
        implicit none

        ! Input arguments
        integer, intent(in) :: residue_type
        real(real64), intent(out) :: self_energy_1

        ! Local variables
        integer :: atom_index_1
        real(real64) :: charge_1

        self_energy_1 = 0.0_real64
        do atom_index_1 = 1, nb%atom_in_residue(residue_type)
            charge_1 = primary%atom_charges(residue_type, atom_index_1)
            if (abs(charge_1) < 1.0D-10) cycle
            self_energy_1 = self_energy_1 - ewald%alpha / SQRTPI * charge_1**2
        end do

        ! Scale by constants EPS0_INV_eVA and KB_eVK
        self_energy_1 = self_energy_1 * EPS0_INV_eVA / KB_eVK 

        return

    end subroutine ComputeEwaldSelfInteractionSingleMol

    !------------------------------------------------------------------------------
    ! ComputeIntraResidueRealCoulombEnergySingleMol
    ! Computes real-space intramolecular Coulomb energy for a molecule.
    !------------------------------------------------------------------------------
    subroutine ComputeIntraResidueRealCoulombEnergySingleMol(residue_type, molecule_index, u_intraCoulomb_1)

        implicit none

        ! Input arguments
        integer, intent(in) :: residue_type
        integer, intent(in) :: molecule_index
        real(real64), intent(out) :: u_intraCoulomb_1

        ! Local variables
        integer :: atom_index_1, atom_index_2
        real(real64) :: distance
        real(real64) :: charge_1, charge_2

        u_intraCoulomb_1 = 0

        do atom_index_1 = 1, nb%atom_in_residue(residue_type)-1
            charge_1 = primary%atom_charges(residue_type, atom_index_1)

            do atom_index_2 = atom_index_1+1, nb%atom_in_residue(residue_type)
                charge_2 = primary%atom_charges(residue_type, atom_index_2)

                distance = ComputeDistance(primary, residue_type, molecule_index, atom_index_1, &
                            residue_type, molecule_index, atom_index_2)

                if (distance > 1.0d-10) then
                    u_intraCoulomb_1 = u_intraCoulomb_1 + &
                        charge_1 * charge_2 * (erfc(ewald%alpha * distance) - 1.0d0) / distance
                end if

            end do
        end do

        ! Scale by constants EPS0_INV_eVA and KB_eVK
        u_intraCoulomb_1 = u_intraCoulomb_1 * EPS0_INV_eVA / KB_eVK

    end subroutine ComputeIntraResidueRealCoulombEnergySingleMol

    ! Precompute complex exponential factors to be used repeatedly in
    ! reciprocal-space calculations.
    subroutine SingleMolFourierTerms(residue_type, molecule_index)

        implicit none

        ! Input arguments
        integer, intent(in) :: residue_type
        integer, intent(in) :: molecule_index

        ! Local variables
        integer :: atom_index_1                 ! Atom index
        real(real64), dimension(3) :: atom      ! Atom coordinates in real space
        real(real64), dimension(3) :: phase     ! Phase factors for Fourier terms
        integer :: kmax_x, kmax_y, kmax_z       ! max k indices for each direction
        complex(real64), allocatable :: temp_x(:), temp_y(:), temp_z(:)  ! temp arrays for phase factors

        ! Determine kmax for each direction
        kmax_x = ewald%kmax(1)
        kmax_y = ewald%kmax(2)
        kmax_z = ewald%kmax(3)

        ! Allocate temporary arrays once
        allocate(temp_x(-kmax_x:kmax_x))
        allocate(temp_y(-kmax_y:kmax_y))
        allocate(temp_z(-kmax_z:kmax_z))

        do atom_index_1 = 1, nb%atom_in_residue(residue_type)

            atom = primary%mol_com(:, residue_type, molecule_index) + &
                primary%site_offset(:, residue_type, molecule_index, atom_index_1)

            ! Compute the phase vector components as the dot product of the atom position
            ! with each reciprocal lattice vector (columns of reciprocal_box), scaled by 2π.
            call ComputeAtomPhase(atom, primary%reciprocal, phase)

            ! Precompute the complex exponential (phase) factors for this atom
            ! along each Cartesian direction. These factors will be used repeatedly
            ! in the reciprocal-space sum for the Ewald energy.
            call ComputePhaseFactors1D(temp_x, phase(1), kmax_x)
            ewald%phase_factor_x(residue_type, molecule_index, atom_index_1, :) = temp_x

            call ComputePhaseFactors1D(temp_y, phase(2), kmax_y)
            ewald%phase_factor_y(residue_type, molecule_index, atom_index_1, :) = temp_y

            call ComputePhaseFactors1D(temp_z, phase(3), kmax_z)
            ewald%phase_factor_z(residue_type, molecule_index, atom_index_1, :) = temp_z

        end do

    end subroutine SingleMolFourierTerms

    subroutine ReplaceFourierTerms_singlemol(residue_type, index_1, index_2)

        implicit none

        ! Input arguments
        integer, intent(in) :: residue_type
        integer, intent(in) :: index_1
        integer, intent(in) :: index_2

        ! Local variables
        integer :: kx_idx, ky_idx, kz_idx
        integer :: atom_index_1
        integer :: number_of_kpoints

        ! Initialize
        number_of_kpoints = 0

        ! Restore IKX, IKY, IKZ from backups
        do atom_index_1 = 1, nb%atom_in_residue(residue_type)

            ! Restore IKX
            do kx_idx = 0, ewald%kmax(1)
                ewald%phase_factor_x(residue_type, index_1, atom_index_1, kx_idx) = &
                    ewald%phase_factor_x(residue_type, index_2, atom_index_1, kx_idx)
            end do

            ! Restore IKY (include negative KY only if non-zero)
            do ky_idx = 0, ewald%kmax(2)
                ewald%phase_factor_y(residue_type, index_1, atom_index_1, ky_idx) = &
                    ewald%phase_factor_y(residue_type, index_2, atom_index_1, ky_idx)
                if (ky_idx /= 0) then
                    ewald%phase_factor_y(residue_type, index_1, atom_index_1, -ky_idx) = &
                        ewald%phase_factor_y(residue_type, index_2, atom_index_1, -ky_idx)
                end if
            end do

            ! Restore IKZ (include negative KZ only if non-zero)
            do kz_idx = 0, ewald%kmax(3)
                ewald%phase_factor_z(residue_type, index_1, atom_index_1, kz_idx) = &
                    ewald%phase_factor_z(residue_type, index_2, atom_index_1, kz_idx)
                if (kz_idx /= 0) then
                    ewald%phase_factor_z(residue_type, index_1, atom_index_1, -kz_idx) = &
                        ewald%phase_factor_z(residue_type, index_2, atom_index_1, -kz_idx)
                end if
            end do
        end do

    end subroutine ReplaceFourierTerms_singlemol

end module ewald_utils
