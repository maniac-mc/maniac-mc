module ewald_phase

    use simulation_state
    use output_utils
    use constants
    use geometry_utils
    use, intrinsic :: iso_fortran_env, only: real64

    implicit none

contains

    !--------------------------------------------------------------------
    ! ComputeAtomPhase
    !
    ! Purpose:
    !   Computes the phase vector (k · r) for a single atom in reciprocal space.
    !
    ! Mathematical Formulation:
    !   For each reciprocal direction i ∈ {1,2,3}:
    !
    !       phase(i) = 2π * Σ_j [ reciprocal_box(j,i) * atom_pos(j) ]
    !
    !   or in vector form:
    !
    !       phase = 2π * (reciprocal_boxᵀ · atom_pos)
    !
    ! Arguments:
    !   atom_pos        [in]  : Real(3) position of the atom in real space.
    !   reciprocal_box  [in]  : Real(3,3) reciprocal lattice vectors (columns).
    !
    ! Returns:
    !   phase(3) : Real(3) phase vector (in radians) used in exp(i * k · r)
    !
    ! Notes:
    !   This produces the phase argument (in radians) for the exponential term
    !   exp(i * k · r) used in the Ewald summation.
    !   This function is pure and side-effect-free. It can be safely used in
    !   parallel sections or array expressions.
    !--------------------------------------------------------------------
    pure function ComputeAtomPhase(atom_pos, reciprocal_box) result(phase)

        ! Input arguments
        real(real64), intent(in)  :: atom_pos(3)         ! Cartesian position of the atom in real space (x, y, z)
        real(real64), intent(in)  :: reciprocal_box(3,3) ! Reciprocal lattice vectors as columns of a 3×3 matrix

        ! Function result
        real(real64) :: phase(3)                        ! Returned phase components (2π * k·r)

        ! Local variables
        integer :: i, j                                 ! Loop indices over spatial dimensions

        do i = 1, 3
            phase(i) = zero
            do j = 1, 3
                ! Accumulate dot product of atom position (r_j)
                ! with i-th reciprocal lattice vector (k_i = column i of reciprocal_box)
                phase(i) = phase(i) + reciprocal_box(j,i) * atom_pos(j)
            end do
            ! Multiply by 2π to convert to the proper phase angle (k·r in radians)
            phase(i) = TWOPI * phase(i)
        end do

    end function ComputeAtomPhase

    !--------------------------------------------------------------------
    ! ComputePhaseFactors1D
    !
    ! Purpose:
    !   Computes 1D complex exponential phase factors exp(i * k * θ)
    !   along a single reciprocal axis for one atom.
    !
    ! Mathematical Formulation:
    !   For k ∈ [-kmax, kmax]:
    !       phase_factor_axis(k) = exp(i * k * phase_component)
    !
    ! Arguments:
    !   phase_factor_axis [inout] : Complex array (-kmax:kmax) storing computed
    !                               phase factors for all integer k values.
    !   phase_component   [in]    : Real phase angle (θ) for the atom along this axis.
    !   kmax              [in]    : Maximum positive reciprocal index (defines array bounds).
    !
    ! Notes:
    !   Exploits conjugate symmetry:
    !       exp(-i * k * θ) = conjg(exp(i * k * θ))
    !   to avoid redundant trigonometric calculations for negative k-values.
    !   The resulting array is reused during the reciprocal-space energy evaluation
    !   to minimize repeated sine/cosine computations.
    !--------------------------------------------------------------------
    pure subroutine ComputePhaseFactors1D(phase_factor_axis, phase_component, kmax)

        ! Input arguments
        complex(real64), intent(inout) :: phase_factor_axis(-kmax:kmax) ! Array of complex phase factors for all k indices along one axis
        real(real64), intent(in) :: phase_component                  ! Phase angle for this atom along the current axis
        integer, intent(in) :: kmax                                  ! Maximum k-index in the positive direction
        ! Local arguments
        integer :: k                                                 ! Loop index over k-values
        complex(real64) :: phase_exp                                 ! Temporary complex exponential for current k

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

    !--------------------------------------------------------------------
    ! SaveSingleMolFourierTerms
    !
    ! Purpose:
    !   Saves the current Fourier-space phase factors and reciprocal amplitudes
    !   for a specific molecule or residue. This enables rollback after a
    !   rejected Monte Carlo or molecular dynamics move.
    !
    ! Arguments:
    !   residue_type   [in] : Integer residue type identifier.
    !   molecule_index [in] : Integer molecule index within this residue type.
    !
    ! Description:
    !   1. Saves per-atom 1D phase factor arrays (IKX, IKY, IKZ) for all k-indices.
    !      Handles both positive and negative k-values where applicable.
    !   2. Saves all reciprocal amplitudes A(k) into backup arrays.
    !
    ! Notes:
    !   The backup arrays (`_old`) preserve the previous Fourier state, allowing
    !   full restoration without recomputation if a trial move is rejected.
    !--------------------------------------------------------------------
    subroutine SaveSingleMolFourierTerms(residue_type, molecule_index)

        implicit none

        ! Input arguments
        integer, intent(in) :: residue_type    ! Residue type identifier
        integer, intent(in) :: molecule_index  ! Index of the molecule to save
        ! Local variables
        integer :: atom_index                  ! Atom index within the residue
        integer :: kx_idx, ky_idx, kz_idx      ! Reciprocal vector indices
        integer :: idx                         ! Index over precomputed k-vectors

        do atom_index = 1, nb%atom_in_residue(residue_type)

            ! Save IKX terms (kx_idx from 0 to kmax(1))
            do kx_idx = 0, ewald%kmax(1)
                ewald%phase_factor_x_old(atom_index, kx_idx) = &
                    ewald%phase_factor_x(residue_type, molecule_index, atom_index, kx_idx)
            end do

            ! Save IKY terms (KY from 0 to kmax(2)), include negative KY only if non-zero
            do ky_idx = 0, ewald%kmax(2)
                ewald%phase_factor_y_old(atom_index, ky_idx) = &
                    ewald%phase_factor_y(residue_type, molecule_index, atom_index, ky_idx)
                if (ky_idx /= 0) then
                    ewald%phase_factor_y_old(atom_index, -ky_idx) = &
                        ewald%phase_factor_y(residue_type, molecule_index, atom_index, -ky_idx)
                end if
            end do

            ! Save IKZ terms (KZ from 0 to kmax(3)), include negative KZ only if non-zero
            do kz_idx = 0, ewald%kmax(3)
                ewald%phase_factor_z_old(atom_index, kz_idx) = &
                    ewald%phase_factor_z(residue_type, molecule_index, atom_index, kz_idx)
                if (kz_idx /= 0) then
                    ewald%phase_factor_z_old(atom_index, -kz_idx) = &
                        ewald%phase_factor_z(residue_type, molecule_index, atom_index, -kz_idx)
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
    ! RestoreSingleMolFourier
    !
    ! Purpose:
    !   Restores previously saved Fourier-space phase factors and reciprocal
    !   amplitudes for a molecule or residue after a rejected move.
    !
    ! Arguments:
    !   residue_type   [in] : Integer residue type identifier.
    !   molecule_index [in] : Integer molecule index within this residue type.
    !
    ! Description:
    !   1. Restores all 1D phase factor arrays (IKX, IKY, IKZ) from backups.
    !      Includes negative k-values for ky and kz when nonzero.
    !   2. Restores reciprocal amplitudes A(k) from backup arrays.
    !
    ! Notes:
    !   Ensures that the reciprocal-space data structures are consistent with
    !   the accepted configuration before the rejected move attempt.
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

    !--------------------------------------------------------------------
    ! ReplaceFourierTermsSingleMol
    !
    ! Purpose:
    !   Replaces the Fourier-space phase factors of one molecule (`index_1`)
    !   with those from another molecule (`index_2`) of the same residue type.
    !   Used when molecule identities or positions are swapped (e.g., in
    !   exchange Monte Carlo or symmetry operations).
    !
    ! Arguments:
    !   residue_type [in] : Integer residue type identifier.
    !   index_1      [in] : Destination molecule index (to be updated).
    !   index_2      [in] : Source molecule index (providing data).
    !
    ! Notes:
    !   Copies 1D phase factor arrays (IKX, IKY, IKZ) for all atoms and all
    !   reciprocal indices, including negative k-values where nonzero.
    !   Reciprocal amplitudes are not modified by this routine.
    !--------------------------------------------------------------------
    subroutine ReplaceFourierTermsSingleMol(residue_type, index_1, index_2)

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

    end subroutine ReplaceFourierTermsSingleMol

    !--------------------------------------------------------------------
    ! ComputeAllFourierTerms
    !
    ! Purpose:
    !   Computes and stores the Fourier structure factors exp(i * k · r)
    !   for all atoms in all molecules across all residue types.
    !
    ! Description:
    !   Loops over every residue type and every molecule within that type,
    !   invoking SingleMolFourierTerms to precompute and cache phase factors
    !   used in the reciprocal-space portion of the Ewald summation.
    !
    ! Notes:
    !   This routine ensures that all required exp(i * k · r) terms are
    !   available before evaluating reciprocal-space energies or forces.
    !--------------------------------------------------------------------
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
    ! SingleMolFourierTerms
    !
    ! Purpose:
    !   Computes the per-atom Fourier-space phase factors exp(i * k · r)
    !   for all k-indices along each Cartesian direction for a given molecule.
    !
    ! Arguments:
    !   residue_type   [in] : Integer residue type identifier.
    !   molecule_index [in] : Integer molecule index within this residue type.
    !
    ! Description:
    !   For each atom in the molecule:
    !     1. Compute its real-space position.
    !     2. Compute its phase vector: phase = 2π * (reciprocal_boxᵀ · r_atom).
    !     3. Precompute exp(i * k * phase_component) arrays for x, y, and z axes.
    !
    ! Notes:
    !   These precomputed 1D phase factors are used to construct the full
    !   reciprocal-space structure factors efficiently during the Ewald sum.
    !--------------------------------------------------------------------
    subroutine SingleMolFourierTerms(residue_type, molecule_index)

        implicit none

        ! Input arguments
        integer, intent(in) :: residue_type
        integer, intent(in) :: molecule_index

        ! Local variables
        integer :: atom_index_1                 ! Atom index
        real(real64), dimension(3) :: atom      ! Atom coordinates in real space
        real(real64), dimension(3) :: phase     ! Phase factors for Fourier terms
        integer :: kmax_x, kmax_y, kmax_z       ! Max k indices for each direction
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
            ! Compute the phase vector (2π * reciprocal_boxᵀ · atom)
            phase = ComputeAtomPhase(atom, primary%reciprocal)

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

end module ewald_phase
