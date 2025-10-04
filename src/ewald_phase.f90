module ewald_phase

    use simulation_state
    use output_utils
    use constants
    use geometry_utils
    use, intrinsic :: iso_fortran_env, only: real64

    implicit none

contains

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

end module ewald_phase
