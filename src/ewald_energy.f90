module ewald_energy

    use simulation_state
    use output_utils
    use constants
    use geometry_utils
    use, intrinsic :: iso_fortran_env, only: real64

    implicit none

contains

    !--------------------------------------------------------------------
    ! ComputeRecipAmplitude
    !
    ! Purpose:
    !   Computes the reciprocal-space structure factor amplitude A(k) for
    !   a specified reciprocal lattice vector (kx, ky, kz).
    !
    ! Description:
    !   The structure factor is obtained by summing over all atoms in all
    !   molecules of all residue types:
    !
    !       A(k) = Σ q_i * exp(i k · r_i)
    !
    !   where q_i is the atomic charge and the exponential is represented
    !   by precomputed phase factors in x, y, and z directions.
    !
    ! Input:
    !   kx_idx, ky_idx, kz_idx  - indices of the reciprocal lattice vector
    !
    ! Output:
    !   amplitude (complex)     - total structure factor amplitude A(k)
    !
    ! Notes:
    !   - The returned amplitude is dimensionless and used later in the
    !     reciprocal-space energy calculation.
    !   - All atoms in the system contribute to the sum.
    !--------------------------------------------------------------------
    pure function ComputeRecipAmplitude(kx_idx, ky_idx, kz_idx) result(amplitude)

        implicit none

        ! Input arguments
        integer, intent(in) :: kx_idx, ky_idx, kz_idx       ! Reciprocal lattice vector indices
        ! Internal variables
        complex(real64) :: amplitude                        ! Accumulated structure factor amplitude A(k)
        complex(real64) :: phase                            ! Phase factor product for a single atom
        real(real64) :: charges                             ! Partial charge of the current atom
        integer :: residue_type                             ! Index of the current residue type
        integer :: molecule_index                           ! Index of the current molecule
        integer :: atom_index                               ! Index of the current atom

        ! Initialize amplitude to zero (complex)
        amplitude = (zero, zero)

        ! Loop over all residue types
        do residue_type = 1, nb%type_residue
            ! Loop over all molecules of this residue type
            do molecule_index = 1, primary%num_residues(residue_type)
                ! Loop over sites in molecule
                do atom_index = 1, nb%atom_in_residue(residue_type)

                    ! Extract charge and phase
                    charges = primary%atom_charges(residue_type, atom_index)
                    phase = ewald%phase_factor_x(residue_type, molecule_index, atom_index, kx_idx) * &
                            ewald%phase_factor_y(residue_type, molecule_index, atom_index, ky_idx) * &
                            ewald%phase_factor_z(residue_type, molecule_index, atom_index, kz_idx)

                    ! Accumulate contribution from this atom:
                    ! charge * exp(i k · r) factor in x, y, z directions
                    amplitude = amplitude + charges * phase

                end do
            end do
        end do
    end function ComputeRecipAmplitude

    !--------------------------------------------------------------------
    ! ComputeReciprocalEnergy
    !
    ! Purpose:
    !   Computes the total reciprocal-space Coulomb energy for the system.
    !
    ! Description:
    !   Loops over all precomputed k-vectors and evaluates:
    !
    !       E_k = form_factor * W(k) * |A(k)|^2
    !
    !   where:
    !     - A(k) is the structure factor amplitude
    !     - W(k) is the reciprocal-space weight
    !     - form_factor accounts for k vs -k symmetry
    !
    ! Input:
    !   None
    !
    ! Output:
    !   u_recipCoulomb : total reciprocal Coulomb energy (eV per molecule)
    !
    ! Notes:
    !   - Uses precomputed reciprocal vectors and phase factors.
    !   - Converts the final sum to physical units at the end.
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
    ! ComputeEwaldSelfInteractionSingleMol
    !
    ! Purpose:
    !   Computes the self-interaction correction for the reciprocal-space
    !   Coulomb energy of a single molecule/residue in the Ewald summation.
    !
    ! Background:
    !   In Ewald summation, each charge interacts with all periodic images
    !   of the system. This includes a spurious interaction of a charge
    !   with itself, which must be subtracted to obtain the correct energy.
    !
    !   The self-energy for a single charge q_i is:
    !
    !       E_self(i) = - (α / √π) * q_i^2
    !
    !   where α is the Ewald screening parameter. The total self-energy
    !   for a molecule/residue is the sum over all atoms in that molecule:
    !
    !       E_self = Σ_i E_self(i) = - (α / √π) * Σ_i q_i^2
    !
    ! Inputs:
    !   residue_type : integer index identifying the molecule/residue
    !
    ! Outputs:
    !   self_energy_1 : real(real64)
    !       Reciprocal-space self-energy correction for this molecule (in eV)
    !
    ! Notes:
    !   - Charges below 1e-10 are ignored to prevent unnecessary computations.
    !   - The result is scaled to eV per molecule using EPS0_INV_eVA / KB_eVK.
    !--------------------------------------------------------------------
    subroutine ComputeEwaldSelfInteractionSingleMol(residue_type, self_energy)
        
        implicit none

        ! Input arguments
        integer, intent(in) :: residue_type
        real(real64), intent(out) :: self_energy
        ! Local variables
        integer :: atom_index_1
        real(real64) :: charge_1

        ! Initialize self-energy accumulator
        self_energy = zero

        ! Loop over all atoms in the molecule/residue
        do atom_index_1 = 1, nb%atom_in_residue(residue_type)
            charge_1 = primary%atom_charges(residue_type, atom_index_1)

            ! Skip atoms with negligible charge
            if (abs(charge_1) < error) cycle

            ! Add the individual atomic self-energy contribution
            self_energy = self_energy - ewald%alpha / SQRTPI * charge_1**2
        end do

        ! Scale by constants EPS0_INV_eVA and KB_eVK
        self_energy = self_energy * EPS0_INV_eVA / KB_eVK 

    end subroutine ComputeEwaldSelfInteractionSingleMol

    !------------------------------------------------------------------------------
    ! ComputeIntraResidueRealCoulombEnergySingleMol
    !
    ! Purpose:
    !   Computes the **real-space intramolecular Coulomb energy** for a single
    !   molecule or residue using the Ewald summation method.
    !
    ! Background:
    !   In Ewald summation, the real-space contribution for a pair of charges
    !   q_i and q_j separated by distance r_ij is:
    !
    !       E_ij_real = q_i * q_j * (erfc(α * r_ij) / r_ij)
    !
    !   For intramolecular interactions, a correction term (-q_i*q_j/r_ij) is
    !   sometimes included to avoid double counting or to remove the singularity
    !   at r -> 0. This implementation computes:
    !
    !       E_intra = Σ_{i<j} q_i * q_j * (erfc(α * r_ij) - 1) / r_ij
    !
    ! Inputs:
    !   residue_type   : integer
    !       Index identifying the residue type / molecule
    !   molecule_index : integer
    !       Index of the molecule in the system
    !
    ! Outputs:
    !   u_intraCoulomb_1 : real(real64)
    !       Real-space intramolecular Coulomb energy (in eV)
    !
    ! Notes:
    !   - Pairs of atoms with distance < 1e-10 are skipped to avoid singularity.
    !   - The result is scaled to physical units using EPS0_INV_eVA / KB_eVK.
    !------------------------------------------------------------------------------
    subroutine ComputeIntraResidueRealCoulombEnergySingleMol(residue_type, molecule_index, u_intraCoulomb)

        implicit none

        ! Input arguments
        integer, intent(in) :: residue_type
        integer, intent(in) :: molecule_index
        real(real64), intent(out) :: u_intraCoulomb
        ! Local variables
        integer :: atom_index_1, atom_index_2
        real(real64) :: distance
        real(real64) :: charge_1, charge_2

        ! Initialize energy accumulator
        u_intraCoulomb = 0

        ! Loop over all unique atom pairs in the molecule
        do atom_index_1 = 1, nb%atom_in_residue(residue_type)-1
            charge_1 = primary%atom_charges(residue_type, atom_index_1)

            do atom_index_2 = atom_index_1+1, nb%atom_in_residue(residue_type)
                charge_2 = primary%atom_charges(residue_type, atom_index_2)

                ! Compute interatomic distance
                distance = ComputeDistance(primary, residue_type, molecule_index, atom_index_1, &
                            residue_type, molecule_index, atom_index_2)

                ! Skip extremely small distances to avoid singularity
                if (distance > error) then
                    ! Add the real-space intramolecular contribution:
                    ! (erfc(α * r) - 1) / r
                    u_intraCoulomb = u_intraCoulomb + charge_1 * charge_2 * (erfc(ewald%alpha * distance) - one) / distance
                end if

            end do
        end do

        ! Scale by constants EPS0_INV_eVA and KB_eVK
        u_intraCoulomb = u_intraCoulomb * EPS0_INV_eVA / KB_eVK

    end subroutine ComputeIntraResidueRealCoulombEnergySingleMol

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

end module ewald_energy
