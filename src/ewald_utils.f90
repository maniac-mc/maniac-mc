module ewald_utils

    use simulation_state
    use output_utils 
    use constants
    use geometry_utils
    use, intrinsic :: iso_fortran_env, only: real64

    implicit none

contains

    !------------------------------------------------------------------------------
    ! Precomputes the reciprocal-space weighting factors for the Ewald summation.
    ! Uses a 3x3 matrix to store the reciprocal lattice vectors and computes the
    ! magnitude of each k-vector efficiently. Only done once at the start of the simulation.
    !------------------------------------------------------------------------------
    subroutine InitializeReciprocalWeights()

        implicit none

        ! Local variables
        real(real64) :: k_squared      ! Normalized squared k-vector for unit-sphere check
        real(real64) :: k_squared_mag  ! Squared magnitude of the k-vector
        integer :: number_of_kpoints   ! Counter for valid k-vectors
        integer :: kx_idx, ky_idx, kz_idx  ! Reciprocal lattice indices
        real(real64) :: kvec(3)        ! 3D reciprocal lattice vector
        real(real64) :: kidx_scaled(3) ! Scaled k-vector indices for unit-sphere filter
        real(real64), dimension(3,3) :: kvec_matrix  ! Columns are reciprocal lattice vectors b1, b2, b3

        ! Store reciprocal lattice vectors as columns of a 3x3 matrix
        kvec_matrix = TWOPI * reshape(primary%reciprocal, shape(kvec_matrix))

        ! Computes a weighting factor for each reciprocal lattice vector
        ! Only needs to be done once at the start of the simulations
        number_of_kpoints = 0

        ! Loop over all reciprocal lattice indices within kmax limits
        do kx_idx = 0, ewald%kmax(1)
            do ky_idx = -ewald%kmax(2), ewald%kmax(2) 
                do kz_idx = -ewald%kmax(3), ewald%kmax(3) 

                    ! Skip k=0 vector
                    if (kx_idx == 0 .and. ky_idx == 0 .and. kz_idx == 0) cycle

                    ! Compute the k-vector
                    kvec = dble(kx_idx)*kvec_matrix(:,1) + &
                        dble(ky_idx)*kvec_matrix(:,2) + &
                        dble(kz_idx)*kvec_matrix(:,3)

                    ! Normalized k^2 to filter vectors outside the unit sphere
                    kidx_scaled = [dble(kx_idx)/dble(ewald%kmax(1)), &
                                dble(ky_idx)/dble(ewald%kmax(2)), &
                                dble(kz_idx)/dble(ewald%kmax(3))]

                    k_squared = sum(kidx_scaled**2)

                    ! Skip k-vector if it lies outside the unit sphere in normalized space
                    if (k_squared > one) cycle

                    ! Weighting factor for this k-vector
                    k_squared_mag = dot_product(kvec, kvec)

                    ! Weight for this k-vector: damps short-range contributions with Gaussian and scales by 1/k^2
                    number_of_kpoints = number_of_kpoints + 1
                    ewald%recip_constants(number_of_kpoints) = exp(-k_squared_mag/(four*ewald%alpha**2)) / k_squared_mag

                end do
            end do
        end do

        ! Sanity check
        if (number_of_kpoints > ewald%num_recip_vectors) then
            call AbortRun("Too many reciprocal vectors")
        end if

    end subroutine InitializeReciprocalWeights

    ! Precompute complex exponential factors to be used repeatedly in
    ! reciprocal-space calculations.
    subroutine SingleMolFourierTerms(residue_type, molecule_index)

        implicit none

        ! Input arguments
        integer, intent(in) :: residue_type
        integer, intent(in) :: molecule_index

        ! Local variables
        integer :: atom_index_1
        integer :: kx_idx, ky_idx, kz_idx
        real(real64), dimension(3) :: atom ! Atom coordinates in real space
        real(real64), dimension(3) :: phase ! Phase factors for Fourier terms
        integer :: i, j ! Integer for looping

        do atom_index_1 = 1, nb%atom_in_residue(residue_type)

            atom = primary%mol_com(:, residue_type, molecule_index) + &
                primary%site_offset(:, residue_type, molecule_index, atom_index_1)

            ! Compute the phase vector components as the dot product of the atom position
            ! with each reciprocal lattice vector (columns of reciprocal_box), scaled by 2π.
            do i = 1, 3
                phase(i) = 0.0d0
                do j = 1, 3
                    phase(i) = phase(i) + primary%reciprocal(j, i) * atom(j)
                end do
                phase(i) = TWOPI * phase(i)
            end do

            do kx_idx = 0, ewald%kmax(1)
                ewald%phase_factor_x(residue_type, molecule_index, atom_index_1, kx_idx) = &
                    cmplx(dcos(kx_idx*phase(1)), dsin(kx_idx*phase(1)), KIND=8)
            end do

            do ky_idx = 0, ewald%kmax(2)
                ewald%phase_factor_y(residue_type, molecule_index, atom_index_1, ky_idx)  = &
                    cmplx(dcos(ky_idx*phase(2)), dsin(ky_idx*phase(2)), KIND=8)
                if (ky_idx /= 0) then
                    ewald%phase_factor_y(residue_type, molecule_index, atom_index_1, -ky_idx) = &
                        cmplx(dcos(-ky_idx*phase(2)), dsin(-ky_idx*phase(2)), KIND=8)
                end if
            end do

            do kz_idx = 0, ewald%kmax(3)
                ewald%phase_factor_z(residue_type, molecule_index, atom_index_1, kz_idx) = &
                    cmplx(dcos(kz_idx*phase(3)), dsin(kz_idx*phase(3)), KIND=8)
                if (kz_idx /= 0) then
                    ewald%phase_factor_z(residue_type, molecule_index, atom_index_1, -kz_idx) = &
                        cmplx(dcos(-kz_idx*phase(3)), dsin(-kz_idx*phase(3)), KIND=8)
                end if
            end do
        end do

    end subroutine SingleMolFourierTerms

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
    ! The energy is computed using Ewald summation in reciprocal space by
    ! looping over allowed k-vectors (excluding zero and vectors outside
    ! the unit sphere in normalized k-space) and summing the squared
    ! structure factor amplitudes multiplied by the reciprocal constants.
    !--------------------------------------------------------------------
    subroutine ComputeReciprocalEnergy(u_recipCoulomb)
        implicit none

        ! Input arguments
        real(real64), intent(out) :: u_recipCoulomb

        ! Internal variables
        integer :: number_of_kpoints
        integer :: kx_idx, ky_idx, kz_idx          ! Index for reciprocal lattice vector x-y-z-component
        integer :: atom_index_1, molecule_index, residue_type
        real(real64) :: form_factor, k_squared

        ! Initialize
        number_of_kpoints = 0
        u_recipCoulomb = 0.0

        ! Compute reciprocal-space energy sum over allowed k-vectors
        ! excluding zero and outside cutoff
        do kx_idx = 0, ewald%kmax(1)

            ! Compute form factor for current kx index: 1 for zero, 2 otherwise
            form_factor =  FormFactor(kx_idx)

            do ky_idx = -ewald%kmax(2), ewald%kmax(2)
                do kz_idx = -ewald%kmax(3), ewald%kmax(3)

                    ! Compute normalized squared magnitude of k-vector
                    k_squared = (dble(kx_idx)/dble(ewald%kmax(1)))**2 + &
                                (dble(ky_idx)/dble(ewald%kmax(2)))**2 + &
                                (dble(kz_idx)/dble(ewald%kmax(3)))**2

                    ! Skip k=0 vector and vectors outside the unit sphere in normalized space
                    if (abs(k_squared) < error .or. k_squared > one) cycle

                    ! Increment the number of valid k-points
                    number_of_kpoints = number_of_kpoints + 1

                    ! Compute the complex reciprocal amplitude (structure factor) for this k-vector
                    ewald%recip_amplitude(number_of_kpoints) = computeRecipAmplitude(kx_idx, ky_idx, kz_idx)

                    ! Accumulate reciprocal-space energy:
                    ! form_factor * reciprocal constant * |amplitude|^2
                    u_recipCoulomb = u_recipCoulomb + form_factor * ewald%recip_constants(number_of_kpoints) * &
                                    real(ewald%recip_amplitude(number_of_kpoints) * &
                                    conjg(ewald%recip_amplitude(number_of_kpoints)), KIND=8)

                end do
            end do
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

    ! Returns form factor: 1 if index is zero, 2 otherwise
    pure function FormFactor(idx) result(factor)
        integer, intent(in) :: idx
        real(real64) :: factor

        if (idx == 0) then
            factor = 1.0_real64
        else
            factor = 2.0_real64
        end if
    end function FormFactor







    ! subroutine ComputeReciprocalEnergy(u_recipCoulomb)

    !     implicit none

    !     ! Input arguments
    !     real(real64), intent(out) :: u_recipCoulomb

    !     ! Internal variables
    !     integer :: number_of_kpoints
    !     integer :: kx_idx, ky_idx, kz_idx ! Index for reciprocal lattice vector x-y-z-component
    !     integer :: atom_index_1, molecule_index, residue_type
    !     real(real64) :: form_factor, k_squared

    !     ! Compute reciprocal-space energy sum over allowed k-vectors excluding zero and outside cutoff
    !     number_of_kpoints = 0
    !     u_recipCoulomb = 0
    !     do kx_idx = 0, ewald%kmax(1)

    !         ! ! Weight for positive/negative k
    !         ! form_factor = merge(two, one, kx_idx==0)

    !         if (kx_idx == 0) then
    !             form_factor = one
    !         else
    !             form_factor = two
    !         end if

    !         do ky_idx = -ewald%kmax(2), ewald%kmax(2)
    !             do kz_idx = -ewald%kmax(3), ewald%kmax(3)

    !                 ! Skip k=0 vector
    !                 if (kx_idx==0 .and. ky_idx==0 .and. kz_idx==0) cycle

    !                 k_squared = (dble(kx_idx)/dble(ewald%kmax(1)))**2 + &
    !                                 (dble(ky_idx)/dble(ewald%kmax(2)))**2 + &
    !                                 (dble(kz_idx)/dble(ewald%kmax(3)))**2

    !                 ! Skip k=0 vector and vectors outside the unit sphere in normalized space
    !                 if (abs(k_squared) < error .OR. k_squared > one) cycle

    !                 number_of_kpoints = number_of_kpoints + 1
    !                 ewald%recip_amplitude(number_of_kpoints) = ComputeRecipAmplitude(kx_idx, ky_idx, kz_idx)
    !                 u_recipCoulomb = u_recipCoulomb + form_factor * ewald%recip_constants(number_of_kpoints) * &
    !                     real(ewald%recip_amplitude(number_of_kpoints) * CONJG(ewald%recip_amplitude(number_of_kpoints)), KIND=8)
    !             end do
    !         end do
    !     end do

    !     u_recipCoulomb = u_recipCoulomb * EPS0_INV_eVA / KB_eVK * TWOPI / primary%volume

    ! end subroutine ComputeReciprocalEnergy

    !------------------------------------------------------------------------------
    ! ComputeIntraResidueRealCoulombEnergy_singlemol
    ! Computes real-space intramolecular Coulomb energy for a molecule.
    !------------------------------------------------------------------------------
    subroutine ComputeIntraResidueRealCoulombEnergy_singlemol(residue_type, molecule_index, u_intraCoulomb_1)

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

    end subroutine ComputeIntraResidueRealCoulombEnergy_singlemol

    subroutine ComputeReciprocalEnergy_singlemol(residue_type, molecule_index, u_recipCoulomb_new)

        implicit none

        ! Input arguments
        integer, intent(in) :: residue_type
        integer, intent(in) :: molecule_index
        real(real64), intent(out) :: u_recipCoulomb_new

        ! Local variables
        integer :: kx_idx, ky_idx, kz_idx ! Index for reciprocal lattice vector x-y-z-component
        integer :: number_of_kpoints ! Index counter for Fourier components (reciprocal space)
        real(real64) :: form_factor ! Multiplicative factor used in energy calculations
        real(real64) :: k_squared

        u_recipCoulomb_new = 0.0_real64
        number_of_kpoints = 0

        do kx_idx = 0, ewald%kmax(1)

            if (kx_idx == 0) then
                form_factor = 1.0_real64
            else
                form_factor = 2.0_real64
            end if

            do ky_idx = -ewald%kmax(2), ewald%kmax(2)
                do kz_idx = -ewald%kmax(3), ewald%kmax(3)

                    ! Calculate normalized squared magnitude of k-vector
                    k_squared = (dble(kx_idx)/dble(ewald%kmax(1)))**2 + &
                                    (dble(ky_idx)/dble(ewald%kmax(2)))**2 + &
                                    (dble(kz_idx)/dble(ewald%kmax(3)))**2

                    ! Skip k=0 vector and vectors outside the unit sphere in normalized space
                    if (abs(k_squared) < 1.0D-12 .OR. k_squared > 1._real64) cycle

                    number_of_kpoints = number_of_kpoints + 1 ! Increment Fourier component index

                    ! Update Fourier coefficient (recip_amplitude) for current k-vector
                    ewald%recip_amplitude(number_of_kpoints) = ewald%recip_amplitude(number_of_kpoints) + &
                        sum(primary%atom_charges(residue_type, 1:nb%atom_in_residue(residue_type)) * &
                            (ewald%phase_factor_x(residue_type, molecule_index, 1:nb%atom_in_residue(residue_type), kx_idx) * &
                            ewald%phase_factor_y(residue_type, molecule_index, 1:nb%atom_in_residue(residue_type), ky_idx) * &
                            ewald%phase_factor_z(residue_type, molecule_index, 1:nb%atom_in_residue(residue_type), kz_idx) - &
                            ewald%phase_factor_x_old(1:nb%atom_in_residue(residue_type), kx_idx) * &
                            ewald%phase_factor_y_old(1:nb%atom_in_residue(residue_type), ky_idx) * &
                            ewald%phase_factor_z_old(1:nb%atom_in_residue(residue_type), kz_idx)))


                    u_recipCoulomb_new = u_recipCoulomb_new + form_factor * ewald%recip_constants(number_of_kpoints) * &
                                        real(ewald%recip_amplitude(number_of_kpoints) * &
                                        CONJG(ewald%recip_amplitude(number_of_kpoints)), KIND=8)
                end do
            end do
        end do

        u_recipCoulomb_new = u_recipCoulomb_new * EPS0_INV_eVA / KB_eVK * TWOPI / primary%volume

    end subroutine ComputeReciprocalEnergy_singlemol

    subroutine UpdateReciprocalEnergy_creation(residue_type, molecule_index, u_recipCoulomb_new)

        implicit none

        ! Input arguments
        integer, intent(in) :: residue_type
        integer, intent(in) :: molecule_index
        real(real64), intent(out) :: u_recipCoulomb_new

        ! Local variables
        integer :: kx_idx, ky_idx, kz_idx
        integer :: number_of_kpoints
        real(real64) :: form_factor
        real(real64) :: k_squared

        ! Initialize
        u_recipCoulomb_new = 0.0_real64
        number_of_kpoints = 0

        ! Update recip_amplitude by adding (is_creation=.true.) or subtracting (is_creation=.false.) the molecule's contribution
        do kx_idx = 0, ewald%kmax(1)

            if (kx_idx == 0) then
                form_factor = 1.0_real64
            else
                form_factor = 2.0_real64
            end if

            do ky_idx = -ewald%kmax(2), ewald%kmax(2)
                do kz_idx = -ewald%kmax(3), ewald%kmax(3)

                    k_squared = (dble(kx_idx)/dble(ewald%kmax(1)))**2 + &
                                (dble(ky_idx)/dble(ewald%kmax(2)))**2 + &
                                (dble(kz_idx)/dble(ewald%kmax(3)))**2

                    if (abs(k_squared) < 1.0D-12 .OR. k_squared > 1.0_real64) cycle

                    number_of_kpoints = number_of_kpoints + 1

                    ! Update recip_amplitude(number_of_kpoints)
                    ewald%recip_amplitude(number_of_kpoints) = ewald%recip_amplitude(number_of_kpoints) + &
                        sum(primary%atom_charges(residue_type, 1:nb%atom_in_residue(residue_type)) * &
                            ewald%phase_factor_x(residue_type, molecule_index, 1:nb%atom_in_residue(residue_type), kx_idx) * &
                            ewald%phase_factor_y(residue_type, molecule_index, 1:nb%atom_in_residue(residue_type), ky_idx) * &
                            ewald%phase_factor_z(residue_type, molecule_index, 1:nb%atom_in_residue(residue_type), kz_idx))

                    u_recipCoulomb_new = u_recipCoulomb_new + form_factor * ewald%recip_constants(number_of_kpoints) * &
                        real(ewald%recip_amplitude(number_of_kpoints) * CONJG(ewald%recip_amplitude(number_of_kpoints)), KIND=8)
                
                end do
            end do
        end do

        ! Scale energy
        u_recipCoulomb_new = u_recipCoulomb_new * EPS0_INV_eVA / KB_eVK * TWOPI / primary%volume

    end subroutine UpdateReciprocalEnergy_creation

    subroutine UpdateReciprocalEnergy_deletion(residue_type, u_recipCoulomb_new)

        implicit none

        ! Input arguments
        integer, intent(in) :: residue_type
        real(real64), intent(out) :: u_recipCoulomb_new

        ! Local variables
        integer :: kx_idx, ky_idx, kz_idx
        integer :: number_of_kpoints
        real(real64) :: form_factor, k_squared

        ! Initialize
        u_recipCoulomb_new = 0.0_real64
        number_of_kpoints = 0

        ! Update recip_amplitude by adding (is_creation=.true.) or subtracting (is_creation=.false.) the molecule's contribution
        do kx_idx = 0, ewald%kmax(1)

            if (kx_idx == 0) then
                form_factor = 1.0_real64
            else
                form_factor = 2.0_real64
            end if

            do ky_idx = -ewald%kmax(2), ewald%kmax(2)
                do kz_idx = -ewald%kmax(3), ewald%kmax(3)

                    k_squared = (dble(kx_idx)/dble(ewald%kmax(1)))**2 + &
                                (dble(ky_idx)/dble(ewald%kmax(2)))**2 + &
                                (dble(kz_idx)/dble(ewald%kmax(3)))**2

                    if (abs(k_squared) < 1.0D-12 .OR. k_squared > 1.0_real64) cycle

                    number_of_kpoints = number_of_kpoints + 1

                    ! Update recip_amplitude(number_of_kpoints)
                    ewald%recip_amplitude(number_of_kpoints) = ewald%recip_amplitude(number_of_kpoints) - &
                        sum(primary%atom_charges(residue_type, 1:nb%atom_in_residue(residue_type)) * &
                            ewald%phase_factor_x_old(1:nb%atom_in_residue(residue_type), kx_idx) * &
                            ewald%phase_factor_y_old(1:nb%atom_in_residue(residue_type), ky_idx) * &
                            ewald%phase_factor_z_old(1:nb%atom_in_residue(residue_type), kz_idx))

                    u_recipCoulomb_new = u_recipCoulomb_new + form_factor * ewald%recip_constants(number_of_kpoints) * &
                        real(ewald%recip_amplitude(number_of_kpoints) * CONJG(ewald%recip_amplitude(number_of_kpoints)), KIND=8)
                
                end do
            end do
        end do

        ! Scale energy
        u_recipCoulomb_new = u_recipCoulomb_new * EPS0_INV_eVA / KB_eVK * TWOPI / primary%volume

    end subroutine UpdateReciprocalEnergy_deletion

    ! Saves the current Fourier components for a given molecule.
    subroutine SaveFourierTerms_singlemol(residue_type, molecule_index)

        implicit none

        ! Input arguments
        integer, intent(in) :: residue_type
        integer, intent(in) :: molecule_index

        ! Local variables
        integer :: atom_index_1
        integer :: kx_idx, ky_idx, kz_idx
        integer :: number_of_kpoints ! Index counter for Fourier components (reciprocal space)
        real(real64) :: k_squared

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

        ! Save recip_amplitude_old from recip_amplitude
        number_of_kpoints = 0
        do kx_idx = 0, ewald%kmax(1)
            do ky_idx = -ewald%kmax(2), ewald%kmax(2)
                do kz_idx = -ewald%kmax(3), ewald%kmax(3)

                    k_squared = (dble(kx_idx)/dble(ewald%kmax(1)))**2 + &
                                (dble(ky_idx)/dble(ewald%kmax(2)))**2 + &
                                (dble(kz_idx)/dble(ewald%kmax(3)))**2

                    if (abs(k_squared) < 1.0D-12 .OR. k_squared > 1.0_real64) cycle

                    number_of_kpoints = number_of_kpoints + 1
                    ewald%recip_amplitude_old(number_of_kpoints) = ewald%recip_amplitude(number_of_kpoints)
                
                end do
            end do
        end do

    end subroutine SaveFourierTerms_singlemol

    ! Restores the previously saved Fourier components for a given molecule.
    subroutine RestoreFourierState_singlemol(residue_type, molecule_index)

        implicit none

        ! Input arguments
        integer, intent(in) :: residue_type
        integer, intent(in) :: molecule_index

        ! Local variables
        integer :: atom_index_1
        integer :: kx_idx, ky_idx, kz_idx
        integer :: number_of_kpoints ! Index counter for Fourier components
        real(real64) :: k_squared

        ! Restore IKX, IKY, IKZ from backups
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

        ! Restore recip_amplitude from recip_amplitude_old for all valid k-vectors
        number_of_kpoints = 0
        do kx_idx = 0, ewald%kmax(1)
            do ky_idx = -ewald%kmax(2), ewald%kmax(2)
                do kz_idx = -ewald%kmax(3), ewald%kmax(3)

                    k_squared = (dble(kx_idx)/dble(ewald%kmax(1)))**2 + &
                                (dble(ky_idx)/dble(ewald%kmax(2)))**2 + &
                                (dble(kz_idx)/dble(ewald%kmax(3)))**2

                    if (abs(k_squared) < 1.0D-12 .OR. k_squared > 1.0_real64) cycle

                    number_of_kpoints = number_of_kpoints + 1
                    ewald%recip_amplitude(number_of_kpoints) = ewald%recip_amplitude_old(number_of_kpoints)
                
                end do
            end do
        end do

    end subroutine RestoreFourierState_singlemol

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









