module ewald_utils

    use simulation_state
    use output_utils 
    use constants
    use geometry_utils
    use, intrinsic :: iso_fortran_env, only: real64

    implicit none

contains

    subroutine InitializeReciprocalWeights()

        implicit none

        real(real64) :: k_squared
        integer :: number_of_kpoints
        real(real64) :: kx_vector_x, kx_vector_y, kx_vector_z
        real(real64) :: ky_vector_x, ky_vector_y, ky_vector_z
        real(real64) :: kz_vector_x, kz_vector_y, kz_vector_z
        integer :: kx_idx, ky_idx, kz_idx ! Index for reciprocal lattice vector x-y-z-component

        ! Computes a weighting factor for each reciprocal lattice vector
        ! Only needs to be done once at the start of the simulations
        number_of_kpoints = 0
        do kx_idx = 0, ewald%kmax(1)
            kx_vector_x = TWOPI * dble(kx_idx) * primary%reciprocal(1,1)
            kx_vector_y = TWOPI * dble(kx_idx) * primary%reciprocal(2,1)
            kx_vector_z = TWOPI * dble(kx_idx) * primary%reciprocal(3,1)
            do ky_idx = -ewald%kmax(2), ewald%kmax(2) 
                ky_vector_x = kx_vector_x + TWOPI * dble(ky_idx) * primary%reciprocal(1,2)
                ky_vector_y = kx_vector_y + TWOPI * dble(ky_idx) * primary%reciprocal(2,2)
                ky_vector_z = kx_vector_z + TWOPI * dble(ky_idx) * primary%reciprocal(3,2)
                do kz_idx = -ewald%kmax(3), ewald%kmax(3) 
                    kz_vector_x = ky_vector_x + TWOPI * dble(kz_idx) * primary%reciprocal(1,3)
                    kz_vector_y = ky_vector_y + TWOPI * dble(kz_idx) * primary%reciprocal(2,3)
                    kz_vector_z = ky_vector_z + TWOPI * dble(kz_idx) * primary%reciprocal(3,3)

                    ! Square root of k
                    k_squared = (dble(kx_idx)/dble(ewald%kmax(1)))**2 + &
                                (dble(ky_idx)/dble(ewald%kmax(2)))**2 + &
                                (dble(kz_idx)/dble(ewald%kmax(3)))**2

                    ! Skip k=0 vector and vectors outside the unit sphere in normalized space
                    if (abs(k_squared) < 1.0D-12 .OR. k_squared > 1._real64) cycle

                    number_of_kpoints = number_of_kpoints + 1
                    ewald%recip_constants(number_of_kpoints) = &
                        EXP(-(kz_vector_x**2 + kz_vector_y**2 + kz_vector_z**2) / &
                        (4.0_real64 * ewald%alpha**2)) / (kz_vector_x**2 + kz_vector_y**2 + kz_vector_z**2)
                end do
            end do
        end do

        if (number_of_kpoints > ewald%num_recip_vectors) stop

        return

    end subroutine InitializeReciprocalWeights

    ! Precompute complex exponential factors to be used repeatedly in
    ! reciprocal-space calculations.
    subroutine ComputeFourierTerms_singlemol(residue_type, molecule_index)

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

    end subroutine ComputeFourierTerms_singlemol

    subroutine ComputeReciprocalEnergy(u_recipCoulomb)

        implicit none

        ! Input arguments
        real(real64), intent(out) :: u_recipCoulomb

        ! Internal variables
        integer :: number_of_kpoints
        integer :: kx_idx, ky_idx, kz_idx ! Index for reciprocal lattice vector x-y-z-component
        integer :: atom_index_1, molecule_index, residue_type
        real(real64) :: form_factor, k_squared

        ! Compute reciprocal-space energy sum over allowed k-vectors excluding zero and outside cutoff
        number_of_kpoints = 0
        u_recipCoulomb = 0
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

                    ! Skip k=0 vector and vectors outside the unit sphere in normalized space
                    if (abs(k_squared) < 1.0D-12 .OR. k_squared > 1._real64) cycle

                    number_of_kpoints = number_of_kpoints + 1
                    ewald%recip_amplitude(number_of_kpoints) = (0.0_real64, 0.0_real64)

                    ! Loop over all molecule types 1
                    do residue_type = 1, nb%type_residue
                        ! Loop over all molecule index 1
                        do molecule_index = 1, primary%num_residues(residue_type)
                            ! Loop over sites in molecule residue_type
                            do atom_index_1 = 1, nb%atom_in_residue(residue_type) 
                                ewald%recip_amplitude(number_of_kpoints) = ewald%recip_amplitude(number_of_kpoints) &
                                    + primary%atom_charges(residue_type, atom_index_1) * &
                                    ewald%phase_factor_x(residue_type, molecule_index, atom_index_1, kx_idx) * &
                                    ewald%phase_factor_y(residue_type, molecule_index, atom_index_1, ky_idx) * &
                                    ewald%phase_factor_z(residue_type, molecule_index, atom_index_1, kz_idx)
                            end do
                        end do
                    end do

                    u_recipCoulomb = u_recipCoulomb + form_factor * ewald%recip_constants(number_of_kpoints) * &
                        real(ewald%recip_amplitude(number_of_kpoints) * CONJG(ewald%recip_amplitude(number_of_kpoints)), KIND=8)
                end do
            end do
        end do

        u_recipCoulomb = u_recipCoulomb * EPS0_INV_eVA / KB_eVK * TWOPI / primary%volume

    end subroutine ComputeReciprocalEnergy

    !------------------------------------------------------------------------------
    ! ComputeEwaldSelfInteraction_singlemol
    ! Computes the Ewald self-energy correction for a single molecule.
    !------------------------------------------------------------------------------
    subroutine ComputeEwaldSelfInteraction_singlemol(residue_type, self_energy_1)
        
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

    end subroutine ComputeEwaldSelfInteraction_singlemol

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









