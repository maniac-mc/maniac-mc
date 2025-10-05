module prepare_utils

    use simulation_state
    use output_utils
    use ewald_kvectors
    use ewald_phase
    use ewald_energy
    use, intrinsic :: iso_fortran_env, only: real64

    implicit none

contains

    !-----------------------------------------------------------
    ! Subroutine: PrepareSimulationParameters
    ! Purpose: Convert fugacities to dimensionless activities and
    !          initialize Ewald summation parameters for the simulation.
    !-----------------------------------------------------------
    subroutine PrepareSimulationParameters()

        implicit none

        ! Step 1: Convert fugacities to dimensionless activities
        call ConvertFugacity()

        ! Step 2: Initialize Ewald summation parameters
        !         (cutoff, precision, reciprocal space, etc.)
        call SetupEwald()

        ! Step 3: Allocate memory for arrays needed by the Ewald
        !         method (reciprocal vectors, coefficients, etc.)
        call AllocateEwaldArray()

        ! Step 4: Precompute k vectors
        call PrecomputeValidReciprocalVectors()

        ! Step 5: Log the Ewald parameters and settings for
        !         reproducibility and debugging
        call LogEwaldParameters()

    end subroutine PrepareSimulationParameters

    !--------------------------------------------------------------------
    ! Subroutine: ConvertFugacity
    !   Converts the fugacity of each active residue from units of atm
    !   into a dimensionless activity per cubic ångström (Å⁻³).
    !--------------------------------------------------------------------
    subroutine ConvertFugacity()

        use constants

        implicit none

        ! Local variables
        real(real64) :: thermal_energy ! kB * T in Joule
        integer :: id_residue

        thermal_energy = KB_JK * input%temp_K

        do id_residue = 1, nb%type_residue
            if (input%is_active(id_residue) == 1) then

                if (input%fugacity(id_residue) <= 0.0_real64) then
                    call AbortRun("Invalid fugacity for active residue with ID = " // trim(res%names_1d(id_residue)))
                end if

                ! Convert fugacity from atm → Pa → activity in Å⁻³
                input%fugacity(id_residue) = input%fugacity(id_residue) * ATM_TO_PA * A3_TO_M3 / thermal_energy ! atm → Pa

            end if
        end do

    end subroutine ConvertFugacity

    subroutine LogEwaldParameters()

        implicit none

        character(200) :: formatted_msg

        write(formatted_msg, '(A, F10.4)') 'Real-space cutoff (Å): ', input%real_space_cutoff
        call LogMessage(formatted_msg)
        write(formatted_msg, '(A, ES12.5)') 'Ewald accuracy tolerance: ', input%ewald_tolerance
        call LogMessage(formatted_msg)
        write(formatted_msg, '(A, F10.4)') 'Screening factor (dimensionless): ', ewald%screening_factor
        call LogMessage(formatted_msg)
        write(formatted_msg, '(A, F10.4)') 'Ewald damping parameter alpha (1/Å): ', ewald%alpha
        call LogMessage(formatted_msg)
        write(formatted_msg, '(A, F10.4)') 'Fourier-space precision parameter: ', ewald%fourier_precision
        call LogMessage(formatted_msg)
        write(formatted_msg, '(A, I5, A, I5, A, I5)') 'Max Fourier index (kmax(1), kmax(2), kmax(3)): ', &
            ewald%kmax(1), ', ', ewald%kmax(2), ', ', ewald%kmax(3)
        call LogMessage(formatted_msg)
        write(formatted_msg, '(A, I10)') 'Total reciprocal lattice vectors: ', ewald%num_kvectors
        call LogMessage(formatted_msg)

    end subroutine LogEwaldParameters

    !--------------------------------------------------------------------
    ! Subroutine: SetupEwald
    ! Initialize parameters for the Ewald summation method.
    !--------------------------------------------------------------------
    subroutine SetupEwald(verbose)

        implicit none

        ! Input parameters
        logical, optional, intent(in) :: verbose   ! optional flag for logging
        logical :: do_log                          ! Internal flag: whether to log messages in this routine

        ! Default: log unless explicitly disabled
        do_log = .true.; if (present(verbose)) do_log = verbose

        ! Step 1: Adjust cutoff if too large
        call AdjustRealSpaceCutoff(do_log)

        ! Step 2: Clamp tolerance
        call ClampTolerance()

        ! Step 3: Compute Ewald parameters
        call ComputeEwaldParameters()

        ! Step 4: Compute Fourier indices
        call ComputeFourierIndices()

    end subroutine SetupEwald

    !--------------------------------------------------------------------
    ! Subroutine: AdjustRealSpaceCutoff
    ! Purpose   : Ensure that the user-specified real-space cutoff fits
    !             within the periodic simulation box. If it is too large,
    !             reduce it to half the minimum box dimension.
    !--------------------------------------------------------------------
    subroutine AdjustRealSpaceCutoff(do_log)

        logical, intent(in) :: do_log

        character(200) :: msg           ! Buffer for logging
        real(real64) :: safe_cutoff     ! Adjusted real-space cutoff length to fit inside the simulation box safely

        if (input%real_space_cutoff > primary%metrics(1) &
                .OR. input%real_space_cutoff > primary%metrics(2) &
                .OR. input%real_space_cutoff > primary%metrics(3)) then
            if (do_log) then
                write(msg, '(A)') 'WARNING: real_space_cutoff too large for box. Reducing to safe value.'
                call LogMessage(msg)
            end if
            safe_cutoff = min(primary%metrics(1), primary%metrics(2), primary%metrics(3)) / 2.0_real64
            input%real_space_cutoff = safe_cutoff
        end if
    end subroutine AdjustRealSpaceCutoff

    !--------------------------------------------------------------------
    ! Subroutine: ClampTolerance
    ! Purpose   : Limit the Ewald accuracy tolerance to a maximum value of 0.5.
    !--------------------------------------------------------------------
    subroutine ClampTolerance()
        ! Clamp accuracy tolerance to max 0.5
        input%ewald_tolerance = min(abs(input%ewald_tolerance), 0.5_real64)
    end subroutine ClampTolerance

    !--------------------------------------------------------------------
    ! Subroutine: ComputeEwaldParameters
    ! Purpose   : Compute the main Ewald summation parameters:
    !             - screening factor
    !             - damping parameter (alpha)
    !             - Fourier-space precision
    !--------------------------------------------------------------------
    subroutine ComputeEwaldParameters()
        ! Intermediate tolerance factor for screening width
        ewald%screening_factor = sqrt(abs(log(input%ewald_tolerance * input%real_space_cutoff)))

        ! Compute Ewald damping parameter
        ewald%alpha = sqrt(abs(log(input%ewald_tolerance * input%real_space_cutoff * ewald%screening_factor))) / &
                    input%real_space_cutoff

        ! Estimate needed Fourier-space precision
        ewald%fourier_precision = sqrt(-log(input%ewald_tolerance * input%real_space_cutoff * &
                                (2.0_real64 * ewald%screening_factor * ewald%alpha)**2))
    end subroutine ComputeEwaldParameters

    !--------------------------------------------------------------------
    ! Subroutine: ComputeFourierIndices
    ! Purpose   : Compute the maximum Fourier indices (kmax) in the X, Y, Z
    !             directions, and the total number of reciprocal lattice vectors.
    !--------------------------------------------------------------------
    subroutine ComputeFourierIndices()

        implicit none

        integer :: kx_idx, ky_idx, kz_idx   ! Loop indices for each k-component
        integer :: count                    ! Counter for valid k-vectors
        real(real64) :: k_squared           ! Normalized squared magnitude of k-vector

        ! Compute maximum Fourier indices in X, Y, Z directions
        ewald%kmax = nint(0.25_real64 + primary%metrics(1:3) * ewald%alpha * ewald%fourier_precision / PI)

        ! Count the total number of valid reciprocal vectors
        count = 0
        do kx_idx = 0, ewald%kmax(1)
            do ky_idx = -ewald%kmax(2), ewald%kmax(2)
                do kz_idx = -ewald%kmax(3), ewald%kmax(3)
                    if (kx_idx == 0 .and. ky_idx == 0 .and. kz_idx == 0) cycle
                    k_squared = NormalizedKSquared(kx_idx, ky_idx, kz_idx, ewald%kmax)
                    if (CheckValidReciprocalVector(k_squared)) then
                        count = count + 1
                    end if
                end do
            end do
        end do

        ewald%num_kvectors = count

    end subroutine ComputeFourierIndices

    !-------------------------------------------------------------------
    ! Subroutine: ALLOC_TAB
    ! Purpose: Allocate arrays used for Fourier components and related data
    !-------------------------------------------------------------------
    subroutine AllocateEwaldArray()

        implicit none

        ! Allocate real arrays for coefficients (dimension ewald%num_kvectors)
        allocate(ewald%recip_constants(1:ewald%num_kvectors))
        allocate(ewald%recip_amplitude(1:ewald%num_kvectors))
        allocate(ewald%recip_amplitude_old(1:ewald%num_kvectors))
        allocate(ewald%form_factor(ewald%num_kvectors)) ! Note, there is no need for such as large vector

        ! Allocate complex arrays for wave vector components
        allocate(ewald%phase_factor_x(1:nb%type_residue, 0:NB_MAX_MOLECULE, 1:nb%max_atom_in_residue, -ewald%kmax(1):ewald%kmax(1)))
        allocate(ewald%phase_factor_y(1:nb%type_residue, 0:NB_MAX_MOLECULE, 1:nb%max_atom_in_residue, -ewald%kmax(2):ewald%kmax(2)))
        allocate(ewald%phase_factor_z(1:nb%type_residue, 0:NB_MAX_MOLECULE, 1:nb%max_atom_in_residue, -ewald%kmax(3):ewald%kmax(3)))
        allocate(ewald%phase_factor_x_old(1:nb%max_atom_in_residue, -ewald%kmax(1):ewald%kmax(1)))
        allocate(ewald%phase_factor_y_old(1:nb%max_atom_in_residue, -ewald%kmax(2):ewald%kmax(2)))
        allocate(ewald%phase_factor_z_old(1:nb%max_atom_in_residue, -ewald%kmax(3):ewald%kmax(3)))

        ! Allocate kvectors
        allocate(ewald%kvectors(ewald%num_kvectors))

    end subroutine AllocateEwaldArray

end module prepare_utils
