module monte_carlo_utils

    use simulation_state
    use random_utils
    use constants
    use output_utils
    use, intrinsic :: iso_fortran_env, only: real64

    implicit none

contains

    subroutine ApplyRandomRotation(residue_type, molecule_index, full_rotation)

        implicit none

        ! Input arguments
        integer, intent(in) :: residue_type, molecule_index
        logical, intent(in), optional :: full_rotation

        ! Local variables
        integer :: rotation_axis
        logical :: use_full_rotation
        real(real64) :: rotation_matrix(3, 3)
        real(real64) :: theta
        real(real64) :: cos_theta, sin_theta
        character(len=32) :: axis_str

        ! Handle optional argument: default = .false.
        use_full_rotation = .false.
        if (PRESENT(full_rotation)) use_full_rotation = full_rotation

        ! Exit if single-atom residue (nothing to rotate)
        if (nb%atom_in_residue(residue_type) == 1) return

        if (.not. use_full_rotation) then
            ! For small-step mode, make sure rotation_step_angle is reasonable
            if (input%rotation_step_angle <= 0.0_real64 .OR. input%rotation_step_angle > 2.0_real64 * PI) then
                stop 'Invalid rotation_step_angle in ApplyRandomRotation'
            end if
            ! Use small rotation
            theta = (rand_uniform() - 0.5_real64) * input%rotation_step_angle
        else ! Use large roration
            theta = rand_uniform() * 2.0_real64 * PI
        end if

        ! Choose random axis (1=X, 2=Y, 3=Z)
        rotation_axis = INT(rand_uniform() * 3.0_real64) + 1 ! Random integer in [1,3]
        
        ! Compute trigonometric values
        cos_theta = cos(theta)
        sin_theta = sin(theta)

        ! Initialize rotation matrix (identity matrix)
        rotation_matrix = 0.0_real64
        rotation_matrix(1, 1) = 1.0_real64
        rotation_matrix(2, 2) = 1.0_real64
        rotation_matrix(3, 3) = 1.0_real64

        ! Set rotation matrix based on axis
        select case (rotation_axis)
        case (1)  ! X-axis
            rotation_matrix(2, 2) = cos_theta
            rotation_matrix(2, 3) = -sin_theta
            rotation_matrix(3, 2) = sin_theta
            rotation_matrix(3, 3) = cos_theta
        case (2)  ! Y-axis
            rotation_matrix(1, 1) = cos_theta
            rotation_matrix(1, 3) = sin_theta
            rotation_matrix(3, 1) = -sin_theta
            rotation_matrix(3, 3) = cos_theta
        case (3)  ! Z-axis
            rotation_matrix(1, 1) = cos_theta
            rotation_matrix(1, 2) = -sin_theta
            rotation_matrix(2, 1) = sin_theta
            rotation_matrix(2, 2) = cos_theta
        case default
            write(axis_str, '(I0)') rotation_axis   ! convert to string if integer
            call AbortRun("Invalid rotation_axis: " // trim(axis_str))
        end select

        ! Apply rotation to all atoms
        primary%site_offset(:, residue_type, molecule_index, 1:nb%atom_in_residue(residue_type)) = &
            matmul(rotation_matrix, &
                   primary%site_offset(:, residue_type, molecule_index, 1:nb%atom_in_residue(residue_type)))

    end subroutine ApplyRandomRotation

    !> Adjusts the Monte Carlo translation and rotation step sizes dynamically
    !! This subroutine recalibrates the translational and rotational move steps
    !! based on the observed acceptance ratios of recent MC trials. If the
    !! acceptance is too high, the step size is increased (up to a defined maximum).
    !! If the acceptance is too low, the step size is decreased (down to a defined minimum).
    subroutine AdjustMoveStepSizes()

        implicit none

        real(real64) :: acc_trans, acc_rot

        if (input%recalibrate_moves) then

            ! Adjust translation step
            if (counter%trial_translations > MIN_TRIALS_FOR_RECALIBRATION) then
                acc_trans = real(counter%translations) / real(counter%trial_translations)
                if (acc_trans - TARGET_ACCEPTANCE > TOL_ACCEPTANCE) then
                    input%translation_step = min(input%translation_step * 1.05d0, MAX_TRANSLATION_STEP)
                else if (acc_trans - TARGET_ACCEPTANCE < TOL_ACCEPTANCE) then
                    input%translation_step = max(input%translation_step * 0.95d0, MIN_TRANSLATION_STEP)
                end if

            end if

            ! Adjust rotational step
            if (counter%trial_rotations > MIN_TRIALS_FOR_RECALIBRATION) then
                acc_rot = real(counter%rotations) / real(counter%trial_rotations)
                if (acc_rot - TARGET_ACCEPTANCE > TOL_ACCEPTANCE) then
                    input%rotation_step_angle = min(input%rotation_step_angle * 1.05d0, MAX_ROTATION_ANGLE)
                else if (acc_rot - TARGET_ACCEPTANCE < TOL_ACCEPTANCE) then
                    input%rotation_step_angle = min(input%rotation_step_angle * 1.95d0, MIN_ROTATION_ANGLE)
                end if

            end if
        end if

    end subroutine AdjustMoveStepSizes

    !----------------------------------------------------------------------
    ! PickRandomResidueType: randomly selects an active residue type from the
    ! available types in [1, active_residue_count].
    !----------------------------------------------------------------------
    function PickRandomResidueType(is_active) result(residue_type)
        integer, dimension(:), intent(in) :: is_active
        integer :: residue_type
        integer :: i, n_active
        integer, allocatable :: active_indices(:)

        ! Count active residues
        n_active = count(is_active == 1)

        if (n_active == 0) then
            residue_type = 0  ! or some error code
            return
        end if

        ! Collect indices of active residues
        allocate(active_indices(n_active))
        n_active = 0
        do i = 1, size(is_active)
            if (is_active(i) == 1) then
                n_active = n_active + 1
                active_indices(n_active) = i
            end if
        end do

        ! Pick a random index among active ones
        residue_type = active_indices(INT(rand_uniform() * n_active) + 1)

        deallocate(active_indices)
    end function PickRandomResidueType

    !----------------------------------------------------------------------
    ! PickRandomMoleculeIndex: randomly selects a molecule index within a given
    ! residue type. Returns 0 if there are no molecules of that type.
    !----------------------------------------------------------------------
    function PickRandomMoleculeIndex(residue_count_for_type) result(molecule_index)

        integer, intent(in) :: residue_count_for_type
        integer :: molecule_index

        if (residue_count_for_type == 0) then
            molecule_index = 0
        else
            molecule_index = INT(rand_uniform() * residue_count_for_type) + 1
            if (molecule_index > residue_count_for_type) molecule_index = residue_count_for_type
        end if

    end function PickRandomMoleculeIndex

    function mc_acceptance_probability(old, new, residue_type, move_type) result(probability)

        implicit none

        type(energy_state), intent(in) :: old   ! Previous energy states
        type(energy_state), intent(in) :: new   ! New energy states
        integer, intent(in) :: move_type        ! MC move type: TYPE_CREATION, TYPE_DELETION, TYPE_TRANSLATION, TYPE_ROTATION
        integer, intent(in) :: residue_type     ! Index of the residue type
        real(real64) :: probability             ! Calculated acceptance probability (0 <= P <= 1)
        real(real64) :: one = 1.0_real64        ! Constant 1.0, used for min() comparison
        real(real64) :: N                       ! Number of residues of this type (local copy)
        real(real64) :: V                       ! Simulation box volume (local copy)
        real(real64) :: phi                     ! Fugacity of the residue type (local copy)
        real(real64) :: T                       ! Temperature in units of kB*T (optional local copy for clarity)
        real(real64) :: delta_e                 ! Energy difference ΔE between trial and current state

        N   = real(primary%num_residues(residue_type))
        V   = primary%volume
        phi = input%fugacity(residue_type)
        T = input%temp_K

        delta_e = new%total - old%total

        ! Compute factor based on move type
        select case (move_type)
            case (TYPE_CREATION)

                probability = min(one, (phi * V / N) * exp(-delta_e / T))

            case (TYPE_DELETION)

                probability = min(one, ((N + one) / (phi * V)) * exp(-delta_e / T))

            case (TYPE_TRANSLATION, TYPE_ROTATION)

                probability = min(one, exp(-delta_e / T))

            case default
                
                call AbortRun("Unknown move_type in mc_acceptance_probability!", 1)

        end select

    end function mc_acceptance_probability

    function mc_acceptance_probability_swap(old, new, type_old, type_new) result(probability)

        use, intrinsic :: iso_fortran_env, only: real64
        implicit none

        ! Arguments
        type(energy_state), intent(in) :: old   ! Energy of system with old molecule
        type(energy_state), intent(in) :: new   ! Energy of system with new molecule
        integer, intent(in) :: type_old         ! Residue type being removed
        integer, intent(in) :: type_new         ! Residue type being inserted
        integer :: N_old, N_new                 ! Number of molecule per types

        ! Return value
        real(real64) :: probability             ! Acceptance probability (0 <= P <= 1)

        ! Locals
        real(real64) :: delta_e                 ! Energy difference
        real(real64) :: phi_old, phi_new        ! Fugacities of species
        real(real64) :: T                       ! Temperature in units of kB*T
        real(real64), parameter :: one = 1.0_real64
        real(real64) :: combinatorial           ! Prefactor N_A / (N_B + 1)

        N_new = primary%num_residues(type_new)
        N_old = primary%num_residues(type_old)

        ! Compute energy difference
        delta_e = new%total - old%total

        ! Fugacities
        phi_old = input%fugacity(type_old)
        phi_new = input%fugacity(type_new)

        ! Temperature
        T = input%temp_K

        ! combinatorial factor N_A / (N_B + 1)
        combinatorial = real(N_old, real64) / real(N_new + one, real64)

        ! Swap acceptance probability
        probability = min(one, (phi_new / phi_old) * combinatorial * exp(-delta_e / T))

    end function mc_acceptance_probability_swap


end module monte_carlo_utils
