module monte_carlo_utils

    use simulation_state
    use random_utils
    use constants
    use ewald_utils
    use output_utils
    use energy_utils
    use helper_utils
    use, intrinsic :: iso_fortran_env, only: real64

    implicit none

contains

    !========================================================
    ! Subroutine: ApplyRandomRotation
    !
    ! Rotates all atoms of a residue around a randomly chosen
    ! axis (X, Y, or Z) by a random angle. Angle can be a 
    ! small perturbation or a full 0-2π rotation.
    !
    ! Inputs:
    !   res_type - integer, index of the residue type
    !   mol_index - integer, index of the molecule
    !   full_rotation - optional logical, if .true. use full rotation
    !========================================================
    subroutine ApplyRandomRotation(res_type, mol_index, full_rotation)

        implicit none

        ! Input arguments
        integer, intent(in) :: res_type       ! Index of the residue type
        integer, intent(in) :: mol_index     ! Index of the molecule
        logical, intent(in), optional :: full_rotation ! Flag for full or small-step rotation
        ! Local variables
        integer :: rotation_axis                  ! Chosen rotation axis (1=X, 2=Y, 3=Z)
        integer :: n_atoms                        ! Number of atoms in residue
        logical :: use_full_rotation              ! Actual value of full_rotation (defaults to .false.)
        real(real64) :: rotation_matrix(3,3)      ! 3x3 rotation matrix for the rotation
        real(real64) :: theta                     ! Rotation angle in radians

        ! Handle optional argument: default = .false.
        use_full_rotation = .false.
        if (present(full_rotation)) use_full_rotation = full_rotation

        ! Exit if single-atom residue (nothing to rotate)
        n_atoms = nb%atom_in_residue(res_type)
        if (n_atoms == 1) return

        ! Choose rotation angle
        theta = ChooseRotationAngle(use_full_rotation)

        ! Choose random axis (1=X, 2=Y, 3=Z)
        rotation_axis = int(rand_uniform() * three) + 1 ! Random integer in [1,3]

        ! Set rotation matrix based on axis
        rotation_matrix = RotationMatrix(rotation_axis, theta)

        ! Apply rotation to all atoms in the residue
        primary%site_offset(:, res_type, mol_index, 1:n_atoms) = &
            matmul(rotation_matrix, primary%site_offset(:, res_type, mol_index, 1:n_atoms))

    end subroutine ApplyRandomRotation

    ! Returns a random rotation angle (radians); small-step if use_full_rotation=.false., full [0,2π] if .true.
    function ChooseRotationAngle(use_full_rotation) result(theta)

        implicit none

        ! Input variables
        logical, intent(in) :: use_full_rotation
        ! Output variables
        real(real64) :: theta

        if (.not. use_full_rotation) then

            ! For small-step mode, make sure rotation_step_angle is reasonable
            if (input%rotation_step_angle <= zero .or. input%rotation_step_angle > TWOPI) then
                call AbortRun('Invalid rotation_step_angle in ChooseRotationAngle')
            end if

            ! Use small rotation
            theta = (rand_uniform() - half) * input%rotation_step_angle
        else
            ! Use large rotation
            theta = rand_uniform() * TWOPI
        end if

    end function ChooseRotationAngle

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

    !---------------------------------------------------------------------------
    ! Purpose:
    !   Compute the updated energy of a single molecule after a trial move 
    !   (translation or rotation) for use in the Monte Carlo acceptance test.
    !---------------------------------------------------------------------------
    subroutine ComputeNewEnergy(residue_type, molecule_index, new)

        implicit none

        integer, intent(in) :: residue_type         ! Residue type to be moved
        integer, intent(in) :: molecule_index       ! Index of the molecule to move
        type(energy_state), intent(out) :: new      ! New energy states

        ! Recompute Fourier terms for the moved molecule
        call SingleMolFourierTerms(residue_type, molecule_index)

        ! Compute new energies
        call ComputeReciprocalEnergy_singlemol(residue_type, molecule_index, new%recip_coulomb)
        call ComputePairInteractionEnergy_singlemol(primary, residue_type, molecule_index, new%non_coulomb, new%coulomb)

        ! Total energy of the molecule
        new%total = new%non_coulomb + new%coulomb + new%recip_coulomb

    end subroutine ComputeNewEnergy

    !---------------------------------------------------------------------------
    ! Purpose:
    !   Compute the current (pre-move) energy of a single molecule for use
    !   in the Monte Carlo Metropolis criterion.
    !
    ! Inputs:
    !   residue_type   - Integer: Type of residue/molecule
    !   molecule_index - Integer: Index of the molecule in the system
    !
    ! Outputs:
    !   old            - Type(energy_state): Energy components before trial move
    !
    ! Notes:
    !   - Computes reciprocal-space Coulomb energy using Ewald summation.
    !   - Computes pairwise interaction energies (van der Waals + real-space Coulomb).
    !   - Translation or rotation does not affect self-energy or intra-molecular Coulomb terms.
    !---------------------------------------------------------------------------
    subroutine ComputeOldEnergy(residue_type, molecule_index, old)

        implicit none

        integer, intent(in) :: residue_type         ! Residue type to be moved
        integer, intent(in) :: molecule_index       ! Index of the molecule to move
        type(energy_state), intent(out) :: old      ! Old energy states

        ! Compute current energy (cf e_old)
        call ComputeReciprocalEnergy_singlemol(residue_type, molecule_index, old%recip_coulomb)
        call ComputePairInteractionEnergy_singlemol(primary, residue_type, molecule_index, old%non_coulomb, old%coulomb)

        ! Note: Translation does not affect ewald_self or intra_coulomb
        old%total = old%non_coulomb + old%coulomb + old%recip_coulomb

    end subroutine ComputeOldEnergy

    !---------------------------------------------------------------------------
    ! Purpose:
    !   Update the global system energy and Monte Carlo counters after 
    !   accepting a trial move (translation or rotation) of a molecule.
    !
    ! Inputs:
    !   old         - Type(energy_state): Energy of the molecule before the move
    !   new         - Type(energy_state): Energy of the molecule after the move
    !
    ! Input/Output:
    !   counter_var - Integer: Monte Carlo counter for successful moves
    !                 (e.g., translations or rotations), incremented if move accepted
    !---------------------------------------------------------------------------
    subroutine AcceptMove(old, new, counter_var)
        
        ! Input arguments
        type(energy_state), intent(in) :: old, new
        integer, intent(inout) :: counter_var
        
        energy%recip_coulomb = energy%recip_coulomb + new%recip_coulomb - old%recip_coulomb
        energy%non_coulomb = energy%non_coulomb + new%non_coulomb - old%non_coulomb
        energy%coulomb = energy%coulomb + new%coulomb - old%coulomb
        energy%total = energy%total + new%total - old%total
        counter_var = counter_var + 1
    
    end subroutine

    !---------------------------------------------------------------------------
    ! Subroutine: SaveMoleculeState
    !
    ! Purpose:
    !   Save the current state of a molecule for a Monte Carlo move.
    !   - For translation: saves center-of-mass (COM)
    !   - For rotation: saves site offsets
    !   - Also saves Fourier terms for later restoration if the move is rejected
    !
    ! Inputs:
    !   residue_type   - Integer: Residue type of the molecule
    !   molecule_index - Integer: Index of the molecule
    !
    ! Outputs:
    !   com_old         - Real(3) array: center-of-mass (for translation)
    !   site_offset_old - Real(3, natoms) array: site offsets (for rotation)
    ! Notes:
    !   - Only one of com_old or site_offset_old is allocated by the caller
    !   - Fourier terms are always saved
    !---------------------------------------------------------------------------
    subroutine SaveMoleculeState(res_type, mol_index, com_old, site_offset_old)

        implicit none

        ! Input arguments
        integer, intent(in) :: res_type         ! Residue type of the molecule
        integer, intent(in) :: mol_index        ! Index of the molecule
        real(real64), intent(out), optional :: com_old(3) ! Center-of-mass before move (for translation)
        real(real64), intent(out), optional :: site_offset_old(:, :) ! Site offsets before move (for rotation)
        ! Local variable
        integer :: natoms

        ! Save Fourier terms
        call SaveSingleMolFourierTerms(res_type, mol_index)

        ! Save center-of-mass if requested (translation)
        if (present(com_old)) then
            com_old(:) = primary%mol_com(:, res_type, mol_index)
        end if

        ! Save site offsets if requested (rotation)
        if (present(site_offset_old)) then
            natoms = nb%atom_in_residue(res_type)
            site_offset_old(:, 1:natoms) = primary%site_offset(:, res_type, mol_index, 1:natoms)
        end if

    end subroutine SaveMoleculeState

    !---------------------------------------------------------------------------
    ! Subroutine: RejectMoleculeMove
    !
    ! Purpose:
    !   Restore a molecule's previous state (COM or site offsets) and Fourier terms
    !   if a Monte Carlo move is rejected.
    !
    ! Inputs:
    !   res_type - Residue type of the molecule
    !   mol_index - Index of the molecule
    !   com_old - Previous COM (for translation, optional)
    !   site_offset_old - Previous site offsets (for rotation, optional)
    !---------------------------------------------------------------------------
    subroutine RejectMoleculeMove(res_type, mol_index, com_old, site_offset_old)
        implicit none

        integer, intent(in) :: res_type       ! Residue type of the molecule
        integer, intent(in) :: mol_index      ! Index of the molecule
        real(real64), intent(in), optional :: com_old(3) ! Previous COM (translation)
        real(real64), intent(in), optional :: site_offset_old(:, :) ! Previous site offsets (rotation)
        ! Local variable
        integer :: natoms

        ! Restore COM if present (translation)
        if (present(com_old)) then
            primary%mol_com(:, res_type, mol_index) = com_old(:)
        end if

        ! Restore site offsets if present (rotation)
        if (present(site_offset_old)) then
            natoms = nb%atom_in_residue(res_type)
            primary%site_offset(:, res_type, mol_index, 1:natoms) = &
                site_offset_old(:, 1:natoms)
        end if

        ! Restore Fourier states
        call RestoreFourierState_singlemol(res_type, mol_index)
    end subroutine RejectMoleculeMove

end module monte_carlo_utils
