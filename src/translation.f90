module molecule_translation

    use monte_carlo_utils
    use ewald_utils
    use simulation_state
    use random_utils
    use energy_utils
    use geometry_utils
    use constants
    use, intrinsic :: iso_fortran_env, only: real64

    implicit none

    private
    public :: Translation

contains

    !---------------------------------------------------------------------------
    ! Subroutine: Translation
    !
    ! Purpose:
    !   Perform a trial translation of a molecule by a random displacement
    !   vector.
    !
    ! Input arguments:
    !   residue_type   - Residue type to be translated
    !   molecule_index - Molecule ID to move
    !
    ! Notes:
    !   - Uses Monte Carlo Metropolis criterion for acceptance/rejection
    !   - Updates energy and Fourier terms as needed
    !   - Applies minimum image convention to new positions
    !---------------------------------------------------------------------------
    subroutine Translation(residue_type, molecule_index)

        implicit none

        ! Input arguments
        integer, intent(in) :: residue_type         ! Residue type to be moved
        integer, intent(in) :: molecule_index       ! Molecule ID
        ! Local variables
        real(real64), dimension(3) :: com_old       ! Previous X, Y, Z coordinate of molecule
        real(real64) :: probability                 ! Acceptance probability of move

        ! Validate translation step size
        call ValidateInputTranslation(input%translation_step)

        ! Exit early if molecule index is zero
        if (molecule_index == 0) return

        ! Increment trial move counter
        counter%trial_translations = counter%trial_translations + 1

        ! Save the current state of the molecule
        call SaveMoleculeState(residue_type, molecule_index, com_old)

        ! Compute old energy of the molecule/system
        call ComputeOldEnergy(residue_type, molecule_index, old)

        ! Propose a random translation move
        call ProposeTranslationMove(residue_type, molecule_index)

        ! Compute new energy after the proposed move
        call ComputeNewEnergy(residue_type, molecule_index, new)

        ! Compute Metropolis acceptance probability
        ! probability = min(1, exp(-(new - old)/kT))
        probability = mc_acceptance_probability(old, new, residue_type, TYPE_TRANSLATION)

        ! Accept or reject
        if (rand_uniform() <= probability) then

            ! Accept move: update system state
            call AcceptTranslationMove(old, new)
        
        else
        
            ! Reject move: restore previous position
            call RejectTranslationMove(residue_type, molecule_index, com_old)
        
        end if

    end subroutine Translation

    subroutine SaveMoleculeState(residue_type, molecule_index, com_old)

        implicit none

        integer, intent(in) :: residue_type         ! Residue type to be moved
        integer, intent(in) :: molecule_index       ! Molecule ID
        real(real64), intent(out) :: com_old(3)     ! Previous X, Y, Z coordinate of molecule

        ! Save current center-of-mass
        com_old(:) = primary%mol_com(:, residue_type, molecule_index)

        ! Save Fourier terms
        call SaveSingleMolFourierTerms(residue_type, molecule_index)

    end subroutine SaveMoleculeState

    !-----------------------------------------------------------------------
    !> ProposeTranslationMove
    !! 
    !! Generate a trial translation move for a selected molecule of a given
    !! residue type. A random displacement vector is drawn uniformly from
    !! [-translation_step/2, translation_step/2) along each coordinate axis,
    !! added to the molecule's previous center of mass, and wrapped into
    !! the simulation box using the minimum image convention.
    !-----------------------------------------------------------------------
    subroutine ProposeTranslationMove(res_type, mol_index)

        implicit none

        ! Input variables
        integer, intent(in) :: res_type             ! Residue type to be moved
        integer, intent(in) :: mol_index            ! Molecule ID
        ! Local variables
        real(real64) :: trial_pos(3)

        ! Generate random move of max size "translation_step/2"
        trial_pos = rand_symmetric(3) * input%translation_step

        ! Apply translation to previous COM position
        primary%mol_com(:, res_type, mol_index) = primary%mol_com(:, res_type, mol_index) + trial_pos(:)

        ! Apply minimum image convension
        call ApplyPBC(primary%mol_com(:, res_type, mol_index), primary)

    end subroutine ProposeTranslationMove

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

    subroutine RejectTranslationMove(residue_type, molecule_index, mol_com_old)

        implicit none

        ! Input arguments
        integer, intent(in) :: residue_type         ! Residue type to be moved
        integer, intent(in) :: molecule_index       ! Index of the molecule to move
        real(real64), intent(in), dimension(3) :: mol_com_old ! For storing old molecule center-of-mass

        ! Restore center-of-mass position
        primary%mol_com(:, residue_type, molecule_index) = mol_com_old(:)

        ! Restore Fourier states (ik_alloc and dk_alloc)
        call RestoreFourierState_singlemol(residue_type, molecule_index)

    end subroutine RejectTranslationMove

    subroutine AcceptTranslationMove(old, new)

        implicit none

        ! Input arguments
        type(energy_state), intent(in) :: old       ! Previous energy states
        type(energy_state), intent(in) :: new       ! New energy states

        ! Update energy totals
        energy%recip_coulomb = energy%recip_coulomb + new%recip_coulomb - old%recip_coulomb
        energy%non_coulomb   = energy%non_coulomb   + new%non_coulomb   - old%non_coulomb
        energy%coulomb       = energy%coulomb       + new%coulomb       - old%coulomb
        energy%total         = energy%total         + new%total         - old%total

        ! Count successful move
        counter%translations = counter%translations + 1

    end subroutine AcceptTranslationMove

    subroutine ValidateInputTranslation(step)

        implicit none

        ! Input arguments
        real(real64), intent(in) :: step         ! Translation step

        ! Local
        character(len=200) :: formatted_msg

        ! Check translation step
        if (step <= zero) then
            write(formatted_msg, '(A,F6.3)') &
                "ERROR: translation_step must be positive. Provided value: ", step
            call AbortRun(trim(formatted_msg))
        end if

    end subroutine ValidateInputTranslation

end module molecule_translation
