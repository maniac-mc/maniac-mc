module molecule_translation

    use monte_carlo_utils
    use simulation_state
    use random_utils
    use energy_utils
    use ewald_kvectors
    use ewald_phase
    use ewald_energy
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

        ! Exit early if molecule index is zero
        if (molecule_index == 0) return

        ! Increment trial move counter
        counter%trial_translations = counter%trial_translations + 1

        ! Save the current state of the molecule
        call SaveMoleculeState(residue_type, molecule_index, com_old = com_old)

        ! Compute old energy of the molecule/system
        call ComputeOldEnergy(residue_type, molecule_index, old)

        ! Propose a random translation move
        call RandomTranslation(residue_type, molecule_index)

        ! Compute new energy after the proposed move
        call ComputeNewEnergy(residue_type, molecule_index, new)

        ! Compute Metropolis acceptance probability
        ! probability = min(1, exp(-(new - old)/kT))
        probability = mc_acceptance_probability(old, new, residue_type, TYPE_TRANSLATION)

        ! Accept or reject
        if (rand_uniform() <= probability) then

            ! Accept move: update system state
            call AcceptMove(old, new, counter%translations)
        
        else
        
            ! Reject move: restore previous position
            call RejectMoleculeMove(residue_type, molecule_index, com_old = com_old)

        end if

    end subroutine Translation

    !-----------------------------------------------------------------------
    !> RandomTranslation
    !! 
    !! Generate a trial translation move for a selected molecule of a given
    !! residue type. A random displacement vector is drawn uniformly from
    !! [-translation_step/2, translation_step/2) along each coordinate axis,
    !! added to the molecule's previous center of mass, and wrapped into
    !! the simulation box using the minimum image convention.
    !-----------------------------------------------------------------------
    subroutine RandomTranslation(res_type, mol_index)

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

    end subroutine RandomTranslation

end module molecule_translation
