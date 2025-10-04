module molecule_rotation

    use monte_carlo_utils
    use simulation_state
    use ewald_utils
    use random_utils
    use energy_utils
    use constants
    use, intrinsic :: iso_fortran_env, only: real64

    implicit none

contains

    !---------------------------------------------------------------------------
    ! Subroutine: Rotation
    !
    ! Purpose:
    !   Performs a trial rotation of a molecule around a random axis by a random
    !   angle. The move is accepted or rejected based on the Metropolis 
    !   criterion using energy differences.
    !
    ! Inputs:
    !   residue_type   - Integer: Residue type of the molecule
    !   molecule_index - Integer: Index of the molecule
    !
    ! Outputs:
    !   Updates primary%site_offset if move accepted
    !   Updates energy totals if move accepted
    !   Updates Monte Carlo counters
    !---------------------------------------------------------------------------
    subroutine Rotation(residue_type, molecule_index)

        implicit none

        ! Input arguments
        integer, intent(in) :: residue_type        ! Residue type to be moved
        integer, intent(in) :: molecule_index     ! Molecule ID
        ! Local variables
        real(real64), dimension(:, :), allocatable :: site_offset_old    ! Local site X Y Z displacements from molecule center
        real(real64) :: probability                              ! Acceptance probability of rotation move
        
        ! Exit early if molecule cannot rotate
        if ((nb%atom_in_residue(residue_type) == 1) .or. (molecule_index == 0)) return

        ! Allocate site_offset_old
        allocate(site_offset_old(3, nb%max_atom_in_residue))

        ! Count trial move (success + fail)
        counter%trial_rotations = counter%trial_rotations + 1

        call SaveMoleculeState(residue_type, molecule_index, site_offset_old = site_offset_old)

        ! Compute previous energy
        call ComputeOldEnergy(residue_type, molecule_index, old)

        ! Rotate the molecule randomly
        call ApplyRandomRotation(residue_type, molecule_index)

        ! Recompute energies
        call ComputeNewEnergy(residue_type, molecule_index, new)

        ! Compute acceptance probability for the move
        probability = mc_acceptance_probability(old, new, residue_type, TYPE_ROTATION)

        ! Accept or reject
        if (rand_uniform() <= probability) then ! Accept move

            call AcceptMove(old, new, counter%rotations)
        
        else ! Reject move
        
            call RejectMoleculeMove(residue_type, molecule_index, site_offset_old = site_offset_old)
        
        end if
        
    end subroutine Rotation

end module molecule_rotation
