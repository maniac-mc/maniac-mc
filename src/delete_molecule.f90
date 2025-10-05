module molecule_deletion

    use monte_carlo_utils
    use simulation_state
    use geometry_utils
    use random_utils
    use energy_utils
    use ewald_kvectors
    use ewald_phase
    use ewald_energy
    use constants
    use, intrinsic :: iso_fortran_env, only: real64

    implicit none

    private                    ! Hide everything by default
    public :: DeleteMolecule   ! Only expose the main entry point

contains

    !---------------------------------------------------------------------------
    ! Subroutine: DeleteMolecule
    !
    ! Purpose:
    !   Attempts to remove (delete) an existing molecule of a given residue type
    !   from the simulation box, following the Monte Carlo deletion move scheme.
    !   The move is accepted or rejected based on the Metropolis criterion
    !   using the energy difference and the fugacity of the species.
    !
    !---------------------------------------------------------------------------
    subroutine DeleteMolecule(residue_type, molecule_index)

        implicit none

        ! Input arguments
        integer, intent(in) :: residue_type     ! Residue type to be moved
        integer, intent(in) :: molecule_index   ! Molecule ID

        ! Local variables
        real(real64), dimension(:, :), allocatable :: site_offset_old ! For storing old site offsets
        real(real64), dimension(3) :: mol_com_old ! For storing old molecule center-of-mass
        real(real64) :: probability ! Acceptance probability of creation move
        integer :: last_molecule_index ! Index of the last molecule in the primary box
        logical :: is_deletion                  ! Flag indicating creation

        ! Return immediately if no molecules of type residue_type exist
        if (primary%num_residues(residue_type)==0) return

        ! Count trial move (success + fail)
        counter%trial_deletions = counter%trial_deletions + 1

        allocate(site_offset_old(3, nb%max_atom_in_residue))

        ! Energy of the previous configuration
        is_deletion = .true.
        call ComputeOldEnergy(residue_type, molecule_index, old, is_deletion = is_deletion)

        call SaveMoleculeState(residue_type, molecule_index, com_old = mol_com_old, offset_old = site_offset_old)

        ! Record the index of the last molecule
        last_molecule_index = primary%num_residues(residue_type)

        ! Delete molecule
        call RemoveMolecule(residue_type, molecule_index, last_molecule_index)

        ! Update molecule and atom counts
        primary%num_residues(residue_type) = primary%num_residues(residue_type) - 1
        primary%num_atoms = primary%num_atoms - nb%atom_in_residue(residue_type)

        ! Calculate new energy
        call ComputeNewEnergy(residue_type, molecule_index, new, is_deletion = is_deletion)

        ! Compute acceptance probability for the move
        probability = mc_acceptance_probability(old, new, residue_type, TYPE_DELETION)

        ! Accept or reject
        if (rand_uniform() <= probability) then ! Accept move
            call AcceptDeletionMove(residue_type, last_molecule_index, old, new)
        else ! Reject move
            call RejectDeletionMove(residue_type, molecule_index, mol_com_old, site_offset_old)
        end if

    end subroutine DeleteMolecule

    ! subroutine SaveMoleculeState(residue_type, molecule_index, mol_com_old, site_offset_old)

    !     implicit none

    !     integer, intent(in) :: residue_type      ! Residue type to remove
    !     integer, intent(in) :: molecule_index    ! Molecule index to remove
    !     real(real64), intent(out) :: mol_com_old(3) ! For storing old molecule center-of-mass
    !     real(real64), intent(out), dimension(:, :) :: site_offset_old

    !     ! Store positions and site offsets
    !     mol_com_old(:) = primary%mol_com(:, residue_type, molecule_index)
    !     site_offset_old(:, 1:nb%atom_in_residue(residue_type)) = &
    !         primary%site_offset(:, residue_type, molecule_index, 1:nb%atom_in_residue(residue_type))

    !     ! Save Fourier terms
    !     call SaveFourierTerms_singlemol(residue_type, molecule_index)

    ! end subroutine SaveMoleculeState

    subroutine RemoveMolecule(residue_type, molecule_index, last_molecule_index)

        implicit none

        integer, intent(in) :: residue_type      ! Residue type to remove
        integer, intent(in) :: molecule_index    ! Molecule index to remove
        integer, intent(in):: last_molecule_index ! Index of the last molecule in the primary box

        ! Replace with the last molecule
        primary%mol_com(:, residue_type, molecule_index) = &
            primary%mol_com(:, residue_type, last_molecule_index)
        primary%site_offset(:, residue_type, molecule_index, 1:nb%atom_in_residue(residue_type)) = &
            primary%site_offset(:, residue_type, last_molecule_index, 1:nb%atom_in_residue(residue_type))

        ! Replace Fourier terms
        call ReplaceFourierTermsSingleMol(residue_type, molecule_index, last_molecule_index)

    end subroutine RemoveMolecule

    subroutine AcceptDeletionMove(residue_type, last_molecule_index, old, new)

        implicit none

        integer, intent(in) :: residue_type         ! Residue type to be moved
        integer, intent(in) :: last_molecule_index  ! Index of the last molecule in the reservoir
        type(energy_state), intent(in) :: old       ! Previous energy states
        type(energy_state), intent(in) :: new       ! New energy states
        real(real64) :: trial_pos(3)                ! Random numbers for initial molecule position in the box

        ! Update total energies
        energy%recip_coulomb    = new%recip_coulomb
        energy%non_coulomb      = energy%non_coulomb    + new%non_coulomb   - old%non_coulomb
        energy%coulomb          = energy%coulomb        + new%coulomb       - old%coulomb
        energy%ewald_self       = energy%ewald_self     + new%ewald_self    - old%ewald_self
        energy%intra_coulomb    = energy%intra_coulomb  + new%intra_coulomb - old%intra_coulomb
        energy%total            = energy%total          + new%total         - old%total

        ! Count succesful move
        counter%deletions = counter%deletions + 1

        ! Add the molecule to the reservoir
        if (has_reservoir) then

            ! Generate three random numbers in [0,1) and shift to [-0.5,0.5)
            call random_number(trial_pos)
            trial_pos = trial_pos - 0.5_real64

            ! Place the deleted molecule randomly in the reservoir
            reservoir%mol_com(:, residue_type, reservoir%num_residues(residue_type)+1) = &
                trial_pos(1)*reservoir%matrix(:, 1) + &
                trial_pos(2)*reservoir%matrix(:, 2) + &
                trial_pos(3)*reservoir%matrix(:, 3)
            reservoir%site_offset(:, residue_type, reservoir%num_residues(residue_type)+1, &
                1:nb%atom_in_residue(residue_type)) = &
                primary%site_offset(:, residue_type, last_molecule_index, 1:nb%atom_in_residue(residue_type))

            reservoir%num_residues(residue_type) = reservoir%num_residues(residue_type) + 1
            reservoir%num_atoms = reservoir%num_atoms + nb%atom_in_residue(residue_type)

        end if

    end subroutine AcceptDeletionMove

    subroutine RejectDeletionMove(residue_type, molecule_index, mol_com_old, site_offset_old)

        implicit none

        integer, intent(in) :: residue_type     ! Residue type to be moved
        integer, intent(in) :: molecule_index   ! Molecule ID
        real(real64), dimension(3) :: mol_com_old ! For storing old molecule center-of-mass
        real(real64), dimension(:, :) :: site_offset_old

        ! Restore previous residue/atom numbers
        primary%num_residues(residue_type) = primary%num_residues(residue_type) + 1
        primary%num_atoms = primary%num_atoms + nb%atom_in_residue(residue_type)

        ! Restore previous positions and orientation
        primary%mol_com(:, residue_type, molecule_index) = mol_com_old(:)
        primary%site_offset(:, residue_type, molecule_index, 1:nb%atom_in_residue(residue_type)) = &
            site_offset_old(:, 1:nb%atom_in_residue(residue_type))

        ! Restore Fourier states (ik_alloc and dk_alloc)
        call RestoreSingleMolFourier(residue_type, molecule_index)

    end subroutine RejectDeletionMove

end module molecule_deletion
