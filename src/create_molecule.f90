module molecule_creation

    use monte_carlo_utils
    use ewald_kvectors
    use ewald_phase
    use ewald_energy
    use simulation_state
    use random_utils
    use energy_utils
    use geometry_utils
    use constants
    use, intrinsic :: iso_fortran_env, only: real64

    implicit none

    private ! Hide everything by default
    public :: CreateMolecule  ! Only expose the main entry point

contains

    !---------------------------------------------------------------------------
    ! Subroutine: CreateMolecule
    !
    ! Purpose:
    !   Attempts to insert (create) a new molecule of a given residue type
    !   at a random position and orientation inside the simulation box,
    !   following the Monte Carlo insertion move scheme.
    !   The move is accepted or rejected based on the Metropolis criterion
    !   using the energy difference and the fugacity of the species.
    !
    !---------------------------------------------------------------------------
    subroutine CreateMolecule(residue_type, molecule_index)

        implicit none

        ! Input arguments
        integer, intent(in) :: residue_type     ! Residue type to be moved
        integer, intent(in) :: molecule_index   ! Molecule ID

        ! Local variables
        integer :: rand_mol_index               ! Randomly selected molecule index from the reservoir for copying geometry
        real(real64) :: probability             ! Acceptance probability of creation move
        logical :: is_creation                  ! Flag indicating creation

        call CheckMoleculeIndex(molecule_index, NB_MAX_MOLECULE)
        
        ! Count trial move (success + fail)
        counter%trial_creations = counter%trial_creations + 1

        ! Compute old energy
        is_creation = .true.
        call ComputeOldEnergy(residue_type, molecule_index, old, is_creation = is_creation)

        ! Increase the residue and atom counts
        primary%num_residues(residue_type) = primary%num_residues(residue_type) + 1
        primary%num_atoms = primary%num_atoms + nb%atom_in_residue(residue_type)

        ! Save current Fourier terms (should be all zeros here)
        call SaveSingleMolFourierTerms(residue_type, molecule_index) 

        ! Generate random insertion position within the simulation box
        call InsertAndOrientMolecule(residue_type, molecule_index, rand_mol_index)

        ! Compute new energy
        call ComputeNewEnergy(residue_type, molecule_index, new, is_creation = is_creation)

        ! Compute acceptance probability for the move
        probability = mc_acceptance_probability(old, new, residue_type, TYPE_CREATION)

        ! Accept or reject
        if (rand_uniform() <= probability) then ! Accept move
            call AcceptCreationMove(residue_type, rand_mol_index, old, new)
        else ! Reject move
            call RejectCreationMove(residue_type, molecule_index)
        end if

    end subroutine CreateMolecule

    subroutine AcceptCreationMove(residue_type, rand_mol_index, old, new)

        implicit none

        integer, intent(in) :: residue_type     ! Residue type to be moved
        integer :: rand_mol_index               ! Randomly selected molecule index from the reservoir for copying geometry
        type(energy_state), intent(in) :: old   ! Previous energy states
        type(energy_state), intent(in) :: new   ! New energy states
        integer :: last_molecule_index          ! Index of the last molecule in the reservoir (used when removing a molecule)

        ! Update total energies
        energy%recip_coulomb    = new%recip_coulomb
        energy%non_coulomb      = energy%non_coulomb    + new%non_coulomb   - old%non_coulomb
        energy%coulomb          = energy%coulomb        + new%coulomb       - old%coulomb
        energy%ewald_self       = energy%ewald_self     + new%ewald_self    - old%ewald_self
        energy%intra_coulomb    = energy%intra_coulomb  + new%intra_coulomb - old%intra_coulomb
        energy%total            = energy%total          + new%total         - old%total 

        ! Count successful move
        counter%creations = counter%creations + 1

        ! Remove molecule from reservoir if present
        if (has_reservoir) then

            ! Replace molecule_index with the last molecule in the list to maintain continuity
            last_molecule_index = reservoir%num_residues(residue_type)
            reservoir%mol_com(:, residue_type, rand_mol_index) = reservoir%mol_com(:, residue_type, last_molecule_index)
            reservoir%site_offset(:, residue_type, rand_mol_index, 1:nb%atom_in_residue(residue_type)) = &
                reservoir%site_offset(:, residue_type, last_molecule_index, 1:nb%atom_in_residue(residue_type))

            reservoir%num_residues(residue_type) = reservoir%num_residues(residue_type) - 1
            reservoir%num_atoms = reservoir%num_atoms - nb%atom_in_residue(residue_type)
        
        end if

    end subroutine AcceptCreationMove

    subroutine RejectCreationMove(residue_type, molecule_index)

        implicit none

        integer, intent(in) :: residue_type     ! Residue type to be moved
        integer, intent(in) :: molecule_index   ! Molecule ID

        ! Restore previous residue/atom numbers
        primary%num_atoms = primary%num_atoms - nb%atom_in_residue(residue_type)
        primary%num_residues(residue_type) = primary%num_residues(residue_type) - 1

        ! Restore Fourier states (ik_alloc and dk_alloc, all zeros)
        call RestoreSingleMolFourier(residue_type, molecule_index)

    end subroutine RejectCreationMove

    subroutine InsertAndOrientMolecule(residue_type, molecule_index, rand_mol_index)

        implicit none

        ! Input arguments
        integer, intent(in) :: residue_type     ! Residue type to be moved
        integer, intent(in) :: molecule_index   ! Molecule ID
        integer :: rand_mol_index               ! Randomly selected molecule index from the reservoir for copying geometry
        ! Local variables
        logical :: full_rotation                ! Flag indicating whether a full 360° random rotation should be applied
        real(real64) :: random_nmb              ! Uniform random number in [0,1), used for random index selection
        real(real64) :: trial_pos(3)            ! Random numbers for initial molecule position in the box

        ! Generate a random position in the simulation box
        call random_number(trial_pos) ! Random numbers in [0,1)
        primary%mol_com(:, residue_type, molecule_index) = primary%bounds(:,1) &
            + matmul(primary%matrix, trial_pos)

        ! Copy geometry from reservoir or rotate if no reservoir
        if (has_reservoir) then

            ! Pick a random (and existing) molecule in the reservoir
            call random_number(random_nmb) ! generates random_nmb in [0,1)
            rand_mol_index = int(random_nmb * reservoir%num_residues(residue_type)) + 1 ! random integer in [1, N]

            ! Copy site offsets from the chosen molecule
            primary%site_offset(:, residue_type, molecule_index, 1:nb%atom_in_residue(residue_type)) = &
                reservoir%site_offset(:, residue_type, rand_mol_index, 1:nb%atom_in_residue(residue_type))

        else

            ! Copy site offsets from the first molecule
            primary%site_offset(:, residue_type, molecule_index, 1:nb%atom_in_residue(residue_type)) = &
                primary%site_offset(:, residue_type, 1, 1:nb%atom_in_residue(residue_type))

            ! Rotate the new molecule randomly (using full 360° rotation)
            full_rotation = .True.
            call ApplyRandomRotation(residue_type, molecule_index, full_rotation)
        end if

    end subroutine InsertAndOrientMolecule

end module molecule_creation