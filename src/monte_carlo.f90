module montecarlo_module

    use molecule_translation    ! Provides routines to perform translational moves on molecules
    use monte_carlo_utils       ! Provides Monte Carlo-specific utilities (e.g., step adjustment)
    use molecule_creation       ! Provides routines to insert new molecules into the system
    use molecule_deletion       ! Provides routines to remove molecules from the system
    use molecule_rotation       ! Provides routines to perform rotational moves on molecules
    use simulation_state        ! Stores and manages the current state of the simulation
    use output_utils            ! Handles updating files and printing simulation status
    use random_utils            ! Provides random number generation routines
    use write_utils             ! Utilities for writing formatted output or logs
    use constants               ! Defines physical and simulation constants

    use, intrinsic :: iso_fortran_env, only: real64 ! Ensures consistent 64-bit real precision

    implicit none

contains

    !------------------------------------------------------------------------------
    ! subroutine MonteCarloLoop
    !
    ! Performs the main Monte Carlo simulation loop, executing a series of molecular
    ! moves and exchanges over multiple blocks and steps.
    !------------------------------------------------------------------------------
    subroutine MonteCarloLoop()

        implicit none

        ! 1. Local Variables
        real(real64) :: random_draw     ! Random number used for move selection
        integer :: residue_type         ! Index of molecule to be moved
        integer :: molecule_index       ! Index of molecule copy

        ! 2. Initialization
        call LogStartMC()               ! Log starting message
        call UpdateFiles(.false.)       ! Write initial topology

        ! 3. Main Monte Carlo Loop
        do current_block = 1, input%nb_block        ! Loop over blocks
            do current_step = 1, input%nb_step      ! Loop over steps within each block

                ! Pick a molecule type and instance
                residue_type = PickRandomResidueType(input%is_active)
                molecule_index = PickRandomMoleculeIndex(primary%num_residues(residue_type))

                ! Perform Monte Carlo move
                random_draw = rand_uniform()

                if (random_draw <= proba%translation) then

                    ! Case 1: Small translation move
                    call Translation(residue_type, molecule_index)

                else if (random_draw <= proba%rotation+proba%translation) then

                    ! Case 2: Rotation move
                    call Rotation(residue_type, molecule_index)

                else ! Insertion/deletion move

                    ! Case 3: Insertion or deletion move
                    if (rand_uniform() <= PROB_CREATE_DELETE) then

                        ! Attempt to create a molecule of species residue_type
                        molecule_index = primary%num_residues(residue_type) + 1
                        call CreateMolecule(residue_type, molecule_index)

                    else

                        ! Attempt to delete a randomly chosen molecule of species residue_type
                        call DeleteMolecule(residue_type, molecule_index)

                    end if
                end if

             end do

            !----------------------------------------------
            ! Adjust Monte Carlo step sizes & output status
            !----------------------------------------------
            call AdjustMoveStepSizes()          ! Adjust the Monte Carlo move step sizes based on recent acceptance ratios
            call PrintStatus()                  ! Print the current status of the simulation to the log, data, and trajectory files
            call UpdateFiles(.true.)

        end do

    end subroutine MonteCarloLoop

end module montecarlo_module
