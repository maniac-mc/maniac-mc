program MANIAC

    use cli_utils
    use initoutput_utils
    use parameters_parser
    use tabulated_utils
    use prepare_utils
    use energy_utils
    use input_parser
    use data_parser

    implicit none

    ! Step 1 : Program initialization
    call ParseCommandLineArguments()    ! Handle -i, -d, -p, -r, -o options
    call InitOutput()                   ! Open log file and create output directory

    ! Step 2 : Read input files
    call ReadInput()                    ! Read the main MANIAC input file
    call ReadSystemData()               ! Read tpology/data file
    call ReadParameters()               ! Read simulation parameters (Lennard-Jones, etc.)

    ! Step 3 : Simulation preparation
    call PrepareSimulationParameters()  ! Set up MC parameters, initial checks
    call PrecomputeTable()              ! Precompute tables for faster calculation
    call ComputeSystemEnergy(primary)   ! Compute initial total energy

!     ! Step 4 :Monte Carlo simulation
!     call MonteCarloLoop()               ! Main MC loop

    ! Step 5 :Final reporting and cleanup
    call FinalReport()                  ! Print energy and statistics

end program MANIAC
