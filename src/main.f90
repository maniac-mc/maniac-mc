program MANIAC

    use cli_utils
    use initoutput_utils
    use parameters_parser
    use prepare_utils
    use input_parser
    use data_parser

    implicit none

    !==============================
    ! Program initialization
    !==============================
    call ParseCommandLineArguments()    ! Handle -i, -d, -p, -r, -o options
    call InitOutput()                   ! Open log file and create output directory

    !==============================
    ! Read input files
    !==============================
    call ReadInput()                    ! Read the main MANIAC input file
    call ReadSystemData()               ! Read tpology/data file   
    call ReadParameters()               ! Read simulation parameters (Lennard-Jones, etc.)

    !==============================
    ! Simulation preparation
    !==============================
    call PrepareSimulationParameters()  ! Set up MC parameters, initial checks
    ! call ComputeSystemEnergy()          ! Compute initial total energy

!     !==============================
!     ! Monte Carlo simulation
!     !==============================
!     call MonteCarloLoop()               ! Main MC loop

!     !==============================
!     ! Final reporting and cleanup
!     !==============================
!     call FinalReport()                  ! Print energy and statistics
!     call CloseOutput()                  ! Close files and finalize
    
end program MANIAC
