module cli_utils

    use output_utils 
    use, intrinsic :: iso_fortran_env, only: real64

    implicit none

contains

    !-----------------------------------------------------------------------------
    ! ParseCommandLineArguments
    !
    ! Purpose:
    !   Reads and processes command-line arguments for the MANIAC program.
    !   Recognized options:
    !     -i <input_file>      : Main MANIAC input file (mandatory)
    !     -d <data_file>       : Topology/data file (mandatory)
    !     -p <parameters_file> : Parameter include file (mandatory)
    !     -r <reservoir_file>  : Optional reservoir file
    !     -o <output_dir>      : Optional output directory (default: 'outputs/')
    !-----------------------------------------------------------------------------
    subroutine ParseCommandLineArguments()
        call SetDefaultCLIValues()
        call ReadCommandLineArgs()
        call ValidateCLIArguments()
        call NormalizeOutputPath()
    end subroutine ParseCommandLineArguments

    !-----------------------------------------------------------------
    ! Set default CLI values
    !-----------------------------------------------------------------
    subroutine SetDefaultCLIValues()
        maniac_file = ''
        data_file   = ''
        inc_file    = ''
        res_file    = ''
        output_path = 'outputs/'
    end subroutine SetDefaultCLIValues

    !-----------------------------------------------------------------
    ! Read command-line arguments into global variables
    !-----------------------------------------------------------------
    subroutine ReadCommandLineArgs()
        implicit none
        character(len=256) :: arg
        integer :: i, nargs

        nargs = command_argument_count()
        do i = 1, nargs
            call get_command_argument(i, arg)
            select case(trim(arg))
            case ('-i')
                call get_command_argument(i+1, maniac_file)
                maniac_file = trim(maniac_file)
            case ('-d')
                call get_command_argument(i+1, data_file)
                data_file = trim(data_file)
            case ('-p')
                call get_command_argument(i+1, inc_file)
                inc_file = trim(inc_file)
            case ('-r')
                call get_command_argument(i+1, res_file)
                res_file = trim(res_file)
            case ('-o')
                call get_command_argument(i+1, output_path)
                output_path = trim(output_path)
            end select
        end do
    end subroutine ReadCommandLineArgs

    !-----------------------------------------------------------------
    ! Validate CLI arguments and check file existence
    !-----------------------------------------------------------------
    subroutine ValidateCLIArguments()
        if (trim(maniac_file) == '' .or. trim(data_file) == '' .or. trim(inc_file) == '') then
            call AbortRun("Missing mandatory input arguments: -i, -d, -p required.", 1)
        end if

        if (.not. file_exists(maniac_file)) call AbortRun("Input file not found: "//trim(maniac_file), 1)
        if (.not. file_exists(data_file)) call AbortRun("Data file not found: "//trim(data_file), 1)
        if (.not. file_exists(inc_file)) call AbortRun("Parameter file not found: "//trim(inc_file), 1)

        if (trim(res_file) /= '' .and. .not. file_exists(res_file)) call AbortRun("Reservoir file not found: "//trim(res_file), 1)
    end subroutine ValidateCLIArguments

    !-----------------------------------------------------------------
    ! Ensure output path ends with a trailing slash
    !-----------------------------------------------------------------
    subroutine NormalizeOutputPath()
        if (len_trim(output_path) > 0) then
            if (output_path(len_trim(output_path):len_trim(output_path)) /= '/') then
                output_path = trim(output_path) // '/'
            end if
        end if
    end subroutine NormalizeOutputPath

    !-----------------------------------------------------------------
    ! Utility function to check if a file exists
    !
    ! Arguments:
    !   filename - character string with file name
    !
    ! Returns:
    !   .true. if the file exists, .false. otherwise
    !-----------------------------------------------------------------------------
    logical function file_exists(filename)
        character(len=*), intent(in) :: filename
        logical :: exists

        inquire(file=trim(filename), exist=exists)
        file_exists = exists
    end function file_exists

end module cli_utils
