module initoutput_utils

    use simulation_state
    use output_utils
    use helper_utils

    implicit none

contains

    !---------------------------------------------------------------------------
    ! Subroutine: OpenOutput
    ! Purpose   : Open the main output file for writing, creating the output
    !             directory if necessary. Reports errors if operations fail.
    !---------------------------------------------------------------------------
    subroutine InitOutput()

        implicit none

        ! Ensure output directory exists
        call EnsureDirectoryExists(output_path)

        ! Open log file for writing, replacing any existing file
        call OpenLogFile(output_path)

        ! Write header to log
        call WriteHeader()

    end subroutine InitOutput

    !-------------------------------------------------------
    ! Ensure directory exists (create if missing)
    !-------------------------------------------------------
    subroutine EnsureDirectoryExists(path)
    
        character(len=*), intent(in) :: path
        character(len=200) :: command
        integer :: exit_status

        command = 'mkdir -p ' // trim(path)
        call execute_command_line(command, exitstat=exit_status)
        if (exit_status /= 0) then
            call AbortRun("Failed to create output directory: "//trim(path), exit_status)
        end if
    end subroutine EnsureDirectoryExists

    !-------------------------------------------------------
    ! Open log file
    !-------------------------------------------------------
    subroutine OpenLogFile(path)
        character(len=*), intent(in) :: path
        integer :: ios

        open(unit=out_unit, file=trim(path)//'log.maniac', status='replace', iostat=ios)
        if (ios /= 0) call AbortRun("Failed to open log file: "//trim(path)//'log.maniac', ios)
    end subroutine OpenLogFile

    !---------------------------------------------------------------------------
    ! Subroutine: WriteHeader
    ! Purpose   : Write a modern ASCII header with project info (dynamic version)
    !---------------------------------------------------------------------------
    subroutine WriteHeader()

        use version_module
       
        implicit none
        integer, parameter :: box_width = 78 ! width of the box

        ! Blank line before header
        call LogMessage("")

        ! Top border
        call LogMessage("+" // repeat_char("-", box_width-2) // "+")

        ! Version info lines
        call BoxLine("MANIAC-MC - Version " // version, box_width)
        call BoxLine("Git commit: " // git_commit, box_width)

        ! Credits (optional, you can also keep them in version_module if desired)
        call BoxLine("Code written and maintained by Simon Gravelle, LIPhy, CNRS", box_width)

        ! Bottom border
        call LogMessage("+" // repeat_char("-", box_width-2) // "+")

        ! Blank line after header
        call LogMessage("")

    end subroutine WriteHeader

end module initoutput_utils
