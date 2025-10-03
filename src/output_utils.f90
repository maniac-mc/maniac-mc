!-----------------------------------------------------------------------------
! module output_utils
!
! This module handles all output-related operations for the Monte Carlo
! simulation.
!-----------------------------------------------------------------------------

module output_utils

    use simulation_state
    use helper_utils
    use, intrinsic :: iso_fortran_env, only: real64

    implicit none

contains

    !---------------------------------------------------------------------------
    ! Subroutine: CloseOutput
    ! Purpose   : Close the main output file.
    !---------------------------------------------------------------------------
    subroutine CloseOutput()

        implicit none

        call PrintTerminationMessage()

        close(out_unit)

    end subroutine CloseOutput

    !---------------------------------------------------------------------------
    ! Subroutine: LogMessage
    ! Purpose   : Write a custom message to the specified output unit.
    !---------------------------------------------------------------------------
    subroutine LogMessage(msg)

        implicit none

        character(*), intent(in) :: msg

        write(*,*) TRIM(msg)
        write(out_unit,*) TRIM(msg)
        flush(out_unit) ! Forces the data to be written to disk

    end subroutine LogMessage

    ! Subroutine: LogStartMC
    ! Purpose   : Log the start of the Monte Carlo loop with a modern ASCII style
    subroutine LogStartMC()
        implicit none
        integer, parameter :: box_width = 78

        ! Blank line before message
        call LogMessage("")

        ! Top border
        call LogMessage("+" // repeat_char("=", box_width-2) // "+")

        ! Message line centered
        call BoxLine("Started Monte Carlo Loop", box_width)

        ! Bottom border
        call LogMessage("+" // repeat_char("=", box_width-2) // "+")

        ! Blank line after message
        call LogMessage("")

    end subroutine LogStartMC

    !-----------------------------------------------------------------------
    ! Helper: BoxLine
    ! Writes a line inside the ASCII box, padded to box_width
    !-----------------------------------------------------------------------
    subroutine BoxLine(text, width)
        character(len=*), intent(in) :: text
        integer, intent(in) :: width
        character(len=width-4) :: padded

        padded = adjustl(text)
        if (len_trim(padded) < width-4) then
            padded(len_trim(padded)+1:) = ' '  ! pad with spaces
        end if
        call LogMessage("| " // padded(1:width-4) // " |")
    end subroutine BoxLine

    !---------------------------------------------------------------------------
    ! Subroutine: PrintTerminationMessage
    ! Purpose   : Print a formatted termination footer with Monte Carlo summary
    !---------------------------------------------------------------------------
    subroutine PrintTerminationMessage()
        use simulation_state       ! Access counters and output path
        implicit none

        integer, parameter :: box_width = 78
        character(len=256) :: line

        ! Blank line before footer
        call LogMessage("")

        ! Top border
        call LogMessage("+" // repeat_char("-", box_width-2) // "+")

        ! Title
        call BoxLine("MANIAC-MC Simulation Completed", box_width)
        call BoxLine("", box_width)  ! blank line inside box

        ! Summary statistics (Trial / Accepted moves)
        write(line,'(A,I8,A,I8)') "  Translations (Trial/Accepted): ", &
            counter%trial_translations, " / ", counter%translations
        call BoxLine(trim(line), box_width)

        write(line,'(A,I8,A,I8)') "  Rotations    (Trial/Accepted): ", &
            counter%trial_rotations, " / ", counter%rotations
        call BoxLine(trim(line), box_width)

        write(line,'(A,I8,A,I8)') "  Big moves    (Trial/Accepted): ", &
            counter%trial_bigmoves, " / ", counter%big_moves
        call BoxLine(trim(line), box_width)

        write(line,'(A,I8,A,I8)') "  Creations    (Trial/Accepted): ", &
            counter%trial_creations, " / ", counter%creations
        call BoxLine(trim(line), box_width)

        write(line,'(A,I8,A,I8)') "  Deletions    (Trial/Accepted): ", &
            counter%trial_deletions, " / ", counter%deletions
        call BoxLine(trim(line), box_width)

        call BoxLine("", box_width)  ! blank line inside box

        ! Output path information
        call BoxLine("All output files have been written to:", box_width)
        call BoxLine(trim(output_path), box_width)

        ! Bottom border
        call LogMessage("+" // repeat_char("-", box_width-2) // "+")

        ! Blank line after footer
        call LogMessage("")

    end subroutine PrintTerminationMessage

    !------------------------------------------------------------------------------
    ! Subroutine: PrintStatus
    !
    ! Purpose:
    !   Generates a formatted summary of the current Monte Carlo simulation state
    !   and prints it to the log file via the LogMessage routine.
    !
    ! Notes:
    !   - Energies are converted to kcal/mol for reporting.
    !------------------------------------------------------------------------------
    subroutine PrintStatus()

        implicit none

        integer :: nb_type_residue
        real(real64) :: e_tot_kcalmol
        real(real64) :: e_coul_kcalmol
        real(real64) :: e_long_kcalmol
        real(real64) :: e_recip_coulomb_kcalmol
        real(real64) :: e_non_coulomb_kcalmol
        real(real64) :: e_coulomb_kcalmol
        real(real64) :: e_ewald_self_kcalmol
        real(real64) :: e_intra_coulomb_kcalmol
        real(real64) :: e_total_kcalmol
        character(len=64) :: tmp
        character(LEN=1024) :: formatted_msg

        call LogMessage("")

        formatted_msg = "  Energy report | Active molecules: "
        do nb_type_residue = 1, nb%type_residue
            if (primary%num_residues(nb_type_residue) /= 0 .and. input%is_active(nb_type_residue) == 1) then
                write(tmp,'(A,"=",I0)') trim(res%names_1d(nb_type_residue)), primary%num_residues(nb_type_residue)
                if (len_trim(formatted_msg) > 0) then
                    formatted_msg = trim(formatted_msg)//" " // trim(tmp)
                else
                    formatted_msg = trim(tmp)
                end if
            end if
        end do
        call LogMessage(formatted_msg)

        ! Convert energies in kcal/mol
        e_recip_coulomb_kcalmol = energy%recip_coulomb * KB_kcalmol
        e_non_coulomb_kcalmol = energy%non_coulomb * KB_kcalmol
        e_coulomb_kcalmol = energy%coulomb * KB_kcalmol
        e_ewald_self_kcalmol = energy%ewald_self * KB_kcalmol
        e_intra_coulomb_kcalmol = energy%intra_coulomb * KB_kcalmol
        e_total_kcalmol = energy%total * KB_kcalmol

        ! Compute total energies
        e_tot_kcalmol = e_non_coulomb_kcalmol + e_recip_coulomb_kcalmol + e_coulomb_kcalmol + &
                e_ewald_self_kcalmol + e_intra_coulomb_kcalmol
        e_coul_kcalmol = e_coulomb_kcalmol + e_intra_coulomb_kcalmol
        e_long_kcalmol = e_recip_coulomb_kcalmol + e_ewald_self_kcalmol

        write(formatted_msg,'(A10,1X,A14,1X,A14,1X,A14,1X,A14,2X,A10,2X,A10,2X,A20)') &
            'Step','TotEng','E_vdwl','E_coul','E_long','TransStep','RotAngle','MC (acc/trial)'

        call LogMessage(formatted_msg)

        write(formatted_msg,'(I10,1X,F14.4,1X,F14.4,1X,F14.4,1X,F14.4,2X,F10.4,2X,F10.4,2X,'// &
                            ' "T(",I0,"/",I0,") R(",I0,"/",I0,") B(",I0,"/",I0,") C(",I0,"/",I0,") D(",I0,"/",I0,")")') &
            current_block, e_tot_kcalmol, e_non_coulomb_kcalmol, e_coul_kcalmol, e_long_kcalmol, &
            input%translation_step, input%rotation_step_angle, &
            counter%translations, counter%trial_translations, &
            counter%rotations, counter%trial_rotations, &
            counter%big_moves, counter%trial_bigmoves, &
            counter%creations, counter%trial_creations, &
            counter%deletions, counter%trial_deletions
        call LogMessage(formatted_msg)

    end subroutine PrintStatus

    !------------------------------------------------------------------------------
    ! Subroutine: FinalReport
    !------------------------------------------------------------------------------
    subroutine FinalReport()

        implicit none

        ! Energy components in kcal/mol
        real(real64) :: e_tot          ! Total energy for reporting (computed)
        real(real64) :: e_coul         ! Coulombic energy including intra-molecular interactions
        real(real64) :: e_long         ! Long-range Coulombic energy (reciprocal + self)
        real(real64) :: e_recip ! Reciprocal-space Coulombic energy
        real(real64) :: e_non_coulomb  ! Non-Coulombic energy (van der Waals, etc.)
        real(real64) :: e_coulomb       ! Direct-space Coulombic energy
        real(real64) :: e_self    ! Ewald self-energy correction
        real(real64) :: e_intra ! Intra-molecular Coulombic energy
        character(LEN=1024) :: formatted_msg   ! Formatted message for logging
        integer, parameter :: box_width = 78

        ! Convert energies to kcal/mol
        e_recip   = energy%recip_coulomb * KB_kcalmol
        e_non_coulomb     = energy%non_coulomb * KB_kcalmol
        e_coulomb = energy%coulomb * KB_kcalmol
        e_self    = energy%ewald_self * KB_kcalmol
        e_intra   = energy%intra_coulomb * KB_kcalmol

        ! Compute combined components
        e_tot  = e_non_coulomb + e_recip + e_coulomb + e_self + e_intra
        e_coul = e_coulomb + e_intra
        e_long = e_recip + e_self

        ! Blank line before box
        call LogMessage("")

        ! Top border
        call LogMessage("+" // repeat_char("-", box_width-2) // "+")

        ! Box title
        call BoxLine("Final Energy Report", box_width)
        call BoxLine("", box_width)

        ! Column headers
        call BoxLine("  Step        TotEng        E_vdwl        E_coul        E_long", box_width)

        ! Energies line
        write(formatted_msg,'(I10,1X,F15.6,1X,F15.6,1X,F15.6,1X,F15.6)') &
            current_block, e_tot, e_non_coulomb, e_coul, e_long
        call BoxLine(trim(formatted_msg), box_width)

        call BoxLine("", box_width)  ! blank line inside box

        ! Bottom border
        call LogMessage("+" // repeat_char("-", box_width-2) // "+")

        ! Blank line after box
        call LogMessage("")

        call CloseOutput()                  ! Close files and finalize

    end subroutine FinalReport

    subroutine LogParameters(input_file_name, n_pairs, pair1, pair2, epsilons, sigmas)

        implicit none

        character(len=*), intent(in) :: input_file_name         ! Name of the parameter input file being processed
        integer, intent(out) :: n_pairs                         ! Counter for the number of unique atom type pairs
        integer, intent(out) :: pair1(:), pair2(:)              ! Arrays to store the first and second atom type indices of unique pairs
        real(real64), intent(out) :: epsilons(:), sigmas(:)     ! Arrays to store epsilon (kcal/mol) and sigma (Å) for unique pairs

        integer :: i, j, k, l, m, n                             ! Loop indices for iterating over residues and atom types
        integer :: type1, type2, type_m, type_n                 ! Atom type indices for site and atom pairs
        integer :: it1, it2                                     ! Sorted atom type indices (min and max) for unique pair identification
        logical :: found_pair                                   ! Flag to indicate if a pair was found (unused but set for clarity)
        character(len=200) :: formatted_msg                     ! String to hold formatted log messages for output
        logical, allocatable :: printed(:,:)                    ! 2D array to track which atom type pairs have been recorded
        integer :: max_atom_type                                ! Maximum atom type index from atom_types_2d array

        if (has_reservoir) then
            max_atom_type = max(maxval(primary%atom_types(:,:)), maxval(reservoir%atom_types(:,:)) )
        else
            max_atom_type = maxval(primary%atom_types(:,:))
        end if
        allocate(printed(max_atom_type, max_atom_type))
        printed = .false.

        ! Log the start of parameter file import
        call LogMessage("")
        call LogMessage("====== Import parameter file ======")
        call LogMessage("")
        write(formatted_msg, '("Reading file ", A)') TRIM(input_file_name) ! Format message with input file name
        call LogMessage(formatted_msg)                                     ! Log the input file name

        n_pairs = 0 ! Initialize the pair counter
        do i = 1, nb%type_residue
            do j = 1, nb%types_per_residue(i)
                do k = 1, nb%type_residue
                    do l = 1, nb%types_per_residue(k)
                        type1 = res%types_2d(i, j) ! Get atom type for site j in residue i
                        type2 = res%types_2d(k, l) ! Get atom type for site l in residue k
                        found_pair = .false.        ! Reset found_pair flag (unused in logic)
                        if (type1 <= type2) then    ! Process pairs where type1 <= type2 to avoid duplicates
                            do m = 1, nb%atom_in_residue(i)
                                do n = 1, nb%atom_in_residue(k)
                                    type_m = primary%atom_types(i, m) ! Get atom type for atom m in residue i
                                    type_n = primary%atom_types(k, n) ! Get atom type for atom n in residue k
                                    it1 = min(type_m, type_n)    ! Get minimum atom type for pair
                                    it2 = max(type_m, type_n)    ! Get maximum atom type for pair
                                    if (it1 > 0 .and. it2 > 0) then
                                        if (.not. printed(it1, it2)) then ! Check if pair hasn’t been recorded
                                            n_pairs = n_pairs + 1         ! Increment pair counter
                                            pair1(n_pairs) = it1          ! Store first atom type
                                            pair2(n_pairs) = it2          ! Store second atom type
                                            epsilons(n_pairs) = coeff%epsilon(i,k,m,n) * KB_kcalmol ! Store epsilon in kcal/mol
                                            sigmas(n_pairs) = coeff%sigma(i,k,m,n)                  ! Store sigma in Å
                                            printed(it1, it2) = .true.    ! Mark pair as recorded
                                            printed(it2, it1) = .true.    ! Mark symmetric pair as recorded
                                        end if
                                    end if
                                end do
                            end do
                        end if
                    end do
                end do
            end do
        end do
    end subroutine LogParameters

    subroutine LogData(data_file_name, box, is_primary)

        implicit none

        ! Input parameters
        type(type_box), intent(inout) :: box
        logical, intent(in) :: is_primary

        ! Local variables
        character(len=*), intent(in) :: data_file_name
        character(len=200) :: formatted_msg
        integer :: i, atom_id
        integer :: active_molecule_count ! Number of species that will be moved/inserted/deleted

        ! === Step 9: Logging ===
        call LogMessage("")
        call LogMessage("====== Import data file ======")
        write(formatted_msg, '("Reading file ", A)') TRIM(data_file_name)
        call LogMessage(formatted_msg)
        call LogMessage("")
        write(formatted_msg, '("Number of atoms: ", I0)') box%num_atoms
        call LogMessage(formatted_msg)
        write(formatted_msg, '("Number of type of residues: ", I0)') nb%type_residue
        call LogMessage(formatted_msg)
        write(formatted_msg, '("Number of type of atoms: ", I0)') box%num_atomtypes
        call LogMessage(formatted_msg)

        do i = 1, nb%type_residue
            if ((box%num_residues(i) /= 0) .and. (input%is_active(i) == 1)) then
                ! Active residue present in data file
                active_molecule_count = active_molecule_count + box%num_residues(i)
                write(formatted_msg, '("Active residue ", A, " found in the data file: ", I0)') &
                    TRIM(res%names_1d(i)), box%num_residues(i)
            else if ((box%num_residues(i) /= 0) .and. (input%is_active(i) == 0)) then
                ! Inactive residue present in data file
                active_molecule_count = active_molecule_count + box%num_residues(i)
                write(formatted_msg, '("Inactive residue ", A, " found in the data file: ", I0)') &
                    TRIM(res%names_1d(i)), box%num_residues(i)
            else if ((box%num_residues(i) == 0) .and. (input%is_active(i) == 0) .and. (is_primary)) then
                ! Inactive residue defined in input but not present in data file
                call AbortRun("Inactive residue '" // trim(res%names_1d(i)) // "' (ID=" // &
                            trim(adjustl(to_string(i))) // ") defined in input file but not present in data file.", 1)
            end if
            call LogMessage(formatted_msg)
        end do

        call LogMessage("")
        call LogMessage("Simulation box (rows):")
        write(formatted_msg, '(3F12.6)') box%matrix(1,1), box%matrix(1,2), box%matrix(1,3)
        call LogMessage(formatted_msg)
        write(formatted_msg, '(3F12.6)') box%matrix(2,1), box%matrix(2,2), box%matrix(2,3)
        call LogMessage(formatted_msg)
        write(formatted_msg, '(3F12.6)') box%matrix(3,1), box%matrix(3,2), box%matrix(3,3)
        call LogMessage(formatted_msg)

        call LogMessage("")
        call LogMessage("Atoms masses (g/mol):")
        do atom_id = 1, primary%num_atomtypes
            write(formatted_msg, '(I5, 2X, F12.6)') atom_id, box%site_masses_vector(atom_id)
            call LogMessage(formatted_msg)
        end do

    end subroutine LogData

    subroutine LogConnectivity(box)

        implicit none

        ! Input
        type(type_box), intent(inout) :: box

        ! Local
        character(len=200) :: formatted_msg
        integer :: i, j
        integer, parameter :: MAX_PRINT = 6

        if ((box%num_bonds > 0) .or. (box%num_angles > 0)) then

            call LogMessage("")
            call LogMessage("===== Connectivity summary =====")

            ! --- Bonds ---
            call LogMessage("")
            do i = 1, nb%type_residue
                if (box%num_residues(i) > 0) then
                    write(formatted_msg, '("Residue ", A, ": ", I0, " bonds")') &
                        trim(res%names_1d(i)), nb%bonds_per_residue(i)
                    call LogMessage(formatted_msg)

                    ! Print up to MAX_PRINT bonds
                    do j = 1, min(nb%bonds_per_residue(i), MAX_PRINT)
                        write(formatted_msg, '("   bond type ", I0, ": atoms [", I0, ",", I0, "]")') &
                            res%bond_type_2d(i,j,1), res%bond_type_2d(i,j,2), res%bond_type_2d(i,j,3)
                        call LogMessage(formatted_msg)
                    end do

                    ! If more bonds exist, indicate truncation
                    if (nb%bonds_per_residue(i) > MAX_PRINT) then
                        write(formatted_msg, '("   ... ", I0, " more bonds not shown")') &
                            nb%bonds_per_residue(i) - MAX_PRINT
                        call LogMessage(formatted_msg)
                    end if
                end if
            end do

            ! --- Angles ---
            call LogMessage("")
            do i = 1, nb%type_residue

                if (box%num_residues(i) > 0) then
                    write(formatted_msg, '("Residue ", A, ": ", I0, " angles")') &
                        trim(res%names_1d(i)), nb%angles_per_residue(i)
                    call LogMessage(formatted_msg)

                    ! Print up to MAX_PRINT angles
                    do j = 1, min(nb%angles_per_residue(i), MAX_PRINT)
                        write(formatted_msg, '("   angle type ", I0, ": atoms [", I0, ",", I0, ",", I0, "]")') &
                            res%angle_type_2d(i,j,1), res%angle_type_2d(i,j,2), &
                            res%angle_type_2d(i,j,3), res%angle_type_2d(i,j,4)
                        call LogMessage(formatted_msg)
                    end do

                    ! If more angles exist, indicate truncation
                    if (nb%angles_per_residue(i) > MAX_PRINT) then
                        write(formatted_msg, '("   ... ", I0, " more angles not shown")') &
                            nb%angles_per_residue(i) - MAX_PRINT
                        call LogMessage(formatted_msg)
                    end if
                end if
            end do

            ! --- Dihedrals ---
            call LogMessage("")
            do i = 1, nb%type_residue

                if (box%num_residues(i) > 0) then
                    write(formatted_msg, '("Residue ", A, ": ", I0, " dihedrals")') &
                        trim(res%names_1d(i)), nb%dihedrals_per_residue(i)
                    call LogMessage(formatted_msg)

                    ! Print up to MAX_PRINT dihedrals
                    do j = 1, min(nb%dihedrals_per_residue(i), MAX_PRINT)
                        write(formatted_msg, &
                        '("   dihedral type ", I0, ": atoms [", I0, ",", I0, ",", I0, ",", I0, "]")') &
                            res%dihedral_type_2d(i,j,1), res%dihedral_type_2d(i,j,2), &
                            res%dihedral_type_2d(i,j,3), res%dihedral_type_2d(i,j,4), &
                            res%dihedral_type_2d(i,j,5)
                        call LogMessage(formatted_msg)
                    end do

                    ! If more dihedrals exist, indicate truncation
                    if (nb%dihedrals_per_residue(i) > MAX_PRINT) then
                        write(formatted_msg, '("   ... ", I0, " more dihedrals not shown")') &
                            nb%dihedrals_per_residue(i) - MAX_PRINT
                        call LogMessage(formatted_msg)
                    end if
                end if
            end do

            ! --- Impropers ---
            call LogMessage("")
            do i = 1, nb%type_residue

                if (box%num_residues(i) > 0) then
                    write(formatted_msg, '("Residue ", A, ": ", I0, " impropers")') &
                        trim(res%names_1d(i)), nb%impropers_per_residue(i)
                    call LogMessage(formatted_msg)

                    ! Print up to MAX_PRINT impropers
                    do j = 1, min(nb%impropers_per_residue(i), MAX_PRINT)
                        write(formatted_msg, &
                        '("   improper type ", I0, ": atoms [", I0, ",", I0, ",", I0, ",", I0, "]")') &
                            res%improper_type_2d(i,j,1), res%improper_type_2d(i,j,2), &
                            res%improper_type_2d(i,j,3), res%improper_type_2d(i,j,4), &
                            res%improper_type_2d(i,j,5)
                        call LogMessage(formatted_msg)
                    end do

                    ! If more impropers exist, indicate truncation
                    if (nb%impropers_per_residue(i) > MAX_PRINT) then
                        write(formatted_msg, '("   ... ", I0, " more impropers not shown")') &
                            nb%impropers_per_residue(i) - MAX_PRINT
                        call LogMessage(formatted_msg)
                    end if
                end if
            end do

        end if

    end subroutine LogConnectivity

    subroutine AbortRun(error_msg, exit_code)

        use, intrinsic :: iso_fortran_env, only: output_unit

        implicit none

        character(len=*), intent(in) :: error_msg   ! Error description
        integer, intent(in), optional :: exit_code  ! Exit code (default = 1)

        integer :: code

        if (present(exit_code)) then
            code = exit_code
        else
            code = 1
        end if

        call LogMessage("--------------------------------------------------")
        call LogMessage("FATAL ERROR:")
        call LogMessage("" // trim(error_msg))
        call LogMessage("Simulation will now terminate.")
        call LogMessage("--------------------------------------------------")

        ! Ensure buffers/files are flushed if LogMessage writes to files
        call flush(output_unit)

        stop code
    end subroutine AbortRun

    !===============================================================
    ! Subroutine: WarnUser
    ! Purpose:
    !   Print a standardized warning message without terminating
    !===============================================================
    subroutine WarnUser(warn_msg)

        use, intrinsic :: iso_fortran_env, only: output_unit
        implicit none

        character(len=*), intent(in) :: warn_msg   ! Warning description

        call LogMessage("--------------------------------------------------")
        call LogMessage("WARNING:")
        call LogMessage("" // trim(warn_msg))
        call LogMessage("Execution will continue.")
        call LogMessage("--------------------------------------------------")

        ! Ensure message is flushed immediately
        call flush(output_unit)

    end subroutine WarnUser

    !===============================================================
    ! Subroutine: InfoMessage
    ! Purpose:
    !   Print a standardized informational message
    !===============================================================
    subroutine InfoMessage(info_msg)

        use, intrinsic :: iso_fortran_env, only: output_unit
        implicit none

        character(len=*), intent(in) :: info_msg   ! Information description

        call LogMessage("INFO: " // trim(info_msg))
        call flush(output_unit)

    end subroutine InfoMessage

    subroutine CheckMoleculeIndex(molecule_index, NB_MAX_MOLECULE)

        implicit none

        ! Input arguments
        integer, intent(in) :: molecule_index       ! Index of molecule being created
        integer, intent(in) :: NB_MAX_MOLECULE      ! Maximum allowed number of molecules

        ! Local
        character(len=200) :: formatted_msg         ! Buffer for constructing formatted log or error messages

        if (molecule_index > NB_MAX_MOLECULE) then
            write(formatted_msg, '(A, I0)') &
                "Trying to insert a molecule with index = " , molecule_index
            call AbortRun(trim(formatted_msg)//new_line('a')// &
                        "This exceeds the maximum allowed number of molecules = "// &
                        trim(adjustl(to_string(NB_MAX_MOLECULE)))//new_line('a')// &
                        "Increase 'NB_MAX_MOLECULE' in src/parameters.f90 and recompile.")
        end if

    end subroutine CheckMoleculeIndex

    subroutine ValidateInputTranslation(step, move_type)

        implicit none

        ! Input arguments
        real(real64), intent(in) :: step         ! Translation step
        integer, intent(in) :: move_type         ! Move type: 0 or 1

        ! Local
        character(len=200) :: formatted_msg

        ! Check translation step
        if (step <= 0.0_real64) then
            write(formatted_msg, '(A,F6.3)') "ERROR: translation_step must be positive. Provided value: ", step
            call AbortRun(trim(formatted_msg))
        end if

        ! Check move_type
        if (move_type /= TYPE_SMALLMOVE .and. move_type /= TYPE_BIGMOVE) then
            write(formatted_msg, '(A,I0)') "ERROR: move_type must be 0 or 1. Provided value: ", move_type
            call AbortRun(trim(formatted_msg))
        end if

    end subroutine ValidateInputTranslation

    !-----------------------------------------------------------------------------
    ! Function: to_string
    !
    ! Purpose:
    !   Converts an integer input to a trimmed character string.
    !   Returns an allocatable string of minimal length.
    !-----------------------------------------------------------------------------
    pure function to_string(i) result(str)

        implicit none

        integer, intent(in) :: i
        character(len=:), allocatable :: str
        character(len=32) :: tmp                ! big enough for most integers

        write(tmp, '(I0)') i
        str = trim(tmp)

    end function to_string

    !-----------------------------------------------------------------------------
    ! Subroutine: PrintInputSummary
    !
    ! Purpose:
    !   Logs a summary of the parsed input, including global simulation parameters,
    !   Monte Carlo settings, and detailed residue information.
    !-----------------------------------------------------------------------------
    subroutine PrintInputSummary()

        implicit none

        ! Locals
        integer :: i, j
        character(len=200) :: msg, temp

        ! ===== HEADER =====
        call LogMessage("====== Import input file ======")
        call LogMessage("")
        write(msg, '("Reading file ", A)') trim(maniac_file)
        call LogMessage(msg)
        call LogMessage("")

        ! ===== GENERIC PARAMETERS =====
        call LogMessage("=== Generic parameters")
        write(msg, '("Number of blocks: ", I0)') input%nb_block
        call LogMessage(msg)
        write(msg, '("Number of steps: ", I0)') input%nb_step
        call LogMessage(msg)
        write(msg, '("Temperature (K): ", F10.2)') input%temp_K
        call LogMessage(msg)
        call LogMessage("")

        ! ===== ELECTROSTATIC =====
        call LogMessage("=== Electrostatic interactions")
        write(msg, '("Ewald tolerance: ", F15.8)') input%ewald_tolerance
        call LogMessage(msg)
        write(msg, '("Cutoff (Å): ", F10.2)') input%real_space_cutoff
        call LogMessage(msg)
        call LogMessage("")

        ! ===== MONTE CARLO =====
        call LogMessage("=== Monte carlo move")
        write(msg, '("Translation step (Å): ", F10.2)') input%translation_step
        call LogMessage(msg)
        write(msg, '("Rotation step angle (radian): ", F10.2)') input%rotation_step_angle
        call LogMessage(msg)

        write(msg, '("Translation proba: ", F10.2)') proba%translation
        call LogMessage(msg)
        write(msg, '("Rotation proba: ", F10.2)') proba%rotation
        call LogMessage(msg)
        write(msg, '("Big move proba: ", F10.2)') proba%big_move
        call LogMessage(msg)
        write(msg, '("Insertion deletion proba: ", F10.2)') proba%insertion_deletion
        call LogMessage(msg)
        write(msg, '("Swap proba: ", F10.2)') proba%swap
        call LogMessage(msg)
        call LogMessage("")

        ! ===== RESIDUES =====
        call LogMessage("=== Residue information")
        call LogMessage("")
        write(msg, '("Number of type of residue found: ", I0)') nb%type_residue
        call LogMessage(msg)
        call LogMessage("")

        do i = 1, nb%type_residue
            write(msg, '("  Residue ", A)') trim(res%names_1d(i))
            call LogMessage(msg)

            write(msg, '("  Is active: ", A)') merge("yes", "no ", input%is_active(i) == 1)
            call LogMessage(msg)

            if (input%is_active(i) == 1) then
                write(msg, '("  Fugacity (atm): ", F10.2)') input%fugacity(i)
                call LogMessage(msg)
            end if

            write(msg, '("  Number of atoms in residue: ", I0)') nb%atom_in_residue(i)
            call LogMessage(msg)

            write(msg, '("  Number of atom types in residue: ", I0)') nb%types_per_residue(i)
            call LogMessage(msg)

            ! Types
            msg = "  Types: "
            do j = 1, nb%types_per_residue(i)
                write(temp, '(I0)') res%types_2d(i, j)
                msg = trim(msg) // " " // trim(temp)
            end do
            call LogMessage(msg)

            ! Names
            msg = "  Names: "
            do j = 1, nb%types_per_residue(i)
                temp = res%names_2d(i, j)
                msg = trim(msg) // " " // trim(temp)
            end do
            call LogMessage(msg)

            call LogMessage("")
        end do

    end subroutine PrintInputSummary


end module output_utils
