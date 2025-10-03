module parameters_parser

    use parameters
    use simulation_state
    use output_utils
    use check_utils
    use readers_utils
    use, intrinsic :: iso_fortran_env, only: real64

    implicit none

contains

    !---------------------------------------------------------------
    ! Subroutine: ReadParameters
    ! Reads Lennard-Jones parameters (epsilon, sigma) from input file,
    ! assigns them to atom type pairs in the system, and applies
    ! mixing rules (Lorentz-Berthelot) if needed.
    !---------------------------------------------------------------
    subroutine ReadParameters()

        implicit none

        ! Local variables
        integer :: ios                                      ! I/O status code for file operations
        integer :: i, j, k, l                               ! Loop indices for iterating over residues and atoms
        integer :: pos_val                                  ! Position in the input line after the keyword for parsing
        integer :: val_int1, val_int2                       ! Temporary atom type indices read from pair_coeff line
        integer :: type_i, type_j                           ! Atom type indices for matching pair_coeff to residue/atom pairs
        real(real64) :: epsilon, sigma                      ! Lennard-Jones parameters: epsilon (kcal/mol, converted to Kelvin) and sigma (Å)
        integer :: INFILE = 3000                            ! File unit number for opening the parameter input file
        character(len=100) :: line                          ! Buffer to store each line read from the input file
        character(len=100) :: keyword                       ! Keyword extracted from the input line (e.g., "pair_coeff")
        character(len=100) :: rest_line                     ! Remaining part of the input line after the keyword
        logical, allocatable :: printed(:,:)                ! 2D array to track which atom type pairs have been logged
        integer, allocatable :: pair1(:), pair2(:)          ! Arrays to store first and second atom type indices of unique pairs
        real(real64), allocatable :: epsilons(:), sigmas(:) ! Arrays to store epsilon (kcal/mol) and sigma (Å) for unique pairs
        integer :: max_atom_type                            ! Maximum atom type index from atom_types array
        integer :: max_pairs                                ! Maximum possible number of unique atom type pairs
        integer :: n_pairs                                  ! Counter for the number of unique atom type pairs found

        ! Determine the maximum atom type index for array allocation
        if (has_reservoir) then
            max_atom_type = max(maxval(primary%atom_types(:,:)), maxval(reservoir%atom_types(:,:)) )
        else
            max_atom_type = maxval(primary%atom_types(:,:))
        end if
        max_pairs = max_atom_type * max_atom_type           ! Calculate maximum possible pairs for array sizing

        allocate(printed(max_atom_type, max_atom_type))     ! Allocate array to track logged atom type pairs
        allocate(pair1(max_pairs), pair2(max_pairs))        ! Allocate arrays for storing atom type indices of pairs
        allocate(epsilons(max_pairs), sigmas(max_pairs))    ! Allocate arrays for storing epsilon and sigma values
        coeff%sigma = 0.0_real64
        coeff%epsilon = 0.0_real64

        printed = .false.                                   ! Initialize printed array to false (no pairs logged yet)

        open(UNIT=INFILE, FILE=inc_file, STATUS='OLD', ACTION='read')

        do
            read(INFILE, '(A)', IOSTAT=ios) line
            if (ios /= 0) EXIT   ! End of file
            if (line(1:1) == '#' .OR. LEN_TRIM(line) == 0) cycle ! skip comments/blank lines

            read(line, *) keyword

            pos_val = INDEX(line, keyword) + LEN_TRIM(keyword) + 1 ! move past keyword + space
            if (pos_val > LEN(line)) then
                rest_line = ''
            else
                rest_line = ADJUSTL(line(pos_val:))
            end if

            select case (TRIM(ADJUSTL(keyword)))
                case ("pair_coeff")
                read(rest_line, *, IOSTAT=ios) val_int1, val_int2, epsilon, sigma
                if (ios /= 0) then
                    call AbortRun("Failed to read pair_coeff value", 1)
                end if

                ! Convert epsilon from kcal/mol to Kelvin
                epsilon = epsilon / KB_kcalmol
            end select

            do i = 1, nb%type_residue
                do k = 1, nb%atom_in_residue(i)
                    do j = 1, nb%type_residue
                        do l = 1, nb%atom_in_residue(j)
                            type_i = primary%atom_types(i, k)
                            type_j = primary%atom_types(j, l)
                            if ((type_i == val_int1) .AND. (type_j == val_int2)) then
                                ! Record parameters
                                coeff%sigma(i, j, k, l) = sigma
                                coeff%epsilon(i, j, k, l) = epsilon

                                ! Also record the symmetric
                                coeff%sigma(j, i, l, k) = sigma
                                coeff%epsilon(j, i, l, k) = epsilon
                            end if
                        end do
                    end do
                end do
            end do
        end do

        close(INFILE)

        ! If cross parameters are missing, enforce Lorenz-Berthelot rule
        call ApplyLorentzBerthelot()

        ! Print parameters in log
        call LogParameters(inc_file, n_pairs, pair1, pair2, epsilons, sigmas)

    end subroutine ReadParameters

    subroutine ApplyLorentzBerthelot()

        implicit none

        integer :: i, j, k, l               ! Loop indices for iterating over residues and atoms
        integer :: type_i, type_j           ! Atom type indices for the current pair of atoms
        real(real64) :: sigma, epsilon      ! Lennard-Jones parameters: sigma (Å) and epsilon (Kelvin) for the pair interaction
        character(len=200) :: formatted_msg ! String to hold formatted log messages for output
        logical, allocatable :: warned(:,:) ! 2D array to track if a warning has been logged for a pair of atom types
        integer :: max_atom_type            ! Maximum atom type index from atom_types_2d array

        ! Determine the maximum atom type index for array allocation
        if (has_reservoir) then
            max_atom_type = max(maxval(primary%atom_types(:,:)), maxval(reservoir%atom_types(:,:)) )
        else
            max_atom_type = maxval(primary%atom_types(:,:))
        end if
        allocate(warned(max_atom_type, max_atom_type))  ! Allocate warned array to track logged warnings for atom type pairs
        warned = .false.

        do i = 1, nb%type_residue
            do k = 1, nb%atom_in_residue(i)
                do j = 1, nb%type_residue
                    do l = 1, nb%atom_in_residue(j)

                        ! Detect cross coefficients that are zeros
                        if ((abs(coeff%epsilon(i,j,k,l)) < 1.0e-6_real64) .and. abs(coeff%sigma(i,j,k,l)) < 1.0e-6_real64) then

                            ! Detect direct coefficients
                            type_i = primary%atom_types(i, k)
                            type_j = primary%atom_types(j, l)

                            ! Apply Lorentz-Berthelot rule
                            sigma = (coeff%sigma(i, i, k, k) + coeff%sigma(j, j, l, l)) / 2
                            epsilon = sqrt(coeff%epsilon(i, i, k, k) * coeff%epsilon(j, j, l, l))

                            ! If direct coefficients are both positive, then record cross coefficients and print warning
                            if ((sigma > 1.0e-6_real64) .and. (epsilon > 1.0e-6_real64)) then
                                if (.not. warned(type_i, type_j) .and. .not. warned(type_j, type_i)) then

                                    ! Print the "rule enforcement" warning only once globally
                                    if (.not. any(warned)) then
                                        call InfoMessage("Enforcing the Lorentz-Berthelot rule")
                                        call LogMessage("typei typej epsilon sigma")
                                    end if

                                    if (type_i < type_j) then
                                        write(formatted_msg,'("",I3," -",I3," : ",F8.4," Å ",F8.4," kcal/mol")') &
                                            type_i, type_j, sigma, epsilon * KB_kcalmol
                                    else
                                        write(formatted_msg,'("",I3," -",I3," : ",F8.4," Å ",F8.4," kcal/mol")') &
                                            type_j, type_i, sigma, epsilon * KB_kcalmol
                                    end if
                                    call LogMessage(formatted_msg)

                                    warned(type_i, type_j) = .true.
                                    warned(type_j, type_i) = .true.
                                end if
                                coeff%sigma(i,j,k,l) = sigma
                                coeff%epsilon(i,j,k,l) = epsilon
                            end if
                        end if
                    end do
                end do
            end do
        end do
    end subroutine ApplyLorentzBerthelot

end module parameters_parser
