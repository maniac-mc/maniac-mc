module simulation_state

    use constants
    use parameters
    use, intrinsic :: iso_fortran_env, only: real64

    implicit none

    character(len=200) :: output_path   ! Path for saving outputs
    character(len=200) :: maniac_file   ! Main input file
    character(len=200) :: data_file     ! Topology/data file
    character(len=200) :: inc_file      ! Parameters include file
    character(len=200) :: res_file      ! Optional reservoir file
    integer :: current_block            ! Current Monte Carlo block number
    integer :: current_step             ! Current Monte Carlo step within the block
    integer :: out_unit = 10            ! Default log file unit
    logical :: has_reservoir            ! Wether a reservoir was provided or

    type :: counter_type
        integer :: rotations = 0             ! Counter for rotational Monte Carlo moves
        integer :: translations = 0          ! Counter for translational Monte Carlo moves
        integer :: creations = 0             ! Counter for creation moves
        integer :: deletions = 0             ! Counter for deletion moves
        integer :: swap = 0                  ! Counter for swap moves
        integer :: trial_translations = 0    ! Counter for trial translation moves
        integer :: trial_rotations = 0       ! Counter for trial rotation moves
        integer :: trial_creations = 0       ! Counter for trial rotation moves
        integer :: trial_deletions = 0       ! Counter for trial deletion moves
        integer :: trial_swap = 0            ! Counter for swap moves
    end type counter_type
    type(counter_type) :: counter

    ! Monte carlo move probability
    type :: proba_type
        real(real64) :: translation             ! Probability of attempting a translation move
        real(real64) :: rotation                ! Probability of attempting a rotation move
        real(real64) :: insertion_deletion      ! Probability of attempting an insertion/deletion
        real(real64) :: swap                    ! Probability of attempting a swap move
    end type proba_type
    type(proba_type) :: proba

    ! Parameters provided in the input file
    type :: input_type
        real(real64), dimension(:), allocatable :: fugacity ! Fugacity array of the GCMC reservoir
        integer, dimension(:), allocatable :: is_active ! Activity flags or counts for each molecule type
        real(real64) :: temp_K                  ! Temperature in Kelvin
        real(real64) :: translation_step        ! Maximum displacement for MC moves
        real(real64) :: rotation_step_angle     ! Maximum rotation for MC moves
        real(real64) :: ewald_tolerance         ! Numerical accuracy for Ewald summation,
        real(real64) :: real_space_cutoff       ! Cutoff radius - maximum interaction distance in real space
        integer :: seed                         ! Initial seed for the random number generator
        integer :: nb_block                     ! Total desired Monte Carlo block number
        integer :: nb_step                      ! Total desired Monte Carlo step
        logical :: recalibrate_moves            ! Enable automatic recalibration of move steps (true/false)
    end type input_type
    type(input_type) :: input

    ! Energy terms
    type :: energy_type
        real(real64) :: self_interaction          ! Site-site short-range energy
        real(real64) :: coulomb                   ! Charged-electrostatic interaction energy
        real(real64) :: non_coulomb               ! Neutral-charged interaction energy
        real(real64) :: ewald_self                ! Ewald self-interaction energy
        real(real64) :: ke_reciprocal             ! Reciprocal-space (k-space) energy
        real(real64) :: intra_coulomb             ! Intramolecular Coulomb energy (alternative)
        real(real64) :: total_coulomb             ! Total Coulomb energy
        real(real64) :: total_non_coulomb         ! Total non-Coulomb energy
        real(real64) :: recip_coulomb             ! Reciprocal-space Coulomb contribution
        real(real64) :: total                     ! Total system energy
    end type energy_type
    type(energy_type) :: energy

    ! Energy status for Monte Carlo move
    type :: energy_state
        real(real64) :: non_coulomb             ! Old and new non-coulombic energy of molecule
        real(real64) :: coulomb                 ! Old and new coulombic energy of molecule
        real(real64) :: recip_coulomb           ! Old and new reciprocal-space energy
        real(real64) :: ewald_self
        real(real64) :: intra_coulomb
        real(real64) :: total                   ! Total old and new energies
    end type energy_state
    type(energy_state) :: old, new, old_bis, new_bis, intermediate

    ! Simulation box definition
    type :: type_box
        ! Generic box geometry parameters
        real(real64) :: reciprocal(3,3)     ! Reciprocal box matrix
        real(real64) :: matrix(3,3)         ! Simulation box matrix
        real(real64) :: bounds(3,2)         ! Box dimensions (lo, hi)
        real(real64) :: metrics(9)          ! Misc. box properties (e.g., lengths, cosines of angles)
        real(real64) :: tilt(3)             ! Tilt factors (xy, xz, yz)
        real(real64) :: determinant         ! Volume scaling factor of a linear transformation
        real(real64) :: volume              ! Box volume
        ! About box content
        logical :: is_triclinic             ! Indicate if box is triclinic
        integer :: type = 0                 ! 0 = unset, 1 = cubic, 2 = orthorhombic, 3 = triclinic
        integer :: num_atoms                ! Number of atoms in the box
        integer :: num_atomtypes            ! Number of atom types
        integer :: num_bonds                ! Number of bonds
        integer :: num_bondtypes            ! Number of bond types
        integer :: num_angles               ! Number of angles
        integer :: num_angletypes           ! Number of angle types
        integer :: num_dihedrals            ! Number of dihedrals
        integer :: num_dihedraltypes        ! Number of dihedral types
        integer :: num_impropers            ! Number of impropers
        integer :: num_impropertypes        ! Number of improper types
        integer, allocatable :: num_residues(:) ! Number of residue of each type
        ! Atom information
        real(real64), dimension(:), allocatable :: site_masses_vector    ! Mass vector for all atom types in initial inputs
        real(real64), dimension(:, :), allocatable :: atom_charges       ! Partial charges on sites (1,:,:)=system, (2,:,:)=reservoir
        real(real64), dimension(:, :), allocatable :: atom_masses        ! Masses of atoms (1,:,:)=system, (2,:,:)=reservoir
        character(len=10), dimension(:, :), allocatable :: atom_names    ! Atom names for each residue (1,:,:)=system, (2,:,:)=reservoir
        integer, dimension(:, :), allocatable :: atom_types              ! Atom types for each residue (1,:,:)=system, (2,:,:)=reservoir
        integer, dimension(:, :), allocatable :: atom_ids                ! Atom ids for each residue (1,:,:)=system, (2,:,:)=reservoir
        real(real64), dimension(:, :, :), allocatable :: mol_com         ! X Y Z coordinate of molecule centers or atoms (1,:,:,:)=system, (2,:,:,:)=reservoir
        real(real64), dimension(:, :, :, :), allocatable :: site_offset  ! Local site X Y Z displacements from molecule center (1,:,:,:,:)=system, (2,:,:,:,:)=reservoir
    end type type_box
    type(type_box) :: primary, reservoir

    ! Simulation box definition
    type :: type_number
        integer :: type_residue             ! Total number of residues
        integer :: max_atom_in_residue      ! Max number of atoms in the largest residue
        integer :: max_type_per_residue     ! Max number of type per residue
        integer, dimension(:), allocatable :: atom_in_residue ! Number of atoms in the residue
        integer, dimension(:), allocatable :: types_per_residue ! Number of atom types in the residue
        integer, dimension(:), allocatable :: bonds_per_residue ! Number of bonds in the residue
        integer, dimension(:), allocatable :: angles_per_residue ! Number of angles in the residue
        integer, dimension(:), allocatable :: dihedrals_per_residue ! Number of dihedrals in the residue
        integer, dimension(:), allocatable :: impropers_per_residue ! Number of impropers in the residue
        integer, dimension(:, :), allocatable :: types_pattern ! Type pattern in residue (eg, for TIP4P water 1 2 3 3)
    end type type_number
    type(type_number) :: nb

    ! Residues information
    type :: type_residue
        real(real64), dimension(:), allocatable :: masses_1d           ! Array of atoms masses
        character(len=10), dimension(:), allocatable :: names_1d       ! Array of residue names
        character(len=10), dimension(:, :), allocatable :: names_2d    ! Site names for each residue
        integer, dimension(:, :), allocatable :: types_2d              ! Site types for each residue
        integer, dimension(:, :, :), allocatable :: bond_type_2d       ! Site bonds for each residue
        integer, dimension(:, :, :), allocatable :: angle_type_2d      ! Site angles for each residue
        integer, dimension(:, :, :), allocatable :: dihedral_type_2d   ! Site dihedrals for each residue
        integer, dimension(:, :, :), allocatable :: improper_type_2d   ! Site impropers for each residue
    end type type_residue
    type(type_residue) :: res

    ! Interaction arrays
    type :: type_coeff
        real(real64), dimension(:, :, :, :), allocatable :: sigma       ! Lennard-Jones coefficients
        real(real64), dimension(:, :, :, :), allocatable :: epsilon     ! Lennard-Jones oefficients
    end type type_coeff
    type(type_coeff) :: coeff

    ! Type to store precomputed reciprocal vectors
    type :: kvector_type
        integer :: kx
        integer :: ky
        integer :: kz
        real(real64) :: k_squared       ! Normalized k^2 (for validity)
        real(real64) :: k_squared_mag   ! Cartesian squared magnitude (for W(k))
    end type kvector_type

    ! For Ewald calculation
    type :: type_ewald
        integer :: num_kvectors                 ! Number of precomputed vectors
        integer :: kmax(3)                      ! Maximum index for reciprocal lattice vector x y z component
        real(real64) :: alpha                   ! Ewald summation alpha parameter (screening parameter)
        real(real64) :: fourier_precision       ! Estimated precision required for reciprocal (Fourier) space summation
        real(real64) :: screening_factor        ! Intermediate tolerance factor for real-space screening width calculation
        real(real64), dimension(:), allocatable :: recip_constants ! Constants for reciprocal space summations
        complex(real64), dimension(:, :, :, :), allocatable :: phase_factor_x ! Complex exponentials for reciprocal space (x-component)
        complex(real64), dimension(:, :, :, :), allocatable :: phase_factor_y ! Complex exponentials for reciprocal space (y-component)
        complex(real64), dimension(:, :, :, :), allocatable :: phase_factor_z ! Complex exponentials for reciprocal space (z-component)
        complex(real64), dimension(:, :), allocatable :: phase_factor_x_old ! Old Fourier terms along X
        complex(real64), dimension(:, :), allocatable :: phase_factor_y_old ! Old Fourier terms along Y
        complex(real64), dimension(:, :), allocatable :: phase_factor_z_old ! Old Fourier terms along Z
        complex(real64), dimension(:), allocatable :: recip_amplitude ! Fourier coefficients of charge density or potential
        complex(real64), dimension(:), allocatable :: recip_amplitude_old ! Old fourier coefficients of charge density or potential
        real(real64), dimension(:), allocatable :: form_factor ! Factor to account for symmetry (k vs -k)
        type(kvector_type), allocatable :: kvectors(:) ! Precomputed reciprocal vectors
    end type type_ewald
    type(type_ewald) :: ewald

    ! For tabulating potential
    type :: tabulated
        ! Table properties
        real(real64), allocatable :: x(:)       ! Grid points (r, k^2, etc.)
        real(real64), allocatable :: f(:)       ! Function values
        real(real64) :: dx                      ! Grid spacing
        integer :: n                            ! Number of points
        logical :: initialized = .false.        ! Flag to indicate table is ready
    end type tabulated
    type(tabulated) :: erfc_r_table             ! For precomputed erfc(r) / r
    type(tabulated) :: r6_table                 ! For precomputed r**6
    type(tabulated) :: r12_table                ! For precomputed r**12

end module simulation_state
