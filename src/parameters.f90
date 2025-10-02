module parameters

    use, intrinsic :: iso_fortran_env, only: real64

    implicit none

    ! Integer parameters defining maximum sizes and counts
    integer, parameter :: NB_MAX_MOLECULE = 5000    ! Maximum number of molecules of a given type
    integer, parameter :: NB_MAX_BOND = 20000       ! Maximum number of bond per residue
    integer, parameter :: NB_MAX_ANGLE = 20000      ! Maximum number of angle per residue
    integer, parameter :: NB_MAX_DIHEDRAL = 20000   ! Maximum number of dihedral per residue
    integer, parameter :: NB_MAX_IMPROPER = 20000   ! Maximum number of improper per residue

    ! Parameters controlling the MC adjustement
    real(real64), parameter :: TARGET_ACCEPTANCE = 0.40      ! Target acceptance ratio for adaptive MC moves
    real(real64), parameter :: TOL_ACCEPTANCE    = 0.05d0    ! Tolerance window around target acceptance
    real(real64), parameter :: MIN_TRANSLATION_STEP = 1.0d-3 ! Minimum allowed translation displacement (Å)
    real(real64), parameter :: MAX_TRANSLATION_STEP = 3.0d0  ! Maximum allowed translation displacement (Å)
    real(real64), parameter :: MIN_ROTATION_ANGLE  = 1.0d-3  ! Minimum allowed rotation angle (rad)
    real(real64), parameter :: MAX_ROTATION_ANGLE  = 0.78d0  ! Maximum allowed rotation angle (~45° in rad)
    integer, parameter :: MIN_TRIALS_FOR_RECALIBRATION = 500 ! Minimum trials before recalibrating step sizes

    ! Parameters for unit conversion
    real(real64), parameter :: A3_TO_M3 = 1.0D-30    ! Conversion factor: 1 Å³ → 1.0E-30 m³
    real(real64), parameter :: ATM_TO_PA = 1.01325D5 ! Conversion factor: 1 atm → 1.01325 × 10^5 Pa

    ! File names
    character(len=*), parameter :: data_filename = 'topology.data' ! Topology data file

    ! Define constants for move types
    integer, parameter :: TYPE_CREATION = 1 ! Move type: molecule creation
    integer, parameter :: TYPE_DELETION = 2 ! Move type: molecule deletion
    integer, parameter :: TYPE_TRANSLATION = 3 ! Move type: molecular translation
    integer, parameter :: TYPE_ROTATION = 4 ! Move type: molecular rotation
    integer, parameter :: TYPE_SMALLMOVE = 0 ! Subtype: small translation move
    integer, parameter :: TYPE_BIGMOVE = 1 ! Subtype: large "jump" move

end module parameters
