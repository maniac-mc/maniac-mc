module energy_utils

    use simulation_state
    use geometry_utils
    use constants
    ! use ewald_utils

    use, intrinsic :: iso_fortran_env, only: real64
    use, intrinsic :: ieee_arithmetic ! To remove eventually

    implicit none

contains

    subroutine ComputeSystemEnergy(box)

        implicit none

        ! Input arguments
        type(type_box), intent(inout) :: box

        ! Compute each energy components
        call ComputePairwiseEnergy(box)
        call ComputeEwaldSelf(box)
        call ComputeEwaldRecip(box)

    end subroutine ComputeSystemEnergy

    subroutine ComputePairwiseEnergy(box)

        implicit none
        
        ! Input arguments
        type(type_box), intent(inout) :: box

        ! Variables for the residue of interest
        integer :: residue_type_1
        integer :: molecule_index_1
        real(real64) :: e_non_coulomb
        real(real64) :: e_coulomb

        ! First, initialize energies to zero
        energy%non_coulomb = zero
        energy%coulomb = zero

        ! Loop over all residue types
        do residue_type_1 = 1, nb%type_residue
            ! Loop over all molecule of type "residue_type_1"
            do molecule_index_1 = 1, box%num_residues(residue_type_1)

                ! Compute the energy for residue_1, molecule_1
                call SingleMolPairwiseEnergy(box, residue_type_1, &
                    molecule_index_1, e_non_coulomb, e_coulomb)

                ! Add residue_1, molecule_1 energy to the total pairwise energy
                energy%non_coulomb = energy%non_coulomb + e_non_coulomb
                energy%coulomb = energy%coulomb + e_coulomb
            end do
        end do

        write(*,*) energy%non_coulomb, energy%coulomb

    end subroutine ComputePairwiseEnergy

    !------------------------------------------------------------------------------
    ! subroutine SingleMolPairwiseEnergy
    ! Calculates the non-Coulombian and Coulomb (direct space)
    !------------------------------------------------------------------------------
    subroutine SingleMolPairwiseEnergy(box, residue_type_1, molecule_index_1, e_non_coulomb, e_coulomb)

        implicit none

        ! Input arguments
        type(type_box), intent(inout) :: box
        integer, intent(in) :: residue_type_1        ! Residue type to be moved
        integer, intent(in) :: molecule_index_1      ! Molecule ID
        real(real64), intent(out) :: e_non_coulomb
        real(real64), intent(out) :: e_coulomb

        ! Local variables
        integer :: atom_index_1, atom_index_2
        integer :: molecule_index_2
        integer :: residue_type_2
        real(real64) :: distance
        real(real64) :: r6, r12                     ! r^n for LJ potential calculations
        real(real64) :: sigma, epsilon              ! Epsilon and sigma LJ potential calculations
        real(real64) :: charge_1, charge_2          ! Charge for Coulomb interactions

        e_non_coulomb = zero
        e_coulomb = zero

        ! Loop over sites in molecule residue_type
        do atom_index_1 = 1, nb%atom_in_residue(residue_type_1) 

            ! Loop over all molecule types 2
            do residue_type_2 = 1, nb%type_residue

                ! Loop over all molecule index 2
                do molecule_index_2 = 1, box%num_residues(residue_type_2)

                    ! Remove intra molecular contribution
                    if ((molecule_index_1 == molecule_index_2) .and. &
                        (residue_type_1 == residue_type_2)) cycle

                    ! Enforce ordering to avoid double-counting
                    if ((residue_type_2 < residue_type_1) .or. &
                        ((residue_type_2 == residue_type_1) .and. (molecule_index_2 <= molecule_index_1))) cycle

                    ! Loop over all side of the selected molecule 2
                    do atom_index_2 = 1, nb%atom_in_residue(residue_type_2)

                        ! Read pair parameters
                        sigma = coeff%sigma(residue_type_1, residue_type_2, atom_index_1, atom_index_2)                            
                        epsilon = coeff%epsilon(residue_type_1, residue_type_2, atom_index_1, atom_index_2)
                        charge_1 = primary%atom_charges(residue_type_1, atom_index_1)
                        charge_2 = primary%atom_charges(residue_type_2, atom_index_2)

                        distance = ComputeDistance(box, residue_type_1, molecule_index_1, atom_index_1, &
                                   residue_type_2, molecule_index_2, atom_index_2)

                        ! Update non-Coulomb energy
                        e_non_coulomb = e_non_coulomb + LennardJonesEnergy(distance, sigma, epsilon)

                        ! Update Coulomb energy
                        e_coulomb = e_coulomb + CoulombEnergy(distance, charge_1, charge_2)

                    end do
                end do
            end do
        end do

        ! Re-scale energy
        e_coulomb = e_coulomb * EPS0_INV_eVA / KB_eVK

        return

    end subroutine SingleMolPairwiseEnergy

    !------------------------------------------------------------------------------
    ! Function to compute Lennard-Jones interaction energy
    !------------------------------------------------------------------------------
    pure function LennardJonesEnergy(distance, sigma, epsilon) result(energy)

        implicit none

        ! Input variables
        real(real64), intent(in) :: distance   ! Distance between the two atoms
        real(real64), intent(in) :: sigma      ! Lennard-Jones sigma parameter (size, Å)
        real(real64), intent(in) :: epsilon    ! Lennard-Jones epsilon parameter (well depth, eV or kJ/mol depending on units)

        ! Local variables
        real(real64) :: r6      ! (σ / r)^6 term of the LJ potential
        real(real64) :: r12     ! (σ / r)^12 term of the LJ potential
        real(real64) :: energy  ! Computed Lennard-Jones energy contribution

        if (distance >= input%real_space_cutoff) then
            energy = zero
        else
            r6 = (sigma / distance)**6
            r12 = r6 * r6
            energy = four * epsilon * (r12 - r6)
        end if

    end function LennardJonesEnergy

    !------------------------------------------------------------------------------
    ! Function to compute Coulomb interaction energy (Ewald direct-space term)
    !------------------------------------------------------------------------------
    pure function CoulombEnergy(distance, charge1, charge2) result(energy)

        implicit none

        ! Input variables
        real(real64), intent(in) :: distance   ! Distance between the two atoms
        real(real64), intent(in) :: charge1    ! Atomic partial charge of atom 1 (in e)
        real(real64), intent(in) :: charge2    ! Atomic partial charge of atom 2 (in e)

        ! Local variable
        real(real64) :: energy   ! Computed Coulomb energy contribution (unscaled)

        if ((abs(charge1) < error) .or. (abs(charge2) < error)) then
            energy = zero
        else
            ! Direct-space Coulomb potential with Ewald damping:
            ! V(r) = (q1*q2) * erfc(alpha * r) / r
            energy = charge1 * charge2 * (erfc(ewald%alpha * distance)) / distance
        end if

    end function CoulombEnergy

    subroutine ComputeEwaldSelf(box)

        implicit none

        ! Input arguments
        type(type_box), intent(inout) :: box

    end subroutine ComputeEwaldSelf

    subroutine ComputeEwaldRecip(box)

        implicit none

        ! Input arguments
        type(type_box), intent(inout) :: box

    end subroutine ComputeEwaldRecip

end module energy_utils
