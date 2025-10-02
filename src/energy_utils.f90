module energy_utils

    use simulation_state
    use geometry_utils
    ! use constants
    ! use ewald_utils

    use, intrinsic :: iso_fortran_env, only: real64
    use, intrinsic :: ieee_arithmetic ! To remove eventually

    implicit none

contains

    subroutine ComputeSystemEnergy()

        implicit none

        ! Compute each energy components
        call ComputePairwiseEnergy()
        call ComputeEwaldSelf()
        call ComputeEwaldRecip()

    end subroutine ComputeSystemEnergy

    subroutine ComputePairwiseEnergy()

        implicit none
        
        ! Variables for the residue of interest
        integer :: residue_type_1
        integer :: molecule_index_1
        real(real64) :: e_non_coulomb
        real(real64) :: e_coulomb

        energy%non_coulomb = 0.0_real64 
        energy%coulomb = 0.0_real64

        ! Loop over all residue types
        do residue_type_1 = 1, nb%type_residue
            ! Loop over all molecule of type "residue_type_1"
            do molecule_index_1 = 1, primary%num_residues(residue_type_1)
                ! Compute the energy
                call ComputePairInteractionEnergy_singlemol(primary, residue_type_1, molecule_index_1, &
                    e_non_coulomb, e_coulomb)
                
                ! Remove double couting
                e_coulomb = e_coulomb / 2
                e_non_coulomb = e_non_coulomb / 2

                energy%non_coulomb = energy%non_coulomb + e_non_coulomb
                energy%coulomb = energy%coulomb + e_coulomb
            end do
        end do

        write (*,*) energy%non_coulomb, energy%coulomb

    end subroutine ComputePairwiseEnergy

    !------------------------------------------------------------------------------
    ! subroutine ComputePairInteractionEnergy_singlemol
    ! Calculates the non-Coulombian and Coulomb (direct space)
    !------------------------------------------------------------------------------
    subroutine ComputePairInteractionEnergy_singlemol(box, residue_type_1, molecule_index_1, e_non_coulomb, e_coulomb)

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
        ! real(real64), dimension(3) :: delta
        real(real64) :: distance
        real(real64) :: r6, r12 ! r^n for LJ potential calculations
        real(real64) :: sigma, epsilon ! Epsilon and sigma LJ potential calculations
        real(real64) :: charge_1, charge_2 ! Charge for Coulomb interactions

        e_non_coulomb = 0.0_real64
        e_coulomb = 0.0_real64

        ! Loop over sites in molecule residue_type
        do atom_index_1 = 1, nb%atom_in_residue(residue_type_1) 

            ! Loop over all molecule types 2
            do residue_type_2 = 1, nb%type_residue

                ! Loop over all molecule index 2
                do molecule_index_2 = 1, box%num_residues(residue_type_2)

                    ! Remove intra molecular contribution
                    if ((molecule_index_1 == molecule_index_2) .and. &
                        (residue_type_1 == residue_type_2)) cycle

                    ! Loop over all side of the selected molecule 2
                    do atom_index_2 = 1, nb%atom_in_residue(residue_type_2)

                        distance = ComputeDistance(box, residue_type_1, molecule_index_1, atom_index_1, &
                                   residue_type_2, molecule_index_2, atom_index_2)

                        if (distance < input%real_space_cutoff) then
                            ! LJ potential
                            sigma = coeff%sigma(residue_type_1, residue_type_2, atom_index_1, atom_index_2)                            
                            epsilon = coeff%epsilon(residue_type_1, residue_type_2, atom_index_1, atom_index_2)
                            r6 = (sigma / distance)**6
                            r12 = r6 * r6
                            e_non_coulomb = e_non_coulomb + 4.0_real64 * epsilon * (r12 - r6)
                        end if

                        ! Use Coulomb potential
                        charge_1 = primary%atom_charges(residue_type_1, atom_index_1)
                        charge_2 = primary%atom_charges(residue_type_2, atom_index_2)

                        if ((abs(charge_1) < 1.0D-10) .or. (abs(charge_2) < 1.0D-10)) cycle

                        e_coulomb = e_coulomb + charge_1 * charge_2 * (erfc(ewald%alpha * distance)) / distance

                    end do
                end do
            end do
        end do

        ! Re-scale energy
        e_coulomb = e_coulomb * EPS0_INV_eVA / KB_eVK

        return

    end subroutine ComputePairInteractionEnergy_singlemol

    subroutine ComputeEwaldSelf()

        implicit none

    end subroutine ComputeEwaldSelf

    subroutine ComputeEwaldRecip()

        implicit none

    end subroutine ComputeEwaldRecip

end module energy_utils
