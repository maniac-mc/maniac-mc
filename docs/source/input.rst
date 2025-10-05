Input Files
===========

This page describes the structure and parameters of the input file
used by the **MANIAC** Monte Carlo simulation program. Each line in the file
specifies a parameter and its value. Lines starting with ``#`` are comments
and are ignored by the program.

Minimal Example
---------------

Below is a minimal example for a system with two residue types: water molecules
(the active species) and ZIF-8 (the rigid host framework).

.. code-block:: text

    # Input for maniac program

    # Generic parameters
    nb_block 200                     # Number of Monte Carlo blocks
    nb_step 1000                     # Number of Monte Carlo steps per block
    temperature 300                  # Simulation temperature (K)

    # Electrostatic interactions
    ewald_tolerance 1e-5             # Ewald convergence tolerance (smaller = more accurate, slower)
    real_space_cutoff 17             # Real-space cutoff distance for electrostatics (Å)

    # Monte Carlo moves
    translation_step 1               # Maximum translational displacement (Å)
    rotation_step_angle 0.685        # Maximum rotation angle (radians)
    recalibrate_moves true           # Automatically tune step sizes to keep acceptance ~optimal

    translation_proba 0.4            # Probability of attempting a translation move
    rotation_proba 0.4               # Probability of attempting a rotation move
    insertion_deletion_proba 0.2     # Probability of attempting insertion or deletion of molecules

    # Define residue composition and activity
    begin_residue
    name wat                       # Residue name: water
    state actif                    # "actif" = subject to MC moves
    fugacity 10000                 # Reservoir fugacity (driving force for adsorption)
    types 1 2 3                    # Atom type IDs for this residue
    names Ow Mw Hw                 # Atom names: oxygen, virtual site, hydrogens
    nb-atoms 4                     # Number of atoms in this residue
    end_residue

    begin_residue
    name zif                       # Residue name: ZIF-8 framework
    state inactif                  # "inactif" = rigid, no MC moves applied
    types 4 5 6 7 8 9 10 11        # Atom type IDs for the framework
    names Zn C1 C2 H1 C3 H2 H3 N   # Atom names in ZIF-8 (Zn, carbons, hydrogens, nitrogen)
    nb-atoms 2208                  # Total number of atoms in the framework
    end_residue

This example corresponds to water adsorption in a rigid ZIF-8:  

- The first residue block defines a TIP4P water molecule (``wat``).  
- It is marked as **actif**, meaning Monte Carlo moves (translation,
  rotation, insertion, deletion) will be attempted.  
- A high fugacity (``10000 atm``) is used to mimics a rapid adsorption.  
- Atom types (1 2 3) and names (Ow Mw Hw) correspond to the oxygen (Ow),
  a virtual site (Mw) and hydrogen atoms (Hw),
- ``nb-atoms 4`` specifies the number of atoms in a molecule.  

- The second residue block defines the **ZIF-8 framework** (``zif``).  
- It is marked as **inactif**, so it remains rigid and is not subject to Monte
  Carlo moves.  
- Atom types and names (Zn C1 C2 H1 C3 H2 H3 N) represent the crystallographic
  composition of ZIF-8.  
- ``nb-atoms 2208`` specifies the number of atoms in the framework.

The probabilities for translation, rotation, big moves, and insertion/deletion
(`translation_proba`, `rotation_proba`, `big_move_proba`,
`insertion_deletion_proba`) do not need to sum to 1. When all probabilities are
provided, **MANIAC** will automatically rescale them so that their total equals 1.
This ensures that at each Monte Carlo step, one move is chosen according to the
relative weights of the provided probabilities.

Simulation Control
------------------

Number of blocks
################

**nb_block** (Integer, required)  
Number of Monte Carlo **blocks** to run.  
A block is a group of MC steps used for statistical averaging.

Number of steps per block
#########################

**nb_step** (Integer, required)  
Number of Monte Carlo **steps** per block.

Simulation temperature
######################

**temperature** (Float, required, units: K)  
Target simulation temperature.

Random seed
###########

**seed** (Integer, optional, default: random)  
Initial seed for the random number generator. Using the same seed
will reproduce identical simulation trajectories.  
If omitted, the seed is chosen randomly.

Electrostatics
--------------

Ewald summation tolerance
#########################

**ewald_tolerance** (Float, optional, default: 1e-6)  
Convergence tolerance for Ewald summation of electrostatic interactions.  
Smaller values increase accuracy but increase cost.

Real-space cutoff distance
##########################

**real_space_cutoff** (Float, optional, default: 10 Å)  
Cutoff distance in real space for electrostatics.

Monte Carlo Moves
-----------------

Maximum translation step
########################

**translation_step** (Float, required, units: Å)  
Maximum translational displacement allowed in a move.

Maximum rotation angle
######################

**rotation_step_angle** (Float, required, units: radians)  
Maximum rotation angle allowed in a move.

Recalibrate moves automatically
###############################

**recalibrate_moves** (Boolean, optional, default: false)  
If true, automatically adjust step sizes to maintain acceptance ratios.

Move probabilities
##################

- **translation_proba** (Float, optional, default: 0.0)  
  Probability of attempting a translation move.  

- **rotation_proba** (Float, optional, default: 0.0)  
  Probability of attempting a rotation move.  

- **big_move_proba** (Float, optional, default: 0.0)  
  Probability of attempting a "big" displacement move.  

- **insertion_deletion_proba** (Float, optional, default: 0.0)  
  Probability of attempting a molecule insertion or deletion.

Residue Definition
------------------

A residue block defines the composition and properties of each molecule type.  
Each block starts with ``begin_residue`` and ends with ``end_residue``.

Residue parameters:

- **name** (String, required): Residue name (e.g., ``CO2``, ``H2O``, ``MOF``).  
- **state** (String, required): ``actif`` (active, subject to MC moves) or ``inactif`` (rigid).  
- **fugacity** (Float, optional, units: bar): Chemical potential of the residue.  
- **types** (List of Integers, required): Atom type IDs for atoms in the residue.  
- **names** (List of Strings, required): Atom names corresponding to atom types.  
- **nb-atoms** (Integer, required): Number of atoms in one residue.

Notes
-----
- Parameter order is not strict.  
- Each residue must be fully enclosed between ``begin_residue`` and ``end_residue``.  
