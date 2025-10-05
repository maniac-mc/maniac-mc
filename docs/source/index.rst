MANIAC documentation
====================

MANIAC is a lightweight Monte Carlo simulation code written in Fortran,
designed for GCMC and adsorption studies. It reads basic LAMMPS-style topology
files and supports the following Monte Carlo moves:

- Translation  
- Rotation  
- Insertion  
- Deletion  
- Swap  

Why MANIAC?
-----------

The original MANIAC computer (for Mathematical Analyzer, Numerical Integrator, and
Computer) was built in the early 1950s at Los Alamos National Laboratory. It
was one of the first machines used to perform Monte Carlo simulations in
statistical physics and nuclear research.

LAMMPS compatibility
--------------------

MANIAC uses the same `.data` file format as LAMMPS for molecular
topology. MANIAC assumes that the real units system is used, and that
a pair style of the family `lj/cut/coul/long` is used.

.. toctree::
    :maxdepth: 2
    :caption: Contents:
    :hidden:

    input
    build
    unit
    credit
