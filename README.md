# MANIAC

[![Build](https://github.com/maniac-mc/maniac-mc/actions/workflows/tests.yml/badge.svg)](https://github.com/maniac-mc/maniac-mc/actions/workflows/tests.yml)
![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)
![Last Commit](https://img.shields.io/github/last-commit/maniac-mc/maniac-mc)
![GitHub release](https://img.shields.io/github/v/release/maniac-mc/maniac-mc)

<img
    src="https://raw.githubusercontent.com/maniac-mc/mc-visuals/refs/heads/main/gallery/ZIF8-H2O/system.png"
    width="30%" align="right"/>
</a>

MANIAC is a lightweight Monte Carlo simulation code written in Fortran,
designed for GCMC and adsorption studies. It reads basic LAMMPS-style topology
files and supports the following Monte Carlo moves:

- Translation  
- Rotation  
- Insertion  
- Deletion  
- Swap  

## Why MANIAC?

The original MANIAC computer (for Mathematical Analyzer, Numerical Integrator, and
Computer) was built in the early 1950s at Los Alamos National Laboratory. It
was one of the first machines used to perform Monte Carlo simulations in
statistical physics and nuclear research.

## LAMMPS compatibility

MANIAC uses the same `.data` file format as LAMMPS for molecular
topology. MANIAC assumes that the real units system is used, and that
a pair style of the family `lj/cut/coul/long` is used.

## Build the documentation

To build the documentation locally, navigate to the `docs` directory and run:

```bash
cd docs
doxygen doxyfile      # Generates XML and HTML documentation
make html             # Builds the HTML documentation
```

This will create the complete documentation in the `build/html` folder, that can
be open using:

```bash
firefox build/html/index.html
```

The documentation is also visible on [maniac-mc.github.io](https://maniac-mc.github.io).

## Examples and tests

Several example systems are provided in a separate repository:
[mc-topology]([topology-gallery/](https://github.com/maniac-mc/mc-topology)).
The code has been validated against LAMMPS and RASPA for several example cases
located in [mc-topology]([topology-gallery/](https://github.com/maniac-mc/mc-topology)).
Basic adsorption and energy tests are available in the [tests/](tests/) folder.

## Credit

This code was written by Simon Gravelle, who currently maintains the
code, documentation, and the associated [GitHub organization](https://github.com/maniac-mc).
