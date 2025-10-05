Build and Run
=============

Usage
-----

First, clone the repository from GitHub by typing in a terminal:

.. code-block:: bash

   git clone https://github.com/maniac-mc/maniac-mc.git

Then, ``cd`` into the ``maniac-mc/`` folder:

.. code-block:: bash

    cd maniac-mc/

And build the code by runing the following command from a terminal:

.. code-block:: bash

   ./build.sh

This will create the executable in ``./build/maniac``.
The executable ``build.sh`` is a helper that performs the following steps:

1. **Read version information**  
   The script reads the project version from ``version.txt`` and the current
   Git commit ID using ``git rev-parse``.  
   This ensures that every build is traceable.

2. **Generate version module**  
   It processes the template file ``src/version_module.f90.in`` and replaces
   placeholders with the actual version and commit ID, producing
   ``src/version_module.f90``.  
   This Fortran module allows the program to embed its build version.

3. **Clean previous build**  
   Runs ``make distclean`` to remove any previous build artifacts, including
   generated dependencies and include files.

4. **Regenerate dependencies**  
   Executes ``make depend`` so that ``dependencies.d`` is updated. This ensures
   proper compilation order when source files change.

5. **Compile the project**  
   Finally, ``make`` is executed to compile the source code and generate the
   executable in ``./build/maniac``.

Run
---

To run the code, you need to provide three files:

- The main MANIAC input file (``.maniac``) containing the parameters values,
- A LAMMPS data file (``.data``) containing the topology,
- A LAMMPS parameter file (``.inc``), where the ``pair_coeff`` are provided.

To do so, in a terminal, define the paths to each three files:

.. code-block:: bash

   input="input.maniac"
   data="topology.data"
   inc="parameters.inc"

Then, launch **MANIAC** with:

.. code-block:: bash

   path/maniac -i $input -d $data -p $inc

Where ``path`` is the path to the executable ``maniac``. Optionally, you can also
specify:

- An output folder with the ``-o`` flag,
- A molecule reservoir with the ``-r`` flag. The reservoir should also be provided as
a LAMMPS ``.data`` file.

Example:

.. code-block:: bash

    res="reservoir.data"
   ./maniac -i $input -d $data -p $inc -o "outputs/" -r $res

All output files (logs, trajectory, and data) will be written into the specified
folder. When a reservoir is provided, **MANIAC** draws molecules from it and
attempts to insert them into the system. If no reservoir is specified, **MANIAC**
instead replicates the first molecule found in the topology.
