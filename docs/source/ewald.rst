Ewald summation
===============

Overview
--------

This section describes the Ewald summation implementation used for computing
long-range electrostatic interactions in periodic systems. The code is
organized into multiple Fortran modules that handle both real-space and
reciprocal-space contributions.

The Ewald method splits the Coulomb potential into two components:

1. Real-space sum (short-range): rapidly decaying interactions computed
   directly between particles.
2. Reciprocal-space sum (long-range): Fourier-space computation using
   precomputed structure factors and k-vectors.

A self-interaction correction ensures that each particle does not
spuriously interact with its own periodic images.

Short range
###########

The short-range (real-space) contribution for a pair of charges :math:`q_i` and
:math:`q_j` separated by distance :math:`r_{ij}` is given by:

.. math::

    E_{ij}^\text{real} = q_i q_j \frac{\text{erfc}(\alpha r_{ij})}{r_{ij}}

where :math:`\alpha` is the Ewald screening parameter.

The total real-space energy for a single molecule or residue is obtained
by summing over all unique atom pairs within the molecule:

.. math::

   E_\text{real} = \sum_{i<j} q_i q_j \frac{\text{erfc}(\alpha r_{ij})}{r_{ij}}

Self-interactions (when :math:`i=j`) are excluded to prevent spurious
contributions.

Long-range
----------

The long-range (reciprocal-space) contribution is computed in Fourier space.
For a reciprocal lattice vector :math:`\mathbf{k}` with components
:math:`k_x, k_y, k_z`, the structure factor amplitude is:

.. math::

   A(\mathbf{k}) = \sum_i q_i e^{i \mathbf{k} \cdot \mathbf{r}_i}

The reciprocal-space energy is then computed as:

.. math::

   E_\text{recip} = \sum_{\mathbf{k} \neq 0} W(\mathbf{k}) \, |A(\mathbf{k})|^2 \, f_\text{sym}(\mathbf{k})

where 

- :math:`W(\mathbf{k}) = \frac{\exp(-|\mathbf{k}|^2 / 4\alpha^2)}{|\mathbf{k}|^2}` is the precomputed reciprocal-space weight,
- :math:`f_\text{sym}(\mathbf{k})` accounts for symmetry between :math:`\mathbf{k}` and :math:`-\mathbf{k}`.

All structure factors :math:`A(\mathbf{k})` and weights :math:`W(\mathbf{k})` are
precomputed and stored, allowing efficient energy evaluation without repeatedly
looping over all atoms for each k-vector.

Self-interaction correction
--------------------------

Each charge interacts with its periodic images, including itself. This
spurious self-energy is subtracted:

.. math::

   E_\text{self} = - \frac{\alpha}{\sqrt{\pi}} \sum_i q_i^2

Total energy
------------

The total Coulomb energy is the sum of the real-space, reciprocal-space, and
self-interaction terms:

.. math::

   E_\text{Coulomb} = E_\text{real} + E_\text{recip} + E_\text{self}