Unit System
===========

The unit system follows LAMMPS's ``real`` style, with the following conventions:

- Distances are measured in :math:`\text{Å}` (Ångströms).
- Angles are measured in :math:`\text{radians}`.
- Temperature is measured in :math:`\text{K}` (Kelvin).
- Masses are measured in :math:`\text{g·mol}^{-1}`.
- Energies are measured in :math:`\text{kcal·mol}^{-1}`.
- Forces are measured in :math:`\text{kcal·mol}^{-1}\,\text{Å}^{-1}`.
- Charges are expressed in units of the elementary charge :math:`e`.
- Densities are measured in :math:`\text{g·cm}^{-3}`.

All parameter files (e.g. ``.inc`` files specifying ``pair_coeff`` values)
are expected to be written in this **real** unit style. Using parameters from
other unit systems (such as ``metal`` or ``lj``) without conversion will lead
to physically meaningless results.
