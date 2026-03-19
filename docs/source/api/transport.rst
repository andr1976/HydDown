transport - Heat and Mass Transfer
===================================

The ``transport`` module provides calculations for heat and mass transfer, including
dimensionless numbers, heat transfer coefficients, and mass flow rate calculations
for different valve types.

.. automodule:: hyddown.transport
   :members:
   :undoc-members:
   :show-inheritance:

Dimensionless Numbers
---------------------

The module calculates key dimensionless numbers for heat transfer correlations:

- Grashof number (Gr)
- Prandtl number (Pr)
- Nusselt number (Nu)
- Rayleigh number (Ra)

Heat Transfer Coefficients
---------------------------

Functions for calculating heat transfer coefficients for:

- Natural convection
- Forced convection
- Pool boiling
- Film boiling

Mass Flow Rate Calculations
----------------------------

The module provides mass flow rate functions for different valve types. These functions
are called internally by the main HydDown class based on the valve type specified in
the input file.

Valve Types
-----------

The ``valve.type`` parameter determines the mass flow calculation:

- ``orifice``: Compressible flow through orifice (requires ``diameter``, ``discharge_coef``)
- ``control_valve``: Control valve sizing equation (requires ``Cv``, ``N9``)
- ``relief_valve``: API 520/521 relief valve (requires ``diameter``, ``set_pressure``)
- ``mdot``: Constant mass flow rate (requires ``mass_flow``)

Flow direction is set by ``valve.flow``: ``"discharge"`` or ``"filling"``.
