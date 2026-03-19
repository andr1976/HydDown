fire - Fire Heat Load Modeling
===============================

The ``fire`` module provides fire scenario modeling using the Stefan-Boltzmann radiation
approach combined with convection heat transfer.

.. automodule:: hyddown.fire
   :members:
   :undoc-members:
   :show-inheritance:

Predefined Fire Scenarios
--------------------------

The module provides predefined fire scenarios with different heat flux intensities:

Pool Fires
~~~~~~~~~~

- **pool_fire_api521**: 60 kW/m² incident heat flux per API 521
- **pool_fire_scandpower**: 100 kW/m² incident heat flux per Scandpower recommendations

Jet Fires
~~~~~~~~~

- **jet_fire_api521**: 100 kW/m² incident heat flux per API 521
- **jet_fire_scandpower**: 250 kW/m² incident heat flux per Scandpower recommendations

Stefan-Boltzmann Calculation
-----------------------------

Fire heat load is calculated using:

.. math::

   Q = \\varepsilon \\sigma A (T_{fire}^4 - T_{wall}^4) + h_{conv} A (T_{fire} - T_{wall})

where:

- :math:`\\varepsilon` is surface emissivity
- :math:`\\sigma` is the Stefan-Boltzmann constant (5.67×10⁻⁸ W/m²K⁴)
- :math:`A` is the exposed surface area
- :math:`T_{fire}` is the fire temperature
- :math:`T_{wall}` is the vessel wall temperature
- :math:`h_{conv}` is the convective heat transfer coefficient

Usage Example
-------------

To use a fire scenario in your input file:

.. code-block:: yaml

   heat_transfer:
     type: fire
     fire_type: pool_fire_api521
     emissivity: 0.8
