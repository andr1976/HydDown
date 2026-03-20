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

- **api_pool**: ~60 kW/m² incident heat flux per API 521 (pool fire)
- **api_jet**: ~100 kW/m² incident heat flux per API 521 (jet fire)
- **scandpower**: Higher intensity fire per Scandpower guidelines

Stefan-Boltzmann Calculation
-----------------------------

Fire heat load is calculated using:

.. math::

   Q = \\varepsilon \\sigma A (T_{fire}^4 - T_{wall}^4) + h_{conv} A (T_{fire} - T_{wall})

where:

- :math:`\\varepsilon` is surface emissivity (built-in value)
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
     type: "s-b"        # Stefan-Boltzmann radiation model
     fire: "api_pool"   # API 521 pool fire scenario
