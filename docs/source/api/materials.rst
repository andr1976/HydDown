materials - Material Property Database
=======================================

The ``materials`` module provides a database of thermal properties for vessel materials
including metals and composite materials commonly used in pressure vessel construction.

.. automodule:: hyddown.materials
   :members:
   :undoc-members:
   :show-inheritance:

Overview
--------

The materials database includes thermal properties needed for heat transfer calculations:

- Thermal conductivity [W/m·K]
- Heat capacity [J/kg·K]
- Density [kg/m³]
- Emissivity [-]

Some properties are temperature-dependent where data is available.

Available Materials
-------------------

Metals
~~~~~~

- **Steel**: Carbon steel, stainless steel variants
- **Aluminum**: Aluminum alloys (common in cryogenic applications)
- **Other metals**: Copper, titanium, etc.

Composite Materials
~~~~~~~~~~~~~~~~~~~

- **Carbon fiber composites**: For Type III/IV vessels
- **Fiberglass**: Glass fiber reinforced polymers
- **Polymer liners**: HDPE, PA6, etc.

Insulation Materials
~~~~~~~~~~~~~~~~~~~~

- **Foam insulation**: Polyurethane, polystyrene
- **Vacuum insulation**: For cryogenic storage
- **Aerogel**: High-performance insulation

Material Properties
-------------------

Each material is defined with the following properties:

Thermal Conductivity
~~~~~~~~~~~~~~~~~~~~

Thermal conductivity :math:`k` [W/m·K] determines heat conduction through the vessel wall:

.. math::

   q = -k \\frac{dT}{dx}

Lower thermal conductivity materials (composites, insulation) slow heat transfer.

Heat Capacity
~~~~~~~~~~~~~

Specific heat capacity :math:`c_p` [J/kg·K] determines thermal mass:

.. math::

   Q = m c_p \\Delta T

Higher heat capacity means more energy is needed to change temperature.

Density
~~~~~~~

Density :math:`\\rho` [kg/m³] combined with heat capacity gives volumetric heat capacity:

.. math::

   C_v = \\rho c_p

Emissivity
~~~~~~~~~~

Emissivity :math:`\\varepsilon` [-] determines radiative heat transfer (important for fire scenarios):

.. math::

   Q_{rad} = \\varepsilon \\sigma A (T_1^4 - T_2^4)

Typical values:

- Polished metal: 0.05-0.15
- Oxidized metal: 0.6-0.8
- Paint/coating: 0.8-0.95

Temperature Dependence
----------------------

For materials with temperature-dependent properties, the module provides interpolation
functions or correlation equations.

Usage in Input Files
--------------------

Materials are referenced by name in the input file:

.. code-block:: yaml

   vessel:
     material: steel
     wall_thickness: 0.01

For composite vessels with multiple layers:

.. code-block:: yaml

   vessel:
     materials:
       - layer: liner
         material: aluminum
         thickness: 0.003
       - layer: overwrap
         material: carbon_fiber
         thickness: 0.015

Adding Custom Materials
------------------------

To add a custom material, define its properties in the materials database following
the existing structure:

.. code-block:: python

   custom_material = {
       'name': 'custom_alloy',
       'thermal_conductivity': 50.0,  # W/m·K
       'heat_capacity': 500.0,         # J/kg·K
       'density': 7850.0,              # kg/m³
       'emissivity': 0.7               # -
   }
