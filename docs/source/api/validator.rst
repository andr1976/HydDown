validator - Input Validation
=============================

The ``validator`` module provides comprehensive input validation using the Cerberus library.
It defines required and optional parameters for all calculation modes and ensures input
files are correctly structured.

.. automodule:: hyddown.validator
   :members:
   :undoc-members:
   :show-inheritance:

Overview
--------

The validator uses schema-based validation to check:

- Required parameters are present
- Data types are correct
- Values are within acceptable ranges
- Vessel geometry is properly specified
- Initial conditions are physically valid
- Valve parameters are complete for the chosen valve type
- Heat transfer settings match the calculation type

Validation Schema
-----------------

The validation schema is organized into sections matching the input file structure:

- ``vessel``: Geometry and material properties
- ``initial``: Starting pressure, temperature, and fluid
- ``calculation``: Type, time step, and end time
- ``valve``: Flow type, size, and coefficients
- ``heat_transfer``: Heat transfer model (for energybalance calculations)
- ``validation``: Optional validation data for comparison plotting

Input File Sections
-------------------

Vessel Section
~~~~~~~~~~~~~~

Required parameters:

- ``geometry.orientation``: "horizontal" or "vertical"
- ``geometry.length``: Vessel length [m]
- ``geometry.diameter``: Vessel diameter [m]
- ``geometry.type``: "flat-end", "ASME F&D", "DIN", or "Semi-elliptical"
- ``material``: Material name from materials database

Optional parameters:

- ``wall_thickness``: Wall thickness [m]
- ``insulation_thickness``: Insulation thickness [m]

Initial Section
~~~~~~~~~~~~~~~

Required parameters:

- ``pressure``: Initial pressure [Pa]
- ``temperature``: Initial temperature [K]
- ``fluid``: Fluid name (CoolProp format, e.g., "Hydrogen", "HEOS::Hydrogen")

Optional parameters:

- ``liquid_level``: Initial liquid level for two-phase systems [fraction 0-1]

Calculation Section
~~~~~~~~~~~~~~~~~~~

Required parameters:

- ``type``: "isothermal", "isenthalpic", "isentropic", "specified_U", or "energybalance"
- ``time_step``: Time step for integration [s]
- ``end_time``: End time for simulation [s]

Valve Section
~~~~~~~~~~~~~

Required parameters:

- ``type``: "orifice", "control_valve", "relief_valve", or "mdot"
- ``flow``: "discharge" or "filling"

Type-specific parameters:

**Orifice:**

- ``diameter``: Orifice diameter [m]
- ``discharge_coef``: Discharge coefficient [-]

**Control valve:**

- ``Cv``: Valve flow coefficient
- ``N9``: Piping geometry factor

**Relief valve:**

- ``diameter``: Valve diameter [m]
- ``set_pressure``: Set pressure [Pa]
- ``blow_down``: Blow down fraction [-]

**Constant mass flow:**

- ``mass_flow``: Mass flow rate [kg/s]

Heat Transfer Section
~~~~~~~~~~~~~~~~~~~~~

Required for ``calculation.type = energybalance``:

- ``type``: "fixed_U", "fixed_Q", "specified_h", "detailed", or "fire"

Type-specific parameters vary based on the heat transfer model selected.

Error Messages
--------------

The validator provides clear error messages indicating:

- Missing required parameters
- Invalid parameter types
- Out-of-range values
- Incompatible parameter combinations

Example validation error:

.. code-block:: text

   Validation Error: vessel.geometry.orientation must be one of ['horizontal', 'vertical']
   Found: 'diagonal'
