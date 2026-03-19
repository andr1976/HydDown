Parameter Reference
===================

This page explains the less common or industry-specific parameters used in HydDown input files.

Control Valve Parameters
-------------------------

Cv (Flow Coefficient)
~~~~~~~~~~~~~~~~~~~~~

**Parameter**: ``valve.Cv``

**Description**: The flow coefficient (Cv) is a standard measure of control valve capacity defined by ANSI/ISA and IEC 60534 standards.

**Definition**: Cv is the number of US gallons per minute of water at 60°F that will flow through the valve with a pressure drop of 1 psi.

**Units**: US gal/min/psi^0.5

**Typical Values**:
- Small valves: Cv = 0.1 - 10
- Medium valves: Cv = 10 - 100
- Large valves: Cv = 100 - 1000+

**Example**:

.. code-block:: yaml

   valve:
     type: "controlvalve"    # Note: one word, no underscore!
     Cv: 0.1                  # Small control valve
     flow: "filling"
     back_pressure: 20000000  # Required!

**References**:
- ANSI/ISA-75.01 (Flow Equations for Sizing Control Valves)
- IEC 60534-2-1

**Note**: If you don't know the Cv, you can estimate from valve size or use manufacturer data sheets.

**Important**: The valve type must be ``"controlvalve"`` (one word), not ``"control_valve"``.

Fire Heat Transfer Parameters
------------------------------

emissivity
~~~~~~~~~~

**Parameter**: ``heat_transfer.emissivity``

**Description**: Surface emissivity (ε) is a measure of how effectively a surface emits thermal radiation compared to a perfect black body.

**Range**: 0.0 to 1.0 (dimensionless)

**Physical Meaning**:
- ε = 0: Perfect reflector (no radiation emitted)
- ε = 1: Perfect black body (maximum radiation)
- Real surfaces: 0.05 < ε < 0.95

**Typical Values**:

==================  ====================
Surface Type        Emissivity (ε)
==================  ====================
Polished aluminum   0.05 - 0.10
Polished steel      0.07 - 0.15
Oxidized steel      0.60 - 0.80
Painted surfaces    0.80 - 0.95
Carbon steel        0.70 - 0.80
Stainless steel     0.15 - 0.30 (polished)
                    0.40 - 0.60 (oxidized)
==================  ====================

**Used In**: Stefan-Boltzmann radiation heat transfer equation:

.. math::

   Q = \\varepsilon \\sigma A (T_{fire}^4 - T_{wall}^4)

Where:
- ε = emissivity
- σ = Stefan-Boltzmann constant (5.67×10⁻⁸ W/m²K⁴)
- A = surface area
- T = temperature (K)

**Example**:

.. code-block:: yaml

   heat_transfer:
     type: "fire"
     fire_type: "pool_fire_api521"
     emissivity: 0.8    # Painted or oxidized steel

**Recommendation**: Use 0.7-0.9 for most fire scenarios with steel vessels (conservative).

fire_type
~~~~~~~~~

**Parameter**: ``heat_transfer.fire_type``

**Description**: Predefined fire scenario with specific heat flux intensity.

**Allowed Values**:

================================  =================  ====================
Fire Type                         Heat Flux          Standard/Reference
================================  =================  ====================
``pool_fire_api521``              60 kW/m²           API 521
``pool_fire_scandpower``          100 kW/m²          Scandpower
``jet_fire_api521``               100 kW/m²          API 521
``jet_fire_scandpower``           250 kW/m²          Scandpower
================================  =================  ====================

**Pool Fires**: Liquid fuel burning in a pool (lower heat flux)
**Jet Fires**: High-pressure gas release ignited (higher heat flux)

**Example**:

.. code-block:: yaml

   heat_transfer:
     type: "fire"
     fire_type: "pool_fire_api521"    # Conservative 60 kW/m²
     emissivity: 0.8

Valve Back Pressure
-------------------

back_pressure
~~~~~~~~~~~~~

**Parameter**: ``valve.back_pressure``

**Description**: The pressure downstream of the valve (discharge pressure).

**Units**: Pa (Pascals)

**Physical Meaning**:
- For **discharge**: Ambient pressure where gas exits
- For **filling**: Supply pressure feeding the vessel

**Typical Values**:

======================  ====================
Condition               back_pressure (Pa)
======================  ====================
Atmospheric discharge   101325 (1 atm)
Elevated discharge      > 101325
Vacuum discharge        < 101325
Filling from supply     Supply line pressure
======================  ====================

**Example (Discharge to Atmosphere)**:

.. code-block:: yaml

   valve:
     flow: "discharge"
     type: "orifice"
     back_pressure: 101325    # Atmospheric pressure

**Example (Filling)**:

.. code-block:: yaml

   valve:
     flow: "filling"
     type: "control_valve"
     back_pressure: 20000000  # 200 bar supply pressure

**Important**: This is a **required field** for most valve types!

Relief Valve Parameters
-----------------------

set_pressure
~~~~~~~~~~~~

**Parameter**: ``valve.set_pressure``

**Description**: Pressure at which the relief valve opens (set point).

**Units**: Pa (Pascals)

**Typical Range**: 110-120% of vessel design pressure

**Example**:

.. code-block:: yaml

   vessel:
     # Design pressure: 100 bar = 10,000,000 Pa

   valve:
     type: "relief_valve"
     set_pressure: 11000000    # Opens at 110 bar (110% of design)

blow_down (or blowdown)
~~~~~~~~~~~~~~~~~~~~~~~

**Parameter**: ``valve.blow_down`` or ``valve.blowdown``

**Description**: Pressure reduction (as fraction of set pressure) before the valve reseats.

**Range**: 0.0 to 1.0 (typically 0.05 to 0.2)

**Physical Meaning**:
- Relief valve opens at set_pressure
- Closes at: set_pressure × (1 - blow_down)

**Typical Values**:
- Standard PRVs: 0.10 (10% blowdown)
- Fire scenario: 0.05 - 0.15

**Example**:

.. code-block:: yaml

   valve:
     set_pressure: 11000000    # Opens at 110 bar
     blow_down: 0.1            # Closes at 99 bar (110 - 11)

**Calculation**: Reseat pressure = 11 MPa × (1 - 0.1) = 9.9 MPa

Common Units and Conversions
-----------------------------

Pressure
~~~~~~~~

===============  ==================
Unit             Value in Pa
===============  ==================
1 Pa             1
1 kPa            1,000
1 bar            100,000
1 MPa            1,000,000
1 atm            101,325
1 psi            6,894.76
===============  ==================

**Example**: 150 bar = 15,000,000 Pa

Temperature
~~~~~~~~~~~

- Always use **Kelvin (K)** in input files
- Conversion: K = °C + 273.15
- Example: 20°C = 293.15 K

References
----------

For more details on these parameters, see:

- **Control Valves**: ANSI/ISA-75.01, IEC 60534
- **Fire Scenarios**: API 521, Scandpower guidelines
- **Emissivity**: Engineering handbooks, material datasheets
- **Relief Valves**: API 520/521, ASME Section VIII
