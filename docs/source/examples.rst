Examples
========

This page provides examples of common HydDown calculations.

Vessel Depressurization
-----------------------

Isentropic Depressurization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Calculate adiabatic depressurization of a hydrogen vessel through an orifice:

.. code-block:: yaml

   vessel:
     length: 2.0
     diameter: 0.5
     orientation: "vertical"
     type: "Flat-end"

   initial:
     pressure: 15000000
     temperature: 293.15
     fluid: "Hydrogen"

   calculation:
     type: "isentropic"
     time_step: 0.1
     end_time: 100

   valve:
     flow: "discharge"
     type: "orifice"
     diameter: 0.01
     discharge_coef: 0.84
     back_pressure: 101325

Vessel Filling
--------------

Filling with Heat Transfer
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Model vessel filling with heat transfer to surroundings:

.. code-block:: yaml

   vessel:
     length: 2.0
     diameter: 0.5
     thickness: 0.01
     heat_capacity: 500
     density: 7800
     orientation: "vertical"

   initial:
     pressure: 100000
     temperature: 293.15
     fluid: "Hydrogen"

   calculation:
     type: "energybalance"
     time_step: 0.1
     end_time: 300

   valve:
     flow: "filling"
     type: "controlvalve"
     Cv: 0.1
     back_pressure: 20000000

   heat_transfer:
     type: "specified_h"
     h_inner: 100
     h_outer: 10
     temp_ambient: 293.15

Fire Scenario
-------------

Pool Fire Heat Load
~~~~~~~~~~~~~~~~~~~

Simulate vessel response to a pool fire:

.. code-block:: yaml

   vessel:
     length: 5.0
     diameter: 1.0
     thickness: 0.02
     heat_capacity: 500
     density: 7800
     orientation: "horizontal"
     type: "ASME F&D"

   initial:
     pressure: 10000000
     temperature: 293.15
     fluid: "Hydrogen"

   calculation:
     type: "energybalance"
     time_step: 0.1
     end_time: 600

   valve:
     flow: "discharge"
     type: "relief"
     set_pressure: 12000000
     back_pressure: 101325

   heat_transfer:
     type: "s-b"
     fire: "api_pool"

Relief Valve Sizing
-------------------

API 521 Relief Valve
~~~~~~~~~~~~~~~~~~~~

Size a relief valve for fire scenario per API 521:

.. code-block:: yaml

   vessel:
     length: 10.0
     diameter: 2.0
     thickness: 0.03
     heat_capacity: 500
     density: 7800
     orientation: "horizontal"
     type: "Hemispherical"

   initial:
     pressure: 8000000
     temperature: 293.15
     fluid: "Methane"

   calculation:
     type: "energybalance"
     time_step: 0.5
     end_time: 1800

   valve:
     flow: "discharge"
     type: "relief"
     set_pressure: 8500000
     back_pressure: 100000

   heat_transfer:
     type: "s-b"
     fire: "api_jet"

Two-Phase Calculations
----------------------

Liquid Hydrogen Storage
~~~~~~~~~~~~~~~~~~~~~~~~

Model two-phase hydrogen storage with boiling:

.. code-block:: yaml

   vessel:
     length: 3.0
     diameter: 1.0
     thickness: 0.015
     heat_capacity: 900
     density: 2700
     thermal_conductivity: 200
     orientation: "vertical"
     type: "Flat-end"
     liquid_level: 0.5

   initial:
     pressure: 200000
     temperature: 25.0
     fluid: "Hydrogen"

   calculation:
     type: "energybalance"
     time_step: 1.0
     end_time: 3600

   valve:
     flow: "discharge"
     type: "mdot"
     mdot: 0.01
     back_pressure: 101325

   heat_transfer:
     type: "specified_h"
     h_inner: 50
     h_outer: 5
     temp_ambient: 293.15
