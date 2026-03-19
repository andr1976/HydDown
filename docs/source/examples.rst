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
     geometry:
       orientation: vertical
       length: 2.0
       diameter: 0.5
       type: flat-end
     material: steel

   initial:
     pressure: 150e5
     temperature: 293.15
     fluid: Hydrogen

   calculation:
     type: isentropic
     time_step: 0.1
     end_time: 100

   valve:
     type: orifice
     flow: discharge
     diameter: 0.01
     discharge_coef: 0.84

Vessel Filling
--------------

Filling with Heat Transfer
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Model vessel filling with heat transfer to surroundings:

.. code-block:: yaml

   vessel:
     geometry:
       orientation: vertical
       length: 2.0
       diameter: 0.5
       type: flat-end
     material: steel
     wall_thickness: 0.01

   initial:
     pressure: 1e5
     temperature: 293.15
     fluid: Hydrogen

   calculation:
     type: energybalance
     time_step: 0.1
     end_time: 300

   valve:
     type: control_valve
     flow: filling
     Cv: 0.1
     N9: 0.12

   heat_transfer:
     type: specified_h
     h_inner: 100
     h_outer: 10
     T_ambient: 293.15

Fire Scenario
-------------

Pool Fire Heat Load
~~~~~~~~~~~~~~~~~~~

Simulate vessel response to a pool fire:

.. code-block:: yaml

   vessel:
     geometry:
       orientation: horizontal
       length: 5.0
       diameter: 1.0
       type: ASME F&D
     material: steel
     wall_thickness: 0.02

   initial:
     pressure: 100e5
     temperature: 293.15
     fluid: Hydrogen

   calculation:
     type: energybalance
     time_step: 0.1
     end_time: 600

   valve:
     type: relief_valve
     flow: discharge
     diameter: 0.05
     set_pressure: 120e5
     blow_down: 0.1

   heat_transfer:
     type: fire
     fire_type: pool_fire_api521
     emissivity: 0.8

Relief Valve Sizing
-------------------

API 521 Relief Valve
~~~~~~~~~~~~~~~~~~~~

Size a relief valve for fire scenario per API 521:

.. code-block:: yaml

   vessel:
     geometry:
       orientation: horizontal
       length: 10.0
       diameter: 2.0
       type: Semi-elliptical
     material: steel
     wall_thickness: 0.03

   initial:
     pressure: 80e5
     temperature: 293.15
     fluid: Methane

   calculation:
     type: energybalance
     time_step: 0.5
     end_time: 1800

   valve:
     type: relief_valve
     flow: discharge
     diameter: 0.1
     set_pressure: 85e5
     blow_down: 0.1
     back_pressure: 1e5

   heat_transfer:
     type: fire
     fire_type: jet_fire_api521
     emissivity: 0.9

Two-Phase Calculations
----------------------

Liquid Hydrogen Storage
~~~~~~~~~~~~~~~~~~~~~~~~

Model two-phase hydrogen storage with boiling:

.. code-block:: yaml

   vessel:
     geometry:
       orientation: vertical
       length: 3.0
       diameter: 1.0
       type: flat-end
     material: aluminum
     wall_thickness: 0.015

   initial:
     pressure: 2e5
     temperature: 25.0
     fluid: Hydrogen
     liquid_level: 0.5

   calculation:
     type: energybalance
     time_step: 1.0
     end_time: 3600

   valve:
     type: mdot
     flow: discharge
     mass_flow: 0.01

   heat_transfer:
     type: specified_h
     h_inner_gas: 50
     h_inner_liquid: 500
     h_outer: 5
     T_ambient: 293.15
