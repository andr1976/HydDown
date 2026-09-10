.. _introduction:

============
Introduction
============

Motivation
==========

In carbon capture and storage (CCS) with ship or truck transport, CO\ :sub:`2` is
handled as a pressurised, liquefied fluid at vapour/liquid equilibrium. An
accidental leak - at an onshore tank, on a ship, or during transfer - depressurises
the inventory towards atmospheric pressure. Because the atmospheric pressure of
CO\ :sub:`2` (1 atm) lies **below the triple point** (:math:`P_\text{tp} = 5.18` bar,
:math:`T_\text{tp} = -56.6\,^{\circ}`\ C), the expanding stream necessarily passes into
the **vapour + solid** region: it forms **dry ice ("CO2 snow")** both downstream of the
leak and, for a vapour-space release, as a large retained bank inside the vessel.
Dry-ice formation drives long shutdowns (sublimation is slow) and severe low wall
temperatures that govern material selection.

The base HydDown code :cite:`Andreasen2021` uses CoolProp :cite:`Bell2014` as its
thermodynamic backend. CoolProp's CO\ :sub:`2` equation of state
:cite:`Span1996` is **not valid below the triple point**, so it can neither integrate a
homogeneous-equilibrium (HEM) discharge isentrope through a sub-triple throat nor
evaluate the atmospheric dry-ice end state. The CO\ :sub:`2` release add-on therefore
introduces a **thermopack** :cite:`Wilhelmsen2017` backend with a solid CO\ :sub:`2`
model, used only for the parts of the calculation that reach below the triple point,
while CoolProp continues to drive the vessel thermodynamics above it.

Scope
=====

The add-on is tailored to **pure CO\ :sub:`2`** and to two release geometries, matching
the CARDICE experiments (:ref:`experimental`):

* **Gas (vapour-space) release** - a hole in the vapour space; saturated-vapour
  stagnation. The liquid stays in the vessel, cools to the triple point, freezes to
  a large retained dry-ice bank, and the vessel finally blows down on the
  sublimation line. This is the classic "hundreds of kilograms of snow in the vessel"
  case.
* **Liquid (bottom) release** - a hole in the liquid space; saturated-liquid
  stagnation. The liquid drains (and boils off) and is essentially gone *before* the
  vessel reaches the triple point, so **no significant in-vessel dry ice** forms; the
  dry ice appears downstream of the orifice instead.

Both cases require the two-phase (non-equilibrium) vessel model
(``calculation.non_equilibrium: true``, ``calculation.type: energybalance``), a
liquid fill (``vessel.liquid_level``), and the new top-level ``release:`` block.

What the add-on provides
========================

.. list-table:: Physical models added on ``co2-release-hem``
   :widths: 30 70
   :header-rows: 1

   * - Area
     - Model
   * - Thermodynamics
     - thermopack tc-PR :cite:`LeGuennec2016` with a solid CO\ :sub:`2` model; a
       self-consistent triple point; CoolProp :math:`\leftrightarrow` thermopack
       hand-overs (:ref:`thermodynamics`).
   * - Discharge
     - homogeneous-equilibrium (HEM) rate with a solid-aware throat; an optional
       non-equilibrium (HNE) liquid boost for delayed flashing (:ref:`discharge`).
   * - Dry ice
     - atmospheric (downstream) vapour/solid split by an isenthalpic lever; an
       opt-in in-vessel solid model (single-zone and two-zone plateau/descent)
       (:ref:`dry_ice`).
   * - Heat transfer
     - a two-node (gas-contact / wetted) wall; a natural-convection gas-wall HTC
       valid below the triple point; superheated gas and split wall initialisation
       (:ref:`heat_transfer`).
   * - Mass / energy
     - separate gas- and liquid-zone balances with phase-attributed outflow
       enthalpy for a liquid- vs vapour-space release (:ref:`heat_transfer`).
   * - Validation
     - the full CARDICE pure-CO\ :sub:`2` test set (tests 5-10) against 1 Hz
       Ineris data (:ref:`validation`).

Architecture
============

All thermopack work is encapsulated in a standalone module,
``src/hyddown/co2_release.py`` (class :class:`CO2ReleaseModel`); the main
``HydDown`` class in ``hdclass.py`` only calls thin glue, so thermopack is never
imported for a non-release run. The vessel state above the triple point is tracked
by the existing CoolProp non-equilibrium model (NEM); once the tank approaches the
triple point the state is handed to the thermopack solid-in-vessel model.

A minimal input looks like:

.. code-block:: yaml

   vessel:      {length: 0.0, diameter: 1.58, type: Hemispherical,
                 liquid_level: 0.62, thickness: 0.072, ...}
   initial:     {temperature: 253.15, pressure: 2000000, fluid: "CO2",
                 gas_temperature: 256.2}          # superheated gas zone (optional)
   calculation: {type: energybalance, non_equilibrium: true,
                 time_step: 5, end_time: 13000, h_gas_liquid: 2.0}
   valve:       {type: none, flow: discharge, back_pressure: 101325.}
   release:                                       # <-- the CO2 add-on
     type: liquid                                 # or "gas"
     diameter: 0.00266
     discharge_coef: 0.68
     liquid_nonequilibrium: 0.267                 # HNE boost (liquid only)
     eos: tcPR
     solid_in_vessel: true
     back_pressure: 101325.
     atm_pressure: 101325.
   heat_transfer: {type: specified_h, temp_ambient: 283.15,
                   h_outer: 0.4, h_inner: calc}

Worked examples ship in ``src/hyddown/examples/`` (``co2_release.yml``,
``co2_gas_release.yml``, ``co2_solid_in_vessel.yml``) and the validation cases in
``validation/CARDICE_*.yml``.
