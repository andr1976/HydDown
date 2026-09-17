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

The base HydDown code :cite:`Andreasen2021,Andreasen2026` uses CoolProp :cite:`Bell2014`
as its thermodynamic backend. CoolProp's reference CO\ :sub:`2` equation of state
:cite:`Span1996` has **no solid phase** and is guarded at the triple point, so out of the
box it can neither integrate a homogeneous-equilibrium (HEM) discharge isentrope through
a sub-triple throat nor evaluate the atmospheric dry-ice end state. The CO\ :sub:`2`
release add-on removes both limitations **without a second runtime backend**: CoolProp
is made to extrapolate its fluid EoS into the metastable gas region below the triple
point (imposed-phase evaluation), and a small precomputed solid-CO\ :sub:`2` property
table supplies the one thing CoolProp lacks - the solid branch. Below the triple point
pure CO\ :sub:`2` has a single degree of freedom, so every "flash" reduces to lever-rule
arithmetic on the sublimation line (:ref:`thermodynamics`).

Scope
=====

The add-on is tailored to pure CO\ :sub:`2` and to two release geometries:

* **Gas (vapour-space) release** - a hole in the vapour space; saturated-vapour
  stagnation. The liquid stays in the vessel, cools to the triple point, freezes to
  a large retained dry-ice bank, and the vessel finally blows down on the
  sublimation line. This is the classic "hundreds of kilograms of snow in the vessel"
  case.
* **Liquid (bottom) release** - a hole in the liquid space; saturated- or dense-liquid
  stagnation. The liquid drains (and boils off) and is essentially gone *before* the
  vessel reaches the triple point, so **no significant in-vessel dry ice** forms; the
  dry ice appears downstream of the orifice instead.

Both cases require the two-phase (non-equilibrium) vessel model
(``calculation.non_equilibrium: true``, ``calculation.type: energybalance``), a
liquid fill (``vessel.liquid_level``), and the new top-level ``release:`` block. The
add-on is validated against two independent experimental campaigns - the low-pressure
CARDICE sphere and the high-pressure Høydalsvik/Munkejord dense-phase vessel
(:ref:`experimental`).

What the add-on provides
========================

.. list-table:: Physical models added on ``co2-release-hem``
   :widths: 26 74
   :header-rows: 1

   * - Area
     - Model
   * - Thermodynamics
     - all-CoolProp fluid (native above the triple point; imposed-phase extrapolation
       below it) plus a precomputed solid-CO\ :sub:`2` property table; a
       self-consistent triple point and sublimation line; lever-rule flashes
       (:ref:`thermodynamics`).
   * - Discharge
     - homogeneous-equilibrium (HEM) rate with a solid-aware throat; an optional
       non-equilibrium (HNE) liquid boost :math:`N` for delayed flashing, with a
       pressure-scaled fade; phase-split discharge coefficients (:ref:`discharge`).
   * - Dry ice
     - atmospheric (downstream) vapour/solid split by an isenthalpic lever; an
       opt-in in-vessel solid model (two-zone triple-point plateau and sublimation
       descent, with a mass-derived gas/solid interface) (:ref:`dry_ice`).
   * - Heat transfer
     - a two-node (gas-contact / condensate) wall; Churchill-Chu free convection for
       the below-triple gas wall; Cooper/Rohsenow boiling on the wetted wall during
       the plateau; a sublimation-tracking dry-ice wall on the descent
       (:ref:`heat_transfer`).
   * - Mass / energy
     - separate gas- and liquid-zone balances with phase-attributed outflow
       enthalpy for a liquid- vs vapour-space release (:ref:`heat_transfer`).
   * - Validation
     - the full CARDICE pure-CO\ :sub:`2` set (tests 5-10, 10-20 bar) and the
       Høydalsvik/Munkejord dense-phase set (120 bar, riser and no-riser)
       (:ref:`validation`).

Architecture
============

The CO\ :sub:`2` physics lives in a standalone module, ``src/hyddown/co2_release.py``
(the lever/table/flash mathematics, class :class:`CO2ReleaseModel`) with a
**CoolProp-backed implementation** ``co2_release_cp.py`` (:class:`CO2ReleaseModelCP`)
that supplies every thermodynamic call from CoolProp plus the solid table
:mod:`hyddown.co2_solid`. The runtime therefore needs **only CoolProp and numpy**.
The main ``HydDown`` class in ``hdclass.py`` calls thin glue, so the CO\ :sub:`2`
machinery is never imported for a non-release run. The vessel state above the triple
point is tracked by the existing CoolProp non-equilibrium model (NEM); once the tank
approaches the triple point the state is handed to the below-triple two-zone model.

.. note::

   The solid-CO\ :sub:`2` property table (:mod:`hyddown.co2_solid`) is generated
   **offline, once**, by ``scripts/gen_co2_solid.py`` from a solid-CO\ :sub:`2` model in
   thermopack :cite:`Wilhelmsen2017` (GERG-2008 fluid basis), then re-referenced into
   CoolProp's enthalpy/entropy basis via the physical latent heat of sublimation. The
   table is committed to the repository; thermopack is **not** a runtime dependency.

A minimal input looks like:

.. code-block:: yaml

   vessel:      {length: 0.0, diameter: 1.58, type: Hemispherical,
                 liquid_level: 0.62, thickness: 0.072, ...}
   initial:     {temperature: 253.15, pressure: 2000000, fluid: "CO2",
                 gas_temperature: 256.2}          # superheated gas zone (optional)
   calculation: {type: energybalance, non_equilibrium: true,
                 time_step: 20, end_time: 60000, h_gas_liquid: calc_two_sided}
   valve:       {type: none, flow: discharge, back_pressure: 101325.}
   release:                                       # <-- the CO2 add-on
     type: liquid                                 # or "gas"
     diameter: 0.003
     discharge_coef: 0.62                         # liquid (flashing) orifice
     discharge_coef_gas: 0.84                     # gas tail
     liquid_nonequilibrium: 0.174                 # HNE boost (liquid only)
     eos: CoolProp
     solid_in_vessel: true
     solid_h_gas_wall: churchill                  # below-triple gas-wall convection
     solid_h_inner: cooper                        # plateau boiling on the wetted wall
     solid_h_gas_solid: 3.0                       # gas<->dry-ice interphase HTC
     back_pressure: 101325.
     atm_pressure: 101325.
   heat_transfer: {type: specified_h, temp_ambient: 283.15,
                   h_outer: 0.4, h_inner: calc}

Worked validation cases ship in ``validation/reconciled/CARDICE_*.yml`` (Ineris) and
``validation/munkejord/*.yml`` (SINTEF); the figure generators are
``scripts/validate_cardice_reconciled.py`` and ``scripts/validate_munkejord.py``. A
full list of ``release:`` keys and their defaults is in :ref:`nomenclature`.
