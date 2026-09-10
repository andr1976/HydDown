.. _thermodynamics:

===============================
Thermodynamics and Hand-overs
===============================

This chapter describes the thermopack backend, the CO\ :sub:`2` phase regimes, the
self-consistent triple point, and the hand-overs between CoolProp and thermopack.
All of it is implemented in ``co2_release.py`` (class :class:`CO2ReleaseModel`).

Why a second backend
====================

The atmospheric pressure of CO\ :sub:`2` lies below its triple point, so a release
crosses into the **vapour + solid** region. CoolProp's reference CO\ :sub:`2`
equation of state :cite:`Span1996` is not defined there, so it cannot:

* integrate the HEM discharge isentrope through a throat that drops below
  5.18 bar, nor
* evaluate the atmospheric (1 atm) dry-ice end state.

thermopack :cite:`Wilhelmsen2017` is used for exactly these sub-triple tasks. The
translated-consistent Peng-Robinson EoS (tc-PR) :cite:`LeGuennec2016` is initialised
with a solid CO\ :sub:`2` model:

.. code-block:: python

   self.eos = tcPR("CO2")
   self.eos.init_solid("CO2")

Its solid-aware ``two_phase_psflash`` returns a vapour+solid state (phase index 7,
with the **solid mole fraction in the** ``betaL`` **slot**). Keeping every sub-triple
property inside thermopack means a single, consistent reference basis - the code
never mixes CoolProp and thermopack enthalpies/entropies (which differ by a large,
fixed reference offset). ``eos: tcPR`` is the only supported option because the solid
model is calibrated against it.

Phase regimes
=============

A CO\ :sub:`2` release passes through three regimes:

#. **Above the triple point** - liquid + vapour equilibrium. The vessel is tracked by
   the existing CoolProp non-equilibrium model (NEM); discharge stagnation states are
   saturated liquid or vapour.
#. **At the triple point** - solid + liquid + vapour coexist. By Gibbs' phase rule a
   pure three-phase point is invariant: :math:`T` and :math:`P` are pinned while the
   leak and heat only shift the phase *amounts* (the freezing "lever"). This is the
   pressure **plateau** seen in gas-release blowdowns.
#. **Below the triple point** - vapour + solid on the **sublimation line**. Once the
   liquid is exhausted the state rides :math:`T_\text{sub}(P)` down to the back
   pressure (frost point :math:`\approx 194.14` K / :math:`-79\,^{\circ}`\ C at 1 atm).

The self-consistent triple point
================================

The literature triple point is :math:`P_\text{tp} = 5.18` bar,
:math:`T_\text{tp} = 216.59` K, but tc-PR has its **own** self-consistent triple point,
which the code locates so that the solid and fluid models are thermodynamically
consistent (no spurious latent-heat jumps). The triple point is where the solid and
liquid Gibbs energies are equal along the fluid saturation line (on which
:math:`g_\text{liquid} = g_\text{vapour}` already, for a pure component):

.. math::
   :label: triple

   g_\text{solid}(T, P_\text{sat}(T)) = g_\text{liquid}(T, P_\text{sat}(T)),
   \qquad g = h - T s

Solving (``_init_triple_point``) gives :math:`T_\text{tp}^\text{EoS} = 216.59` K,
:math:`P_\text{tp}^\text{EoS} = 5.255` bar. The pure-phase specific volumes,
enthalpies and internal energies of solid, liquid and vapour at this point are stored
as the vertices of the three-phase lever.

Hand-overs
==========

Two distinct hand-overs occur; understanding them is key to the model.

Vessel state: CoolProp :math:`\rightarrow` thermopack (at the triple point)
---------------------------------------------------------------------------

Above the triple point the vessel inventory is integrated by the CoolProp NEM.
Once the tank pressure reaches :math:`P \le 1.05\,P_\text{tp}^\text{EoS}` and the opt-in
``release.solid_in_vessel`` is set, the state is handed to the thermopack
solid-in-vessel model (``hdclass.py``). The NEM's two zones - a warm (superheated)
gas zone and a liquid/solid zone - are re-based into thermopack's energy basis (masses
are basis-free), and a single-zone :math:`(M_\text{vessel}, U_\text{vessel})` state is
also seeded for the sublimation descent. From then on the below-triple model
(:ref:`dry_ice`) advances the vessel.

Discharge: always thermopack
----------------------------

The **discharge rate** is computed by the thermopack HEM at *all* pressures - above
the triple point too - so there is no discontinuity in the mass rate at the
triple-point crossing. Only the *vessel state* hands over; the discharge model does
not. A fine-resolution check across the triple point shows the rate is smooth and
matches a choked-orifice value in both regimes (:ref:`discharge`).

.. note::

   The gas *temperature* is not :math:`C^1`-continuous across the vessel hand-over:
   the below-triple descent cools the residual gas a few degrees more steeply than the
   NEM just above (a :math:`\sim 3\,^{\circ}`\ C kink for liquid-release cases). This is
   small compared with the 0-D stratification gap (:ref:`heat_transfer`) and is a known,
   deferred item.

Precomputed property tables
===========================

The below-triple two-zone model steps the gas zone many times with per-step
root-finds; live thermopack flashes would be far too slow. ``_init_gas_table``
therefore precomputes (once, :math:`\sim` 1 s):

* a **1-D table** of superheated-gas properties (:math:`u, \rho, h`, HEM flux
  :math:`G`) versus :math:`T` at the triple-point pressure - used on the plateau,
  where :math:`P` is fixed;
* a **2-D table** :math:`(T, P) \rightarrow (\rho, u, h)` for the descent, where
  :math:`P` floats below the triple point;
* the **sublimation line** :math:`T_\text{sub}(P)` and the saturated-vapour /
  solid internal energies along it.

The descent leak rate itself is evaluated live per step (there are far fewer descent
steps than table nodes). Interpolants (``gas_T_from_u``, ``_gas2d``, ``_T_sub_of_P``,
:math:`\ldots`) provide the state during the integration.
