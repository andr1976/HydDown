.. _thermodynamics:

==================================
Thermodynamics and the Solid Lever
==================================

This chapter describes how the add-on evaluates CO\ :sub:`2` properties **entirely from
CoolProp** plus a precomputed solid-CO\ :sub:`2` table, the phase regimes, the
self-consistent triple point, the lever-rule flashes below it, and the CoolProp
hand-overs. The lever/table mathematics is in ``co2_release.py``
(:class:`CO2ReleaseModel`); the CoolProp-backed implementation is ``co2_release_cp.py``
(:class:`CO2ReleaseModelCP`), and the solid table is :mod:`hyddown.co2_solid`.

The problem below the triple point
==================================

The atmospheric pressure of CO\ :sub:`2` lies below its triple point, so a release
crosses into the **vapour + solid** region. CoolProp's reference CO\ :sub:`2` equation
of state :cite:`Span1996` cannot be used there directly for two reasons:

* it has **no solid phase** - there is no dry-ice branch to flash onto;
* it is **guarded at the triple point** - a normal ``PropsSI`` call below 5.18 bar in the
  two-phase region either raises or snaps to the saturation boundary.

The add-on solves both without a second runtime backend.

CoolProp gas extrapolation (imposed phase)
------------------------------------------

Pure CO\ :sub:`2` vapour still exists below the triple point (it is the gas in
equilibrium with dry ice). CoolProp's low-level ``AbstractState`` allows the phase to be
**imposed**, which makes the Span-Wagner EoS extrapolate smoothly into the metastable
gas region past the triple-point guard:

.. code-block:: python

   self._gas = AbstractState("HEOS", "CO2")
   self._gas.specify_phase(CP.iphase_gas)     # forced-gas root below the triple point

Every sub-triple **vapour** property (:math:`\rho, u, h, s`, and the HEM flux integrand)
is taken from this forced-gas state. It was verified against an independent GERG-2008
reference to agree to :math:`\sim` 0.01 K in temperature and four decimals in phase
fraction down to :math:`\sim` 1 bar, so the extrapolation is quantitatively sound over
the whole descent.

The solid-CO\ :sub:`2` table
----------------------------

The one thing CoolProp cannot supply is the **solid** branch. It is provided by a small
precomputed table (:mod:`hyddown.co2_solid`) of mass-specific solid properties tabulated
along the sublimation line :math:`P_\text{sub}(T)`:

.. list-table:: Solid-CO\ :sub:`2` table (``co2_solid``)
   :widths: 30 70
   :header-rows: 1

   * - Quantity
     - Role
   * - :math:`P_\text{sub}(T)`, :math:`T_\text{sub}(P)`
     - the sublimation line (solid :math:`+` vapour equilibrium) and its inverse
   * - :math:`v_s(T)`
     - solid specific volume (nearly incompressible; off-line pressure dependence neglected)
   * - :math:`h_s(T),\, s_s(T)`
     - solid enthalpy / entropy in **CoolProp's reference basis** (so they combine
       directly with CoolProp fluid values)
   * - :math:`u_s(T,P) = h_s - P v_s`
     - solid internal energy
   * - :math:`c_{p,s}(T) = \mathrm{d} h_s/\mathrm{d}T`
     - solid heat capacity (for the descent self-cooling term)

The solid values are re-referenced into CoolProp's enthalpy/entropy basis through the
**physical latent heat of sublimation** at each temperature, so a solid property and a
CoolProp gas property can be added on one consistent basis - there is no fixed
reference-offset mismatch. The table is generated offline (see the provenance note in
:ref:`introduction`) and committed to the repository; the **runtime uses only CoolProp
and numpy**.

Why a table is enough: the single degree of freedom
---------------------------------------------------

Below the triple point pure CO\ :sub:`2` is a one-component, two-phase (vapour + solid)
system, so by Gibbs' phase rule it has a **single degree of freedom**: fixing the
pressure fixes the temperature (:math:`T_\text{sub}(P)`) and both pure-phase states.
Every "flash" therefore reduces to **lever-rule arithmetic on the sublimation line** - no
solid-aware equation-of-state solver is needed. This lever reproduces a solid-aware
reference flash essentially exactly, which is why the table + CoolProp combination is
faithful rather than approximate.

Phase regimes
=============

A CO\ :sub:`2` release passes through three regimes:

#. **Above the triple point** - liquid + vapour equilibrium. The vessel is tracked by
   the CoolProp non-equilibrium model (NEM); discharge stagnation states are saturated
   liquid or vapour (native CoolProp).
#. **At the triple point** - solid + liquid + vapour coexist. By Gibbs' phase rule a
   pure three-phase point is invariant: :math:`T` and :math:`P` are pinned while the leak
   and heat only shift the phase *amounts* (the freezing "lever"). This is the pressure
   **plateau** seen in gas-release blowdowns.
#. **Below the triple point** - vapour + solid on the **sublimation line**. Once the
   liquid is exhausted the state rides :math:`T_\text{sub}(P)` down to the back pressure
   (frost point :math:`\approx 194.14` K / :math:`-79\,^{\circ}`\ C at 1 atm).

The self-consistent triple point
================================

The literature triple point is :math:`P_\text{tp} = 5.18` bar,
:math:`T_\text{tp} = 216.59` K. The solid table is built on the **GERG-2008
self-consistent triple point** used to generate it,

.. math::
   :label: triple

   T_\text{tp} = 216.592\ \text{K}, \qquad P_\text{tp} = 5.179\ \text{bar},

so that the solid and fluid branches meet with no spurious latent-heat jump. The
pure-phase specific volumes, enthalpies and internal energies of solid, liquid and vapour
at this point (:math:`v_{g,l,s}`, :math:`u_{g,l,s}`) are stored once as the **vertices of
the three-phase lever** (``_init_triple_point``): the solid vertex from the table, the
liquid and vapour vertices from CoolProp saturation at :math:`P_\text{tp}`. The latent
heat of sublimation follows as :math:`L_\text{sub} = h_g - h_s`.

Hand-overs
==========

Two distinct hand-overs occur.

Vessel state: CoolProp NEM :math:`\rightarrow` two-zone model (at the triple point)
-----------------------------------------------------------------------------------

Above the triple point the vessel inventory is integrated by the CoolProp NEM. Once the
tank pressure reaches :math:`P \le 1.05\,P_\text{tp}` and the opt-in
``release.solid_in_vessel`` is set, the state is handed to the below-triple **two-zone**
model (``hdclass.py``). The NEM's two zones - a warm (superheated) gas zone and a
liquid/solid zone - carry over their masses and internal energies, and a single-zone
:math:`(M_\text{vessel}, U_\text{vessel})` state is also seeded for the terminal
all-solid regime. From then on the below-triple model (:ref:`dry_ice`) advances the
vessel.

Discharge: always CoolProp + solid table
-----------------------------------------

The **discharge rate** is computed on the CoolProp + solid-table basis at *all* pressures
- above the triple point too - so there is no discontinuity in the mass rate at the
triple-point crossing. Only the *vessel state* hands over; the discharge model does not.
A fine-resolution check across the triple point shows the rate is smooth and matches a
choked-orifice value in both regimes (:ref:`discharge`).

.. note::

   The gas *temperature* is not :math:`C^1`-continuous across the vessel hand-over: the
   below-triple descent cools the residual gas a few degrees more steeply than the NEM
   just above (a :math:`\sim 3\,^{\circ}`\ C kink for liquid-release cases). This is small
   compared with the 0-D stratification gap (:ref:`heat_transfer`) and is a known,
   deferred item.

Precomputed property tables
===========================

The below-triple two-zone model steps the gas zone many times with per-step root-finds;
live property calls at every node would be slower than necessary. ``_init_gas_table``
therefore precomputes (once, :math:`\sim` 1 s) from CoolProp + the solid table:

* a **1-D table** of superheated-gas properties (:math:`u, \rho, h`, HEM flux :math:`G`)
  versus :math:`T` at the triple-point pressure - used on the plateau, where :math:`P` is
  fixed;
* a **2-D table** :math:`(T, P) \rightarrow (\rho, u, h)` for the descent, where :math:`P`
  floats below the triple point (the forced-gas root);
* the **sublimation line** :math:`T_\text{sub}(P)` and the saturated-vapour / solid
  internal energies along it.

The descent leak rate itself is evaluated live per step (there are far fewer descent
steps than table nodes). Interpolants (``gas_T_from_u``, ``_gas2d``, ``_T_sub_of_P``,
:math:`\ldots`) provide the state during the integration.
