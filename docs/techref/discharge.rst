.. _discharge:

===============
Discharge Model
===============

The release rate is a homogeneous-equilibrium (HEM) model with a solid-aware throat,
plus an optional non-equilibrium (HNE) boost for a flashing liquid. It is evaluated on
the CoolProp + solid-table basis (:ref:`thermodynamics`) at all pressures. Implemented in
``co2_release.py`` / ``co2_release_cp.py``.

Stagnation state
================

The upstream (stagnation) state feeding the hole follows ``release.type`` and, once the
liquid is exhausted, the phase flip performed in the mass balance
(``release_phase`` :math:`\;\texttt{liquid} \rightarrow \texttt{gas}`):

* **liquid release** - saturated (or dense) liquid at the tank pressure :math:`P_0`
  (``stagnation(P, "liquid")``);
* **gas release** - saturated vapour at :math:`P_0` (``stagnation(P, "gas")``);
* **below the triple point** - an explicit vapour state on the sublimation line
  (``gas_leak_rate(T, P)``), where a saturated stagnation is undefined.

Each returns the mass-specific stagnation enthalpy :math:`h_0`, entropy :math:`s_0`,
density :math:`\rho_0` and temperature :math:`T_0`.

HEM mass flux
=============

The flow is assumed homogeneous (no slip) and isentropic from the stagnation state to the
throat. Along the isentrope :math:`s = s_0`, the mass flux at throat pressure :math:`P` is

.. math::
   :label: hem-flux

   G(P) = \rho(P)\,\sqrt{2\,\bigl(h_0 - h(P)\bigr)},
   \qquad s(P) = s_0 ,

and the discharge maximises :math:`G` over the throat pressure. The maximiser is the
**choked throat**; if :math:`G` keeps rising down to the back pressure the flow is
unchoked and the throat sits at :math:`P_\text{back}`:

.. math::
   :label: hem-mdot

   \dot{m} = C_d\,A\,G_\text{max},
   \qquad
   G_\text{max} = \max_{P_\text{back} \le P \le P_0} G(P) .

The isentrope properties come from the CoolProp + solid-table state (``_iso_props``):
above the triple point the throat is native CoolProp; below it the throat state is a
**vapour + dry-ice** mixture resolved by the sublimation-line lever, and its
density/enthalpy carry the solid fraction - the "dry-ice lever at the throat".
Consequently a liquid or gas stream that expands below the triple point correctly forms
solid at the orifice.

.. note::

   The flux curve :math:`G(P)` can carry **two local maxima** across the triple-point
   kink (one on the vapour+liquid branch above 5.18 bar, one on the vapour+solid branch
   below it). A bare bounded search can lock onto the wrong one, so
   ``hem_rate_from_stagnation`` first locates the global maximum on a coarse grid and then
   refines locally (``scipy.optimize.minimize_scalar``, bounded).

The returned dictionary carries ``mdot``, ``G``, the throat pressure/temperature, the
throat solid fraction, a ``choked`` flag, and the stagnation state.

**Implementation**: ``co2_release.py:hem_rate() / hem_rate_from_stagnation()``.

Verification against theory
---------------------------

At a representative gas point (saturated vapour, 12 bar, 6 mm, :math:`C_d = 0.68`) the
CoolProp + solid HEM flux matches a textbook **choked ideal-gas orifice** to within a few
percent:

.. list-table::
   :widths: 50 25 25
   :header-rows: 1

   * - Method
     - :math:`G` [kg/m\ :sup:`2`\ s]
     - :math:`\dot{m}` [kg/s]
   * - CoolProp + solid HEM (this model)
     - 3878
     - 0.0746
   * - choked ideal-gas orifice
     - 3755
     - 0.0722

So the HEM flux is physically correct; it is the *feed* (which phase, and the effective
discharge coefficient) that carries the modelling choices.

Non-equilibrium liquid discharge (HNE)
======================================

A **saturated** liquid discharged through a **short, sharp orifice** has too little
residence time to flash to equilibrium: the metastable liquid persists past its
saturation pressure and discharges *above* the equilibrium-HEM rate, toward the frozen
all-liquid limit :cite:`Dyer2007,Fauske1985`. The equilibrium HEM (instantaneous
flashing) is the low bound; the frozen incompressible (Bernoulli) flux is the high bound:

.. math::
   :label: gfrozen

   G_\text{frozen} = \sqrt{2\,\rho_L\,(P_0 - P_\text{back})} .

The add-on interpolates between them with a single non-equilibrium factor
:math:`N \in [0,1]` (``release.liquid_nonequilibrium``):

.. math::
   :label: hne-blend

   G_\text{liq} = \sqrt{(1 - N)\,G_\text{HEM}^2 + N\,G_\text{frozen}^2}

so :math:`N = 0` recovers equilibrium HEM and :math:`N = 1` the frozen liquid. The blend
is applied only to a **liquid** stagnation; a gas discharge (single phase, nothing to
flash) is untouched. The default is :math:`N = 0`, so existing runs are unchanged. For
saturated CO\ :sub:`2` at 10-20 bar the frozen/HEM ratio is large (:math:`\sim 4`),
because HEM chokes hard on the low two-phase sound speed while the dense liquid drives a
much higher Bernoulli flux; a modest :math:`N` therefore produces a substantial boost.

This is a deliberately simple metastable interpolation. A fuller non-equilibrium model
for CO\ :sub:`2` releases is the Vianna/Lopes HNM :cite:`Lopes2018`, derived and validated
for **initially sub-cooled / dense-phase liquid at high pressure** (49-159 bar); its
dominant incompressible term vanishes at the saturated bubble point, so it is not directly
applicable to the CARDICE saturated, low-pressure cases, where the single-factor blend
:eq:`hne-blend` is sufficient. The Leung :math:`\omega`-method :cite:`Leung1986` provides
an alternative subcooled/flashing correlation.

**Implementation**: ``co2_release.py:hem_rate()`` (liquid branch).

Orifice and discharge coefficients
==================================

The GHGT-16 paper :cite:`Drescher2022` (Table 1) reports the actual restriction-orifice
sizes: the pure-CO\ :sub:`2` tests use **3, 4 and 6 mm** holes. With the true sizes fixed,
the discharge coefficient is split by phase and calibrated directly against the raw 1 Hz
Ineris data (:ref:`validation`):

* the **gas** discharge is reliable single-phase HEM (it equals a choked orifice), so it
  takes one :math:`C_{d,\text{gas}}`. Fitting the gas tests (5, 7, 9) gives 0.87-0.92, and
  the coefficient is *constant along each pressure decline* - an independent confirmation
  of the HEM + single-\ :math:`C_d` gas model. The canonical sharp-edged value
  :math:`C_{d,\text{gas}} = 0.84` is adopted (low end of the range);
* the **liquid** discharge takes the two-phase / flashing sharp-orifice coefficient
  :math:`C_{d,\text{liq}} = 0.62` (Darby / API 520 :cite:`Darby2004`) together with the
  non-equilibrium factor :math:`N`.

A liquid release therefore uses :math:`C_{d,\text{liq}}` for the liquid and
:math:`C_{d,\text{gas}}` for its gas tail (``release.discharge_coef_gas``).

.. list-table:: Reconciled discharge (true orifices, calibrated to the 1 Hz data)
   :widths: 10 12 12 12 12 12
   :header-rows: 1

   * - Test
     - Phase
     - :math:`P_0` [bar]
     - Orifice
     - :math:`C_d`
     - :math:`N`
   * - 5
     - gas
     - 20
     - 3 mm
     - 0.84
     - --
   * - 6
     - liquid
     - 20
     - 3 mm
     - 0.62
     - 0.174
   * - 7
     - gas
     - 15
     - 4 mm
     - 0.84
     - --
   * - 8
     - liquid
     - 15
     - 5 mm\ :sup:`*`
     - 0.62
     - 0.175
   * - 9
     - gas
     - 10
     - 4 mm
     - 0.84
     - --
   * - 10
     - liquid
     - 10
     - 4 mm
     - 0.62
     - 0.020

\ :sup:`*` Table 1 lists 4 mm for test 8, but its measured 0.314 kg/s exceeds a controlled
49 bar lab test :cite:`Pursell2012` through the *same* 4 mm hole - impossible, since flow
must rise with upstream pressure. An effective 5 mm orifice reconciles it and drops
:math:`N` from 0.53 to 0.18, in line with test 6.

Pressure-dependent non-equilibrium factor
-----------------------------------------

Reconciling :math:`N` at the fixed :math:`C_{d,\text{liq}} = 0.62` against the measured
steady liquid-drain rate shows that, across the saturated-CO\ :sub:`2` regime, :math:`N`
scales linearly with the distance above the triple point:

.. math::
   :label: n-of-p

   N \approx 0.013\,(P_0 - P_\text{tr})\quad[\text{bar}], \qquad R^2 = 0.88

(:math:`N \to 0` at the triple point, rising to :math:`\sim 0.31` at 28 bar). It is applied
through ``release.liquid_ne_pressure_scaled``, which fades each test's calibrated
:math:`N` as

.. math::

   N_\text{eff}(P) = N\,\operatorname{clip}\!\left(\frac{P - P_\text{tr}}{P_0 - P_\text{tr}},\,0,\,1\right),

so :math:`N = N` at the initial pressure and vanishes near the triple point as the vessel
blows down. Physically the metastable boost is a short-residence, above-triple effect:
test 10 at 10 bar already reads :math:`\sim`\ equilibrium because it sits close to the
triple point.

.. _hole-size:

Hole size, residence time, and full-bore releases
=================================================

The non-equilibrium boost is a **short-residence** phenomenon: it appears when the liquid
does not spend long enough in (and just downstream of) the restriction to nucleate and
flash to equilibrium. The mainstream engineering frameworks parameterise this by **flow
length**, not diameter - Fauske's :math:`\sim 0.1` m "equilibrium length"
:cite:`Fauske1985` and the HNE-DS boiling-delay length threshold both make the flow revert
to HEM once the flow path exceeds :math:`\sim` 0.1 m. The direct CO\ :sub:`2` evidence for
a *diameter* effect at fixed length is thin: choked-flow tests through CO\ :sub:`2`
orifices and nozzles :cite:`Hammer2022` show HEM's under-prediction shrinking toward larger
holes, but confounded by the vena-contracta discharge coefficient, and no clean
fixed-length, varying-diameter data exist above :math:`\sim` 13 mm.

The add-on therefore keeps :math:`N` **independent of hole diameter** (calibrated per case,
faded only with pressure). The responsibility to choose :math:`N` consistently with the
geometry rests with the user, with one important rule:

.. warning::

   **Full-bore / large-diameter ruptures: use** :math:`C_d = 1` **and** :math:`N = 0`
   **(HEM); do not stack the small-orifice boost on** :math:`C_d = 1`. The discharge
   coefficient (a geometric vena-contracta reduction) and the non-equilibrium factor (a
   thermodynamic metastability) capture different physics and **anti-correlate with hole
   size**: a small sharp orifice is low-\ :math:`C_d`, high-\ :math:`N` (metastable); a
   full-bore rupture is :math:`C_d \to 1`, :math:`N \to 0` (long flow path, equilibrium).
   Setting :math:`C_d = 1` **and** keeping a small-orifice :math:`N` takes the favourable
   end of both knobs at once and over-predicts the rate - e.g. at 16 barg saturated liquid,
   :math:`C_d = 1` with the calibrated :math:`N \approx 0.175` gives :math:`\approx 1.9\times`
   the :math:`C_d = 1` HEM rate. Reserve :math:`N > 0` for genuine small restrictions where
   :math:`C_{d,\text{liq}} < 1` is also in play.

Two regimes; HEM vs the HNM
---------------------------

Reconciling :math:`N` across the wider CO\ :sub:`2` release literature reveals **two
regimes that do not share one** :math:`N` **law**:

* **low-pressure saturated** CO\ :sub:`2` (10-30 bar; the HydDown regime), where :math:`N`
  follows :eq:`n-of-p`;
* **high-pressure dense / subcooled** liquid (49-159 bar), where :math:`N` is governed by
  subcooling instead.

For the dense regime, plain HEM with a fixed :math:`C_d` reproduces the measured flux **as
well as** the Vianna homogeneous non-equilibrium model (mean absolute error 20.8 % vs
22.2 %), so the HNM/HRM machinery is unnecessary there - HEM plus the low-pressure
:math:`N(P)` boost is sufficient. For **large leaks** the long residence lets the liquid
flash to equilibrium, so :math:`N \to 0` (pure HEM); the discharge coefficient is
geometry-set and essentially size-independent.
