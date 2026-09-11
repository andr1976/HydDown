.. _discharge:

===============
Discharge Model
===============

The release rate is computed by a homogeneous-equilibrium (HEM) model with a
solid-aware throat, plus an optional non-equilibrium (HNE) boost for a flashing
liquid. Implemented in ``co2_release.py``.

Stagnation state
================

The upstream (stagnation) state feeding the hole follows ``release.type`` and,
once the liquid is exhausted, the phase flip performed in the mass balance
(``release_phase`` :math:`\;\texttt{liquid} \rightarrow \texttt{gas}`):

* **liquid release** - saturated liquid at the tank pressure :math:`P_0`
  (``stagnation(P, "liquid")``);
* **gas release** - saturated vapour at :math:`P_0` (``stagnation(P, "gas")``);
* **below the triple point** - an explicit vapour state on the sublimation line
  (``gas_leak_rate(T, P)``), where a saturated stagnation is undefined.

Each returns the mass-specific stagnation enthalpy :math:`h_0`, molar entropy
:math:`s_0`, density :math:`\rho_0` and temperature :math:`T_0`.

HEM mass flux
=============

The flow is assumed homogeneous (no slip) and isentropic from the stagnation state
to the throat. Along the isentrope :math:`s = s_0`, the mass flux at throat pressure
:math:`P` is

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

The isentrope properties come from thermopack's **solid-aware** ``two_phase_psflash``
(``_iso_props``): below the triple pressure the throat state is a vapour + dry-ice
mixture, and its density/enthalpy carry the solid fraction - the "dry-ice lever at the
throat". Consequently a liquid or gas stream that expands below the triple point
correctly forms solid at the orifice.

.. note::

   The flux curve :math:`G(P)` can carry **two local maxima** across the triple-point
   kink (one on the vapour+liquid branch above 5.18 bar, one on the vapour+solid branch
   below it). A bare bounded search can lock onto the wrong one, so
   ``hem_rate_from_stagnation`` first locates the global maximum on a coarse
   40-point grid and then refines locally (``scipy.optimize.minimize_scalar``, bounded).

The returned dictionary carries ``mdot``, ``G``, the throat pressure/temperature,
the throat solid fraction, a ``choked`` flag, and the stagnation state.

**Implementation**: ``co2_release.py:hem_rate() / hem_rate_from_stagnation()``.

Verification against theory
---------------------------

At a representative gas point (saturated vapour, 12 bar, 6 mm, :math:`C_d = 0.68`)
the thermopack HEM, a CoolProp HEOS HEM, and a textbook **choked ideal-gas orifice**
agree to within a few percent:

.. list-table::
   :widths: 40 30 30
   :header-rows: 1

   * - Method
     - :math:`G` [kg/m\ :sup:`2`\ s]
     - :math:`\dot{m}` [kg/s]
   * - thermopack tc-PR HEM
     - 3861
     - 0.0742
   * - CoolProp HEOS HEM
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
saturation pressure and discharges *above* the equilibrium-HEM rate, toward the
frozen all-liquid limit :cite:`Dyer2007,Fauske1985`. The equilibrium HEM
(instantaneous flashing) is the low bound; the frozen incompressible (Bernoulli)
flux is the high bound:

.. math::
   :label: gfrozen

   G_\text{frozen} = \sqrt{2\,\rho_L\,(P_0 - P_\text{back})} .

The add-on interpolates between them with a single non-equilibrium factor
:math:`N \in [0,1]` (``release.liquid_nonequilibrium``):

.. math::
   :label: hne-blend

   G_\text{liq} = \sqrt{(1 - N)\,G_\text{HEM}^2 + N\,G_\text{frozen}^2}

so :math:`N = 0` recovers equilibrium HEM and :math:`N = 1` the frozen liquid. The
blend is applied only to a **liquid** stagnation; a gas discharge (single phase,
nothing to flash) is untouched. The default is :math:`N = 0`, so existing runs are
unchanged.

This is a deliberately simple metastable interpolation. A fuller non-equilibrium
model for CO\ :sub:`2` releases is the Vianna/Lopes HNM :cite:`Lopes2018`, which adds a
bubble-nucleation relaxation factor across three isentropic segments (liquid
:math:`\rightarrow` saturation :math:`\rightarrow` triple point :math:`\rightarrow`
atmosphere). That model is derived and validated for **initially sub-cooled / dense-phase
liquid at high pressure** (49-159 bar) and explicitly excludes gas-space releases; its
dominant incompressible term vanishes at the saturated bubble point. It is therefore not
directly applicable to the CARDICE saturated, low-pressure (10-20 bar) cases, where the
single-factor blend :eq:`hne-blend` is sufficient. The Leung :math:`\omega`-method
:cite:`Leung1986` provides an alternative subcooled/flashing correlation.

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
steady liquid-drain rate shows that, across the saturated-CO\ :sub:`2` regime (CARDICE plus
the Ineris single point :cite:`Drescher2022`), :math:`N` scales linearly with the distance
above the triple point:

.. math::
   :label: n-of-p

   N \approx 0.013\,(P_0 - P_\text{tr})\quad[\text{bar}], \qquad R^2 = 0.88

(a free power exponent returns 1.02, i.e. linear; :math:`N \to 0` at the triple point,
rising to :math:`\sim 0.31` at 28 bar). It is applied through
``release.liquid_ne_pressure_scaled``, which fades each test's calibrated :math:`N` as

.. math::

   N_\text{eff}(P) = N\,\operatorname{clip}\!\left(\frac{P - P_\text{tr}}{P_0 - P_\text{tr}},\,0,\,1\right),

so :math:`N = N` at the initial pressure and vanishes near the triple point as the vessel
blows down. Physically the metastable boost is a short-residence, above-triple effect:
test 10 at 10 bar already reads :math:`\sim`\ equilibrium because it sits close to the
triple point.

Two regimes; HEM vs the HNM
---------------------------

Reconciling :math:`N` across the wider CO\ :sub:`2` release literature - CARDICE, Ineris,
the saturated Toesse jets :cite:`Toesse2013`, the Pursell orifice tests :cite:`Pursell2012`
and the Vianna/Lopes HNM dataset :cite:`Lopes2018` - reveals **two regimes that do not
share one** :math:`N` **law**:

* **low-pressure saturated** CO\ :sub:`2` (10-30 bar; the HydDown regime), where :math:`N`
  follows :eq:`n-of-p`;
* **high-pressure dense / subcooled** liquid (49-159 bar), where :math:`N` is governed by
  subcooling instead.

For the dense regime, plain HEM with a fixed :math:`C_d` reproduces the measured flux **as
well as** the Vianna homogeneous non-equilibrium model (mean absolute error 20.8 % vs
22.2 %), so the HNM/HRM machinery is unnecessary there - HEM plus the low-pressure
:math:`N(P)` boost is sufficient. The residence/diameter dependence of :math:`N` (smaller
orifice :math:`\to` higher :math:`N`) is real but cannot be pinned universally because the
orifice *lengths* are unreported. For **large leaks** the long residence lets the liquid
flash to equilibrium, so :math:`N \to 0` (pure HEM); the discharge coefficient is
geometry-set and essentially size-independent.
