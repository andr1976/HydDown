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

Orifice and discharge coefficient
=================================

Because the per-test orifice size is not reported, it is inferred from the measured
discharge with a single physical coefficient :math:`C_d = 0.68`. The **gas** discharge
is reliable HEM (it equals a choked orifice), so for a gas release, and for the *gas
tail* of a liquid release, the measured gas rate fixes the orifice. The **liquid** rate
of that same test is then matched with the HNE factor :math:`N`. The inferred sizes:

.. list-table:: Estimated orifices (single :math:`C_d = 0.68`)
   :widths: 12 12 14 14 18 30
   :header-rows: 1

   * - Test
     - Phase
     - :math:`P_0` [bar]
     - Orifice
     - :math:`N`
     - Inferred from
   * - 5
     - gas
     - 20
     - 3.44 mm
     - --
     - gas rate
   * - 6
     - liquid
     - 20
     - 2.66 mm
     - 0.267
     - gas tail + liquid rate
   * - 7
     - gas
     - 15
     - 4.38 mm
     - --
     - gas rate
   * - 8
     - liquid
     - 15
     - 4.70 mm
     - 0.195
     - gas tail + liquid rate
   * - 9
     - gas
     - 10
     - 4.45 mm
     - --
     - gas rate
   * - 10
     - liquid
     - 10
     - 4.09 mm
     - 0.0
     - liquid rate

The metastable boost :math:`N` decreases with pressure (0.267 at 20 bar, 0.195 at
15 bar, 0 at 10 bar); at 10 bar the equilibrium HEM liquid flux already meets the
data, so test 10's orifice is taken from the liquid rate directly (the gas-inferred
value would over-drain it). The inferred sizes (2.66-4.70 mm) sit inside the papers'
1 mm-to-full-bore range and near the one documented 4 mm example. A plausible
physical reading is that tests 7-10 share a nominal 4 mm orifice and the spread comes
from using one :math:`C_d`; splitting :math:`C_d` by phase (gas :math:`\sim 0.84`-0.9,
liquid :math:`\sim 0.62`-0.65) would collapse them toward 4 mm - a possible future
refinement.

Why the single-:math:`C_d` + HNE convention
-------------------------------------------

If the orifice is instead inferred from the *liquid* rate assuming equilibrium HEM,
it comes out too large (HEM under-predicts the flashing liquid), and that oversized
orifice then makes the *gas* tail :math:`\sim 2\times` too fast. Inferring the orifice
from the reliable gas discharge and giving the liquid the HNE boost removes this
inconsistency and lets one :math:`C_d` fit both phases (e.g. test 8: liquid
0.311 vs 0.317 kg/s, gas tail 0.030 vs 0.030 kg/s, both matched).
