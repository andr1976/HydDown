.. _heat_transfer:

============================
Heat and Mass Transfer
============================

The two-phase (non-equilibrium) vessel carries **separate gas and liquid zones** with
their own temperatures, and a **two-node wall** (a gas-contact node and a condensate
"wetted" node). This chapter documents the wall model above and below the triple point,
the convection coefficients, the initialisation, and the zone mass/energy balances.
Implemented in ``hdclass.py`` (with correlations in ``transport.py`` and the below-triple
wall in the two-zone steppers).

Two-node wall above the triple point (NEM)
==========================================

Above the triple point the CoolProp non-equilibrium model resolves the wall into an
**unwetted (gas-contact)** node :math:`T_\text{wall}` and a **wetted (liquid-contact)**
node :math:`T_\text{wall}^{w}`, exchanging with the respective fluid zone:

.. math::
   :label: qinner

   \dot{Q}_\text{inner} &= h_\text{in}\,A_\text{unwetted}\,\bigl(T_\text{wall} - T_g\bigr) \\
   \dot{Q}_\text{inner}^{w} &= h_\text{in}^{w}\,A_\text{wetted}\,\bigl(T_\text{wall}^{w} - T_l\bigr)

with the wetted/unwetted areas from the vessel geometry at the current liquid level. The
inner coefficients follow natural-convection correlations (``h_inner: calc``); the gas
side uses the **gas** temperature in the NEM. The external side uses ``h_outer`` (the
Jamois-measured :math:`0.39 \approx 0.4` W/m\ :sup:`2`\ K, :ref:`experimental`).

Natural convection (Geankoplis)
===============================

The above-triple natural-convection HTC uses the Geankoplis :cite:`Geankoplis1993`
Nusselt correlation with the Rayleigh number :math:`\mathrm{Ra} = \mathrm{Gr}\,\mathrm{Pr}`:

.. math::
   :label: nu

   \mathrm{Gr} = \frac{g\,\beta\,|T_\text{wall} - T_\text{fluid}|\,L^3}{\nu^2},
   \qquad
   \mathrm{Nu} =
   \begin{cases}
     0.13\,\mathrm{Ra}^{1/3} & \mathrm{Ra} \ge 10^9 \\
     0.59\,\mathrm{Ra}^{1/4} & 10^4 < \mathrm{Ra} < 10^9 \\
     1.36\,\mathrm{Ra}^{1/5} & \text{otherwise}
   \end{cases}
   \qquad h = \frac{\mathrm{Nu}\,k}{L}

with :math:`L` the vessel height (vertical) or diameter (horizontal) and fluid properties
(:math:`\beta, \nu, k, \mathrm{Pr}`) from CoolProp.

Below the triple point (two-zone wall)
======================================

Below the triple point the two-zone model carries the wall as two nodes, and the physics
switches between the **plateau** (liquid present) and the **descent** (dry ice only).

Gas-contact wall (both regimes)
-------------------------------

The gas-contact wall exchanges with the gas and the ambient (the condensate is adiabatic
to it):

.. math::
   :label: gaswall

   m_w c_{p,w}\,\frac{dT_\text{wall}}{dt} = \dot{Q}_\text{out} - \dot{Q}_{wg},
   \quad
   \dot{Q}_{wg} = h_{gw}\,A_g\,\bigl(T_\text{wall} - T_g\bigr),
   \quad
   \dot{Q}_\text{out} = h_\text{out}\,A_\text{out}\,\bigl(T_\text{amb} - T_\text{wall}\bigr),

where :math:`A_g` is the gas-contact area from the vessel geometry and the condensate
volume, and the wall thermal mass / outer area are partitioned between the two nodes by
area fraction so the total is not double-counted.

Plateau wetted wall: boiling into the liquid
--------------------------------------------

While **liquid** is present the condensate-region wall is **boiling-cooled** by the
liquid, not gas-heated. The wall :math:`\rightarrow` liquid heat uses a nucleate
pool-boiling coefficient (below) and flows **into** the liquid/solid zone, boiling some
liquid off:

.. math::
   :label: wetwall-plateau

   m_w^{w} c_{p,w}\,\frac{dT_\text{wall}^{w}}{dt} = \dot{Q}_\text{out}^{w} - \dot{Q}_{wl},
   \quad
   \dot{Q}_{wl} = h_\text{boil}\,A_\text{wet}\,\bigl(T_\text{wall}^{w} - T_\text{tp}\bigr) .

The wall relaxes toward the cold boiling liquid pinned at the triple point, and
:math:`\dot{Q}_{wl}` enters the plateau lever (:ref:`dry_ice`) as extra vaporisation - the
warm wall's stored sensible heat boils liquid off rather than letting it all freeze.

Descent dry-ice wall: sublimation-tracking
------------------------------------------

Once the liquid is exhausted the condensate-region wall is buried against the cold dry-ice
bed and **tracks the sublimation temperature**:

.. math::
   :label: wetwall-descent

   T_\text{wall}^{w} = T_\text{sub}(P) \quad\text{(while dry ice is present)} .

No wall :math:`\rightarrow` solid heat is applied (the ice bed holds the wall on the
sublimation line), so the retained dry-ice mass is unaffected by the wall. This
ice-pinned wall reproduces the CARDICE bottom-sensor cool-down to :math:`\approx
-75\,^{\circ}`\ C.

Phase-depletion clamp
---------------------

When a phase is depleted, its wall node has no fluid to exchange with, so it **clamps to
the surviving phase's wall** (the condensate wall clamps to the gas wall once the dry ice
is gone). This avoids a stranded node drifting to an unphysical temperature.

Churchill-Chu gas-wall convection
=================================

The default below-triple gas :math:`\leftrightarrow` wall coefficient
(``release.solid_h_gas_wall: "churchill"``) is the Churchill-Chu :cite:`ChurchillChu1975`
free-convection correlation for a vertical plate - the same law the
Høydalsvik/Munkejord reference model :cite:`Hoydalsvik2026` uses for the vapour region:

.. math::
   :label: churchill

   \mathrm{Nu} = \left( 0.825 +
   \frac{0.387\,\mathrm{Ra}_{L_v}^{1/6}}
        {\bigl[\,1 + (0.492/\mathrm{Pr})^{9/16}\,\bigr]^{8/27}} \right)^2,
   \qquad h_{gw} = \frac{\mathrm{Nu}\,k}{L_v} .

Two features matter. First, the properties are the CO\ :sub:`2` **vapour** properties at
the film temperature, taken from CoolProp's imposed-gas root (valid below the triple
point). Second, the characteristic length is the **vapour-column height**
:math:`L_v = (1 - \alpha_\ell)\,L_z` above the settled condensate, which shrinks as dry
ice and liquid accumulate; for a sphere/horizontal cylinder (``length = 0``) the vessel
diameter is used. The single-equation Churchill-Chu form is smooth across the
laminar-turbulent transition, unlike the piecewise Geankoplis form :eq:`nu`, which remains
available as ``"calc"``. In the turbulent regime the two agree; the choice only matters in
the low-Rayleigh tail. A fixed number is also accepted. The coefficient falls back to 15
W/m\ :sup:`2`\ K on any non-finite CoolProp film state or a negligible gas/wall
:math:`\Delta T`.

Nucleate boiling on the wetted wall
===================================

The plateau wall :math:`\rightarrow` liquid coefficient (``release.solid_h_inner``) is a
nucleate pool-boiling correlation evaluated at the triple-point saturated liquid:

* ``"cooper"`` (default) - the Cooper :cite:`Cooper1984` reduced-pressure correlation,
  which carries CO\ :sub:`2`'s high reduced pressure natively (the law the SINTEF reference
  model uses);
* ``"calc"`` - the Rohsenow :cite:`Rohsenow1952` correlation with a CO\ :sub:`2`-on-steel
  surface constant (:math:`C_{sf} = 0.013`, :math:`n = 1.7`);
* a fixed number.

Cooper is the default for a physical reason. Its reduced-property form,

.. math::
   :label: cooper

   h_\text{boil} = 55\,P_r^{\,0.12 - 0.2\log_{10} R_p}\,
                   \bigl(-\log_{10} P_r\bigr)^{-0.55}\,M^{-0.5}\,q^{0.67},

needs only the **reduced pressure** :math:`P_r = P/P_c`, the molar mass :math:`M`, the wall
heat flux :math:`q` and the surface roughness :math:`R_p`. It requires **no surface
tension, no saturated liquid/vapour densities and no latent heat** - exactly the properties
that are ill-defined for CO\ :sub:`2` at and below the triple point, where no stable liquid
exists and the metastable saturated-liquid state is poorly characterised. The Rohsenow form
needs all of them (:math:`\sigma, \rho_\ell, \rho_g, h_{fg}, c_{p,\ell}, \mu_\ell`), which
must be evaluated at the pinned triple-point saturated liquid and carry that extrapolation
uncertainty. Cooper sidesteps it with quantities that are well defined at the triple-point
reference.

Both are capped at 3000 W/m\ :sup:`2`\ K: above that the wall is conduction-limited and its
temperature is insensitive to the exact coefficient (a sweep on the CARDICE cases confirms
the retained mass and wall temperature barely move across ``cooper``/``calc``/fixed). This
conduction-limited insensitivity is also why the exact boiling law is not critical - it sets
how *fast* the wall cools, not the sensible-heat budget that fixes the retained dry ice. The
descent uses no boiling law - the dry-ice wall is sublimation-pinned (above).

Gas-liquid interphase
=====================

``calculation.h_gas_liquid`` accepts ``"calc"``/``"calc_two_sided"`` (natural-convection
correlations across the interface) or a fixed number. For a liquid drain a small value
keeps the warm gas decoupled from the cold boiling liquid (near-adiabatic), reproducing the
measured warm-gas / cold-liquid split. The descent gas :math:`\leftrightarrow` dry-ice
interphase is separate and acts over the **bed-top interface area** (:ref:`dry_ice`), not
the wall area.

Initialisation of the warm gas and split wall
=============================================

Two related initialisation features (both keyed on ``initial.gas_temperature``):

* **Superheated gas zone.** With ``initial.gas_temperature`` the gas zone starts warmer
  than the saturated liquid (a lower density and mass; the fill takes up the difference at
  a matched total mass). The value is set to the **mid-band (median) of the top 3 gas-space
  thermocouples** at :math:`t_0` - a representative warm-gas start rather than the single
  warmest sensor.
* **Split wall init.** With a superheated gas zone the wall has equilibrated with the phase
  it touches before the blowdown, so the **gas-contact wall starts at the gas temperature**
  and the **wetted wall at the liquid temperature** (rather than both at
  ``initial.temperature``). An optional ``initial.wall_temperature`` gives a uniform
  override. Default behaviour (no ``gas_temperature``) is unchanged.

Zone mass and energy balances (NEM)
===================================

Each zone integrates its own internal energy. The outflow enthalpy is **attributed to the
phase the mass actually leaves from**: a liquid-space release draws liquid, a vapour-space
release / standard blowdown draws gas (exactly one of :math:`h_g^\text{out}`,
:math:`h_l^\text{out}` is non-zero):

.. math::
   :label: nem-energy

   \Delta U_g &= -\,\dot{m}\,h_g^\text{out}\,\Delta t + \dot{Q}_\text{inner}\,\Delta t
                 - \dot{Q}_{gl}\,\Delta t + E_\text{evap} - E_\text{cond} \\
   \Delta U_l &= -\,\dot{m}\,h_l^\text{out}\,\Delta t + \dot{Q}_\text{inner}^{w}\,\Delta t
                 + \dot{Q}_{gl}\,\Delta t - E_\text{evap} + E_\text{cond}

where :math:`\dot{Q}_{gl} = h_{gl}\,A_\text{interface}\,(T_g - T_l)` is the interphase heat
and :math:`E_\text{evap}, E_\text{cond}` are the boil-off / condensation phase-transfer
energies handled in the mass balance. Both wall :math:`\rightarrow` gas and wall
:math:`\rightarrow` liquid terms are present, so the wall heats both zones.

The gas-temperature gap (stratification)
========================================

For the **gas** releases the modelled 0-D gas node runs colder than the measured warm top
(e.g. :math:`\sim -49` vs :math:`-27\,^{\circ}`\ C for test 5). The real gas is
**stratified** - warm on top, cold near the interface - which a single 0-D gas node cannot
represent; it averages the two and mixes in cold boil-off vapour. This is a fundamental 0-D
limitation, not a heat-transfer-coefficient error (the wall model reproduces the measured
wall, and the gas-wall HTC is the physical free-convection value). Capturing it would need
a stratified (multi-node) gas zone; Vessfire, another 0-D-gas tool, shows the same gap
:cite:`Vaillant2021`.

:numref:`fig-t9-heat` shows the mechanism for a gas release: the gas expansion-cools during
the fast blowdown, then **re-warms** during the constant-pressure plateau as the wall heats
the near-isobaric gas, while the two wall nodes separate into a warm gas-contact wall and a
cold dry-ice-contact wall.

.. _fig-t9-heat:

.. figure:: figures/cardice_t9_heat_diag.png
   :width: 95%

   Test 9 gas-release diagnostic: internal/gas temperatures (left) and inner-wall sensors
   vs the two modelled wall nodes (right).
