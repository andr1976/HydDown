.. _dry_ice:

====================
Dry-Ice Estimation
====================

Two dry-ice quantities are estimated: the **atmospheric** (downstream) vapour/solid
split of the released stream, which is always computed, and the **in-vessel** solid
that accumulates when a gas release drags the vessel below the triple point, which is
an opt-in model. All of it is in ``co2_release.py`` (atmospheric) and the two-zone /
solid-regime steppers dispatched from ``hdclass.py``.

Atmospheric (downstream) dry ice
================================

A free jet does no external work, so the released stream conserves its stagnation
enthalpy as it flashes to atmospheric pressure. thermopack's ``two_phase_phflash`` is
not solid-aware, so the vapour/solid split is computed with an explicit **isenthalpic
lever** against pre-computed 1-atm endpoint enthalpies (the pure-vapour
:math:`h_\text{gas}^\text{atm}` and pure-solid :math:`h_\text{solid}^\text{atm}` at the
frost point):

.. math::
   :label: atm-lever

   \beta_\text{gas} = \frac{h_0 - h_\text{solid}^\text{atm}}
                            {h_\text{gas}^\text{atm} - h_\text{solid}^\text{atm}},
   \qquad
   x_\text{solid} = 1 - \operatorname{clip}(\beta_\text{gas},\,0,\,1)

If :math:`\beta_\text{gas} \ge 1` the stream is superheated (no dry ice) and a plain
solid-free flash gives its temperature; otherwise the stream is a vapour + dry-ice
mixture at the **frost point** :math:`T_\text{frost} \approx 194.14` K
(:math:`-79\,^{\circ}`\ C), obtained directly from a psflash whose entropy brackets the
gas/solid values (the gas :math:`\leftrightarrow` solid Gibbs-equality point). The
endpoints depend only on the back pressure and are computed once
(``_init_atm_endpoints``).

An **isentropic** variant (``atm_split_isentropic``) uses the solid-aware psflash
directly and gives the available-work / blast-energy bound.

For a liquid release the discharged stream is typically :math:`\sim 50\,\%` dry ice
by mass at 1 atm - the dry ice a liquid leak produces appears **here**, downstream of
the orifice, not inside the vessel.

**Implementation**: ``co2_release.py:atm_split() / atm_split_isentropic()``.

In-vessel dry ice (opt-in)
==========================

When a **gas** release keeps the liquid in the vessel, the liquid cools to the triple
point, freezes, and a large dry-ice bank is retained; the vessel finally blows down
on the sublimation line. This is enabled with ``release.solid_in_vessel: true``. Once
the tank reaches the triple point (:ref:`thermodynamics`) the vessel is advanced by a
thermopack model with two regimes: a **triple-point plateau** while liquid remains, and
a **sublimation descent** afterwards. Two implementations exist; the default (used for
the CARDICE gas cases) is the two-zone model.

Triple-point plateau (three-phase lever)
----------------------------------------

At the invariant triple point the temperature and pressure are pinned and the specific
volumes and internal energies of the three pure phases (:math:`v_g,v_l,v_s`,
:math:`u_g,u_l,u_s`) are constants. The phase masses then follow from conservation of
total mass, volume and internal energy:

.. math::
   :label: plateau-lever

   \begin{aligned}
   m_g + m_l + m_s &= M \\
   m_g v_g + m_l v_l + m_s v_s &= V \\
   m_g u_g + m_l u_l + m_s u_s &= U
   \end{aligned}

As the leak removes mass and (small) heat enters, :math:`M` and :math:`U` change and
the lever redistributes the phases - the liquid steadily converts to solid at fixed
:math:`T,P`. This produces the observed **pressure plateau** at the triple point.

In the **two-zone** plateau (``two_zone_plateau_step``) the gas is carried as a
separate warm (superheated) zone that keeps its own temperature - heated by the wall
and only weakly coupled to the cold liquid/solid - so the model reproduces the measured
gas superheat *and* the freezing lever simultaneously. The gas fills the vapour volume
at the triple-point pressure.

Sublimation descent
-------------------

Once the liquid is exhausted the state rides the sublimation line. In the two-zone
descent (``two_zone_descent_step``) a warm gas zone leaks and depressurises while the
dry-ice bank sublimes only enough to (a) cool itself as :math:`T_\text{sub}(P)` falls
and (b) absorb the gas :math:`\rightarrow` solid interphase heat
:math:`\dot{Q}_{gs} = \mathrm{UA}_{gs}(T_g - T_s)`:

.. math::
   :label: descent-sub

   \dot{m}_\text{sub} = \underbrace{\frac{m_s\,c_{p,s}\,\bigl(T_s^{-} - T_s\bigr)}{L_\text{sub}}}_{\text{self-cooling}}
   \;+\; \underbrace{\frac{\dot{Q}_{gs}\,\Delta t}{L_\text{sub}}}_{\text{interphase}}

The interphase conductance :math:`\mathrm{UA}_{gs}` (from
``release.solid_h_gas_solid``) is the single lever that trades gas superheat against
retained solid: with :math:`\mathrm{UA}_{gs} = 0` the dry ice is adiabatic and mostly
retained. The gas energy balance carries the wall heat, the outflow enthalpy and the
sublimed-vapour enthalpy; the solid mass is decremented by :math:`\dot{m}_\text{sub}`.

Two-zone parameters
-------------------

.. list-table::
   :widths: 30 70
   :header-rows: 1

   * - Parameter
     - Meaning
   * - ``solid_h_gas_wall``
     - gas :math:`\leftrightarrow` wall HTC below the triple point; a number, or
       ``"calc"`` for a natural-convection estimate (:ref:`heat_transfer`).
   * - ``solid_h_gas_liquid``
     - plateau gas :math:`\leftrightarrow` liquid/solid interphase HTC (small - keeps
       the gas warm; :math:`\approx 0` for a flat liquid surface).
   * - ``solid_h_gas_solid``
     - descent gas :math:`\leftrightarrow` dry-ice interphase conductance (trades gas
       superheat vs retained solid; :math:`\approx 3` W/m\ :sup:`2`\ K calibrated).
   * - ``solid_gas_wall_frac``
     - fraction of the inner wall in gas contact (the rest is wetted).
   * - ``solid_h_inner``
     - wall :math:`\leftrightarrow` solid/liquid HTC for the wetted wall node
       (:ref:`heat_transfer`).

Behaviour
=========

The model reproduces the central experimental contrast:

* **Gas releases** (tests 5/7/9) retain a large in-vessel dry-ice bank
  (:math:`\sim` 270/310/400 kg modelled vs :math:`\sim` 217/298/412 kg measured) with a
  multi-hour triple-point plateau.
* **Liquid releases** (tests 6/8/10) retain essentially none (:math:`\sim` 4-5 kg vs
  :math:`\sim 0`): the liquid drains and boils out **before** the vessel reaches the
  triple point, so it never freezes in place. The dry ice appears downstream instead
  (:eq:`atm-lever`).

See :ref:`validation` for the full comparison.
