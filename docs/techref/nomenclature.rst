.. _nomenclature:

=============================
Nomenclature and Input Schema
=============================

Principal symbols
=================

.. list-table::
   :widths: 20 55 25
   :header-rows: 1

   * - Symbol
     - Meaning
     - Units
   * - :math:`G`
     - mass flux (per unit throat area)
     - kg/m\ :sup:`2`/s
   * - :math:`\dot m`
     - mass flow rate
     - kg/s
   * - :math:`C_d`
     - discharge coefficient (:math:`C_{d,\text{gas}}`, :math:`C_{d,\text{liq}}`)
     - --
   * - :math:`N`
     - non-equilibrium (HNE) factor, 0 = HEM, 1 = frozen
     - --
   * - :math:`P_0,\,T_0,\,\rho_0,\,h_0,\,s_0`
     - stagnation pressure / temperature / density / enthalpy / entropy
     - Pa, K, kg/m\ :sup:`3`, J/kg, J/kg/K
   * - :math:`P_\text{tp},\,T_\text{tp}`
     - triple point (5.179 bar, 216.592 K)
     - Pa, K
   * - :math:`T_\text{sub}(P)`
     - sublimation temperature on the vapour+solid line
     - K
   * - :math:`L_\text{sub}`
     - latent heat of sublimation (:math:`h_g - h_s`)
     - J/kg
   * - :math:`m_g,\,m_l,\,m_s`
     - gas / liquid / solid masses in the vessel
     - kg
   * - :math:`v_{g,l,s},\,u_{g,l,s},\,h_{g,l,s}`
     - pure-phase specific volume / internal energy / enthalpy
     - m\ :sup:`3`/kg, J/kg
   * - :math:`\alpha_\ell`
     - liquid volume fraction
     - --
   * - :math:`L_v`
     - vapour-column height above the condensate (Churchill-Chu length)
     - m
   * - :math:`\mathrm{Ra},\mathrm{Nu},\mathrm{Pr},\mathrm{Gr}`
     - Rayleigh / Nusselt / Prandtl / Grashof numbers
     - --
   * - :math:`h_{gw}`
     - below-triple gas :math:`\leftrightarrow` wall HTC
     - W/m\ :sup:`2`/K
   * - :math:`h_\text{boil}`
     - plateau wall :math:`\leftrightarrow` boiling-liquid HTC
     - W/m\ :sup:`2`/K
   * - :math:`h_{gs}`
     - descent gas :math:`\leftrightarrow` dry-ice interphase HTC
     - W/m\ :sup:`2`/K
   * - :math:`A_g,\,A_\text{wet},\,A_{gs}`
     - gas-contact / wetted / gas-solid-interface areas
     - m\ :sup:`2`
   * - :math:`\dot Q_{wg},\,\dot Q_{wl},\,\dot Q_{gs},\,\dot Q_{gl}`
     - wall-gas / wall-liquid / gas-solid / gas-liquid heat rates
     - W
   * - :math:`\dot m_\text{sub}`
     - sublimation mass rate (self-cooling + interphase)
     - kg/s

The ``release:`` input block
============================

The CO\ :sub:`2` add-on is driven by a top-level ``release:`` block (in addition to the
base ``vessel/initial/calculation/valve/heat_transfer`` sections). Its keys:

.. list-table::
   :widths: 26 16 58
   :header-rows: 1

   * - Key
     - Default
     - Meaning
   * - ``type``
     - *required*
     - ``gas`` (vapour-space) or ``liquid`` (bottom) release.
   * - ``diameter``
     - *required*
     - orifice / nozzle throat diameter [m].
   * - ``discharge_coef``
     - 0.62
     - liquid (flashing) discharge coefficient :math:`C_{d,\text{liq}}`.
   * - ``discharge_coef_gas``
     - 0.84
     - gas-tail discharge coefficient :math:`C_{d,\text{gas}}`.
   * - ``liquid_nonequilibrium``
     - 0.0
     - HNE factor :math:`N` (liquid only); 0 = HEM, 1 = frozen.
   * - ``liquid_ne_pressure_scaled``
     - false
     - fade :math:`N` linearly toward the triple point (:eq:`n-of-p` in :ref:`discharge`).
   * - ``liquid_ne_pref``
     - --
     - reference pressure for the :math:`N` fade [Pa].
   * - ``eos``
     - ``CoolProp``
     - thermodynamic backend (CoolProp + solid table).
   * - ``solid_in_vessel``
     - false
     - enable the below-triple in-vessel dry-ice model.
   * - ``solid_h_gas_wall``
     - ``churchill``
     - below-triple gas-wall HTC: ``churchill``, ``calc`` (Geankoplis), or a number.
   * - ``solid_h_inner``
     - ``cooper``
     - plateau wall-liquid boiling HTC: ``cooper``, ``calc`` (Rohsenow), or a number.
   * - ``solid_h_gas_solid``
     - 3.0
     - descent gas :math:`\leftrightarrow` dry-ice interphase HTC [W/m\ :sup:`2`/K].
   * - ``solid_h_gas_liquid``
     - 0.0
     - plateau gas :math:`\leftrightarrow` liquid interphase HTC [W/m\ :sup:`2`/K].
   * - ``solid_gas_wall_frac``
     - --
     - fixed override of the gas-contact wall fraction (default: from geometry).
   * - ``discharge_location``
     - 0.0
     - liquid level [m] at which a liquid draw switches to the gas tail (heel changeover).
   * - ``back_pressure``
     - 101325
     - downstream / choke back pressure [Pa].
   * - ``atm_pressure``
     - 101325
     - atmospheric pressure for the downstream dry-ice flash [Pa].

Related base-model keys used by the add-on: ``initial.gas_temperature`` (superheated gas
zone + split-wall init), ``initial.wall_temperature`` (uniform wall override),
``calculation.non_equilibrium: true``, ``calculation.h_gas_liquid``
(``calc``/``calc_two_sided``/number), and ``heat_transfer.h_outer`` /
``heat_transfer.h_inner``.

.. note::

   For a **full-bore / large-diameter rupture** set ``discharge_coef: 1.0`` **and**
   ``liquid_nonequilibrium: 0.0`` (HEM); see the warning in :ref:`hole-size` on why the
   metastable boost must not be stacked on :math:`C_d = 1`.
