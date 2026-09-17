.. _model_map:

===============================
Model Map (Phase Decision Tree)
===============================

A CO\ :sub:`2` blowdown moves through a sequence of thermodynamic regimes, and the
add-on switches the **discharge model, vessel (non-equilibrium) treatment, heat-transfer
law and phase lever** as it crosses each regime boundary. This chapter is the road map:
the decision tree below shows the phases and their transitions, and the table that follows
lists which model is in force in each phase. The individual models are derived in
:ref:`thermodynamics` (backend), :ref:`discharge` (mass flux), :ref:`dry_ice` (phase lever)
and :ref:`heat_transfer` (wall / interphase).

The starting state can be **dense-phase / sub-cooled liquid** (e.g. the 120 bar
Høydalsvik/Munkejord tests) or **saturated two-phase** liquid + vapour (e.g. the
10-20 bar CARDICE tests); both enter the same tree. The first branch is the **release
location** - a vapour-space hole (gas draw) or a liquid-space hole (liquid draw) - which
selects the left or right spine.

.. only:: latex

   .. raw:: latex

      \begin{figure}[htbp]
      \centering
      \begin{tikzpicture}[node distance=9mm and 8mm]
        \node[phase] (start) {\textbf{Start}\\ dense-phase or saturated\\ two-phase CO$_2$ in vessel};
        \node[branch, below=of start] (br) {release\\ location?};
        % --- gas (vapour-space) spine ---
        \node[phase, below left=11mm and 4mm of br] (g1)
             {\textbf{P1} Gas draw\\ above triple (VLE)};
        \node[phase, below=of g1] (g2) {\textbf{P2} Triple-point\\ plateau: liquid $\to$ dry ice};
        \node[phase, below=of g2] (g3) {\textbf{P3} Gas/solid tail\\ sublimation descent};
        % --- liquid (bottom) spine ---
        \node[phase, below right=11mm and 4mm of br] (l1)
             {\textbf{L1} Liquid draw\\ above triple (VLE)};
        \node[phase, below=of l1] (l2) {\textbf{L2} Two-phase drain\\ to heel / changeover};
        \node[phase, below=of l2] (l3) {\textbf{L3} Gas tail\\ (liquid exhausted)};
        % --- common atmospheric lever ---
        \node[phase, fill=green!10] (atm) at ($(g3)!0.5!(l3) + (0,-2.7cm)$)
             {\textbf{Released stream} (every phase)\\ isenthalpic flash to 1 atm\\
              $\to$ vapour + dry-ice split at frost point};
        % edges
        \draw[flowline] (start) -- (br);
        \draw[flowline] (br) -| (g1) node[pos=0.25, above, font=\scriptsize] {vapour space};
        \draw[flowline] (br) -| (l1) node[pos=0.25, above, font=\scriptsize] {liquid space};
        \draw[flowline] (g1) -- (g2); \draw[flowline] (g2) -- (g3);
        \draw[flowline] (l1) -- (l2); \draw[flowline] (l2) -- (l3);
        \draw[flowline] (g3.south) |- (atm.west);
        \draw[flowline] (l3.south) |- (atm.east);
      \end{tikzpicture}
      \caption{CO\textsubscript{2} release phase decision tree. The left spine (P1--P3)
      is a vapour-space (gas) release; the right spine (L1--L3) a liquid-space release.
      The atmospheric dry-ice lever acts on the released stream in every phase.}
      \label{fig:phasemap}
      \end{figure}

.. only:: html

   .. list-table:: Phase decision tree (text form)
      :widths: 50 50
      :header-rows: 1

      * - Gas (vapour-space) release
        - Liquid (bottom) release
      * - **P1** Gas draw, above triple (VLE)
        - **L1** Liquid draw, above triple (VLE)
      * - **P2** Triple-point plateau: liquid :math:`\to` dry ice
        - **L2** Two-phase drain to heel / changeover
      * - **P3** Gas/solid tail: sublimation descent
        - **L3** Gas tail (liquid exhausted)

   In every phase the released stream is separately flashed isenthalpically to 1 atm,
   giving the downstream vapour + dry-ice split at the frost point.

Models employed in each phase
=============================

.. list-table:: Discharge, vessel (NEM), heat transfer and phase lever by phase
   :widths: 8 23 23 23 23
   :header-rows: 1

   * - Phase
     - Discharge
     - Vessel / NEM
     - Heat transfer
     - Phase lever
   * - **P1**
     - HEM, saturated-vapour stagnation, :math:`C_{d,\text{gas}}`
     - CoolProp two-zone NEM (gas + liquid), native EoS
     - two-node wall: gas wall Geankoplis natural convection, wetted wall Rohsenow boiling
     - vapour/liquid equilibrium (VLE)
   * - **P2**
     - HEM gas on the sub-triple throat (forms solid), :math:`C_{d,\text{gas}}`
     - two-zone plateau: warm gas zone + liquid/solid zone at :math:`P_\text{tp}`
     - gas wall Churchill-Chu; **Cooper/Rohsenow boiling** into the wetted liquid
     - three-phase invariant lever :math:`(m_g,m_l,m_s)` from :math:`M,V,U`
   * - **P3**
     - HEM gas on the sublimation line, :math:`C_{d,\text{gas}}`
     - two-zone descent: gas leaks/depressurises, dry ice sublimes
     - gas wall Churchill-Chu; dry-ice wall tracks :math:`T_\text{sub}(P)`; gas/solid over a mass-derived interface
     - sublimation-line lever + :math:`\dot m_\text{sub}` (self-cool + interphase)
   * - **L1**
     - HEM liquid **+ HNE boost** :math:`N`, :math:`C_{d,\text{liq}}`
     - CoolProp two-zone NEM (gas + liquid), native EoS
     - two-node wall: gas wall Geankoplis natural convection, wetted wall Rohsenow boiling
     - VLE; liquid boils off as :math:`P` falls
   * - **L2**
     - HEM liquid + :math:`N(P)` (fades toward the triple point), :math:`C_{d,\text{liq}}`
     - two-zone NEM to the heel; ``discharge_location`` sets the gas changeover
     - two-node wall: gas wall natural convection, wetted wall Rohsenow boiling
     - VLE; near-adiabatic warm-gas / cold-liquid split
   * - **L3**
     - HEM gas, :math:`C_{d,\text{gas}}`
     - residual-gas blowdown (liquid exhausted, negligible in-vessel solid)
     - gas-contact wall convection
     - single-phase gas
   * - all
     - --
     - --
     - --
     - **atmospheric lever**: isenthalpic flash to 1 atm :math:`\to` vapour + dry ice at the frost point

.. note::

   The **triple-point hand-over** (P1 :math:`\to` P2, or L2 :math:`\to` its tail) is where
   the CoolProp above-triple NEM passes the vessel state to the below-triple two-zone model
   (:ref:`thermodynamics`). The **discharge** model does *not* hand over - the HEM rate is
   computed on the same CoolProp + solid-table basis at all pressures, so the mass rate is
   continuous across the triple point.
