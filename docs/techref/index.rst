.. _index:

============================================================
HydDown CO\ :sub:`2` Release --- Technical Reference
============================================================

.. only:: latex

   .. raw:: latex

      \vspace{1.5cm}
      \begin{center}
      {\Large Physical models, thermodynamics, heat and mass transfer, and validation\\
      of the thermopack-based CO\textsubscript{2} release / dry-ice add-on to HydDown.}
      \end{center}
      \vspace{1cm}

      \begin{tabular}{ll}
      \textbf{Version:}  & 1.0 \\
      \textbf{Date:}     & September 2025 \\
      \textbf{Author:}   & Anders Andreasen \\
      \textbf{Branch:}   & \texttt{co2-release-hem} \\
      \end{tabular}

      \newpage

This technical reference documents the CO\ :sub:`2` **release and dry-ice** add-on
developed on the ``co2-release-hem`` branch of HydDown :cite:`Andreasen2021`. It
covers the additional physical models, the thermopack thermodynamic backend and its
hand-overs with CoolProp, the homogeneous-equilibrium (HEM) and non-equilibrium
(HNE) discharge models, the in-vessel and atmospheric dry-ice estimation, the
heat- and mass-transfer treatment, and the validation against the CARDICE Joint
Industry Project experiments (Ineris 2 m\ :sup:`3` sphere).

It is a companion to the base HydDown *Manual* (``Manual.md``); only the CO\ :sub:`2`
extension is described here. The document follows the structure of the ThermoCalc
technical reference.

.. toctree::
   :maxdepth: 2
   :caption: Contents

   introduction
   experimental
   thermodynamics
   discharge
   dry_ice
   heat_transfer
   validation

.. only:: latex

   Bibliography
   ============

.. bibliography::
   :all:
