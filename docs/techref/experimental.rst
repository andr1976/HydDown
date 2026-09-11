.. _experimental:

============================
Experimental Basis (CARDICE)
============================

The add-on is validated against the **CARDICE** (CARbon Dioxide ICE) Joint
Industry Project, in which Ineris performed pilot-scale CO\ :sub:`2` blowdown
experiments in a 2 m\ :sup:`3` sphere. Two papers document the work: the GHGT-15
paper :cite:`Vaillant2021`, which gives the test matrix and the Vessfire
comparison, and the IJGGC setup paper :cite:`Jamois2023`, which documents the rig,
instrumentation and calibration in detail. The open 1 Hz dataset is on Zenodo
:cite:`CardiceData2024`.

The vessel
==========

.. figure:: figures/jamois_vessel.png
   :width: 90%

   The 2 m\ :sup:`3` steel sphere: (A) bare, (B) top-view drawing, (C) with flanges
   and insulation, (D) with the aluminium cover.

   From :cite:`Jamois2023`.

The vessel characteristics reported by :cite:`Jamois2023` are summarised in
:numref:`tbl-vessel`.

.. _tbl-vessel:

.. list-table:: Vessel characteristics
   :widths: 45 55
   :header-rows: 1

   * - Property
     - Value
   * - Geometry
     - sphere, two welded hemispheres
   * - Internal volume
     - 2.030 :math:`\pm` 0.005 m\ :sup:`3` (geometric)
   * - Inner radius / diameter
     - 778 mm / 1556 mm
   * - Mean wall thickness
     - 54 :math:`\pm` 5 mm
   * - Total vessel mass
     - 4.80 t (steel shell :math:`\approx` 4 t + 3 flanges)
   * - Steel grade
     - P355GH (1.0473)
   * - Heat capacity :math:`c_p`
     - 450 J/(kg\ :math:`\cdot`\ K)
   * - Thermal conductivity :math:`\lambda`
     - 44 W/(m\ :math:`\cdot`\ K)
   * - Density
     - 7800 kg/m\ :sup:`3`
   * - Design
     - 200 bara / :math:`-60\,^{\circ}`\ C

In HydDown the sphere is modelled as two hemispherical heads with zero cylinder
length (``vessel.type: Hemispherical``, ``length: 0``, ``diameter: 1.58``). The wall
thickness is set to 0.072 m rather than the physical 54 mm so that the
shell-corrected wall mass matches the measured 4.80 t total (the extra accounts for
the three heavy flanges/nozzles, which HydDown does not model separately).

Insulation and external heat transfer
=====================================

The sphere is insulated with 10 cm of rubber foam
(:math:`\lambda_\text{fo} = 0.034` W/(m\ :math:`\cdot`\ K)) over an outer area of
11 m\ :sup:`2`. :cite:`Jamois2023` reports the external heat transfer in two ways:

**Theoretical (steady conduction through the foam):**

.. math::
   :label: hext-theo

   h_\text{theo} = \frac{\lambda_\text{fo}\,A_\text{out}}{L_\text{fo}}
   = \frac{0.034 \times 11}{0.10} = 3.75\ \text{W/K}

**Measured (a cooling calorimetry test, their Fig. 2):** the insulated vessel was
warmed and allowed to cool towards the 2 :math:`^{\circ}`\ C ambient. The contents
dropped :math:`\Delta T = 8\,^{\circ}`\ C over :math:`\Delta t = 45` h. The heat-loss
rate is inferred from the vessel thermal mass :math:`C` (dominated by the 4 t of
steel, :math:`4000 \times 450 = 1.8` MJ/K, giving :math:`\approx 89` W; :math:`\approx 107`
W including the contents):

.. math::
   :label: hext-meas

   \dot{Q} = C\,\frac{\Delta T}{\Delta t} \approx 107\ \text{W},
   \qquad
   h = \frac{\dot{Q}}{\overline{\Delta T}} = \frac{107}{25} = 4.3\ \text{W/K}

where :math:`\overline{\Delta T} \approx 25\,^{\circ}`\ C is the mean inside-outside gap
over the run. Note that this ``h`` is a **whole-vessel conductance in W/K**, not a
per-area coefficient; per unit outer area it is
:math:`4.3 / 11 = 0.39` W/(m\ :sup:`2`\ :math:`\cdot`\ K), which is the value used as
``heat_transfer.h_outer`` in the models. The measured 4.3 W/K exceeds the pure-foam
3.75 W/K by :math:`\approx 15\,\%` because the real heat paths also include the flanges,
the release pipe, and the thermocouple wells (8 mm :math:`\times` 250 mm pipes
protruding through the insulation) - thermal bridges the conduction estimate ignores.

.. figure:: figures/jamois_cooling.png
   :width: 75%

   Cooling curves of the insulated vessel used to obtain the measured external
   heat-transfer conductance.

   From :cite:`Jamois2023`.

For the several-hour CO\ :sub:`2` blowdowns this :math:`\sim 100` W ingress is small
against the latent (flashing/sublimation) loads, so the blowdown is effectively
near-adiabatic.

Instrumentation
===============

.. figure:: figures/vaillant_instrumentation.png
   :width: 85%

   Overall instrumentation of the set-up.

   From :cite:`Vaillant2021`.

The measured quantities are:

* absolute **pressure** inside the sphere and in the release pipe;
* **fluid temperature** at 6 heights inside the sphere (``Tc in 1..6``) and in the pipe;
* inner and outer **wall temperature** at 7 heights (``Tc1..Tc7`` / ``Tc F1..F7``);
* **mass flow rate** through the release pipe;
* the **mass of the vessel** (fluid + solids) via 4 load cells (Mettler-Toledo
  0745A, 2200 kg each, :math:`\pm` 100 g; total-mass uncertainty :math:`\pm` 0.4 t on
  4 t; flow-rate uncertainty :math:`\pm` 5 %).

.. figure:: figures/jamois_heattransfer.png
   :width: 90%

   Example measured temperatures during a blowdown: (A) fluid at 6 heights inside the
   vessel and (B) the wall 5 mm from the inner face, showing the stable vertical
   stratification of the vapour phase.

   From :cite:`Jamois2023`.

Two independent flow measurements were used: direct **weighing** (load cells) and an
**ultrasonic** velocity (inline on the gas pipe, clamp-on on the liquid pipe) combined
with the computed density. The load-cell mass, corroborated against the measured
temperature and pressure, is the reference used here for the initial inventory and
the discharge-rate calibration. A camera and a katharometer were also fitted.

The release pipe is 40 mm inner diameter, :math:`\sim` 2 m long, terminated by a
**calibrated orifice** connected to the vessel middle (vapour release) or bottom
(liquid release). The paper gives one pure-CO\ :sub:`2` orifice example of **4 mm**
(range 1 mm to full bore); per-test sizes are not reported, so they are inferred here
from the measured discharge (:ref:`discharge`).

Test matrix
===========

The pure-CO\ :sub:`2` tests are 5-10; tests 1-4 are methane / methane-CO\ :sub:`2`
mixtures and are out of scope. For the low-temperature tests the sphere was filled
with liquid to about the middle. The matrix is reported by :cite:`Vaillant2021`.

.. list-table:: CARDICE pure-CO\ :sub:`2` test matrix
   :widths: 12 22 18 18 15 15
   :header-rows: 1

   * - Test
     - Composition
     - :math:`P_0` [bar]
     - :math:`T_0` [\ :math:`^{\circ}`\ C]
     - Release
     - HydDown file
   * - 5
     - 100 % CO\ :sub:`2`
     - 20
     - :math:`-20`
     - Gas
     - ``CARDICE_test5.yml``
   * - 6
     - 100 % CO\ :sub:`2`
     - 20
     - :math:`-20`
     - Liquid
     - ``CARDICE_batch6.yml``
   * - 7
     - 100 % CO\ :sub:`2`
     - 15
     - :math:`-29`
     - Gas
     - ``CARDICE_test7.yml``
   * - 8
     - 100 % CO\ :sub:`2`
     - 15
     - :math:`-29`
     - Liquid
     - ``CARDICE_test8.yml``
   * - 9
     - 100 % CO\ :sub:`2`
     - 10
     - :math:`-40`
     - Gas
     - ``CARDICE_test9.yml``
   * - 10
     - 100 % CO\ :sub:`2`
     - 10
     - :math:`-40`
     - Liquid
     - ``CARDICE_test10.yml``

Gas releases produced large in-vessel dry-ice banks (over 450 kg reported), as in
:numref:`fig-dryice`; liquid releases produced essentially none inside the vessel.

.. _fig-dryice:

.. figure:: figures/vaillant_dryice.png
   :width: 60%

   Solid CO\ :sub:`2` formed inside the sphere during a **gas** release test (the
   picture is tilted; the metal rod indicates the vertical).

   From :cite:`Vaillant2021`.

.. figure:: figures/jamois_dryice_foam.png
   :width: 60%

   Close-up of the dry-ice "foam-like" structure (floating solid with vapour
   bubbles).

   From :cite:`Jamois2023`.

Raw-data alignment
==================

The 1 Hz files include a pre-release period; the true blowdown start :math:`t_0` is
detected from the onset of the inventory-mass decline (and the pressure drop) and
all comparisons are time-shifted to it (e.g. :math:`t_0 = 438` s for cardice-06).
Folder ``cardice-0N`` corresponds to test :math:`N`.
