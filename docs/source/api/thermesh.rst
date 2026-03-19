thermesh - 1-D Transient Heat Conduction
==========================================

The ``thermesh`` module provides 1-D transient heat conduction modeling for vessel walls
using the finite element method. This is particularly important for Type III/IV vessels
with low thermal conductivity composite materials.

.. automodule:: hyddown.thermesh
   :members:
   :undoc-members:
   :show-inheritance:

Overview
--------

The module is adapted from https://github.com/wjbg/thermesh and implements finite element
analysis for transient heat conduction through vessel walls.

Key Features
------------

- Finite element method for accurate temperature distribution
- Support for composite materials (multi-layer walls)
- Time-dependent boundary conditions
- Integration with fire heat load calculations
- Separate models for wetted/unwetted regions in two-phase flow

Applications
------------

The detailed heat conduction model is essential for:

- Type III/IV composite pressure vessels
- Materials with low thermal conductivity
- Accurate wall temperature prediction during fire scenarios
- Two-phase systems with different heat transfer in gas/liquid contact regions

Wall Heat Conduction Modes
---------------------------

HydDown supports two modes:

1. **Simple**: Uniform wall temperature (lumped capacitance) - fast but approximate
2. **Detailed**: 1-D transient conduction via ``thermesh`` - slower but accurate

The detailed mode is activated by setting appropriate parameters in the ``heat_transfer``
section of the input file.

Composite Materials
-------------------

For multi-layer walls (e.g., Type III vessels with liner, composite overwrap, and outer layer):

1. Define material properties for each layer
2. Specify layer thicknesses
3. The solver handles interfaces automatically

Two-Phase Modeling
------------------

For two-phase systems, the module solves separate heat conduction problems for:

- Wetted region (liquid-contact) - typically higher heat transfer
- Unwetted region (gas-contact) - typically lower heat transfer

The boundary between regions is determined by the liquid level in the vessel.
