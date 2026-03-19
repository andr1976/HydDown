hdclass - Main Calculation Engine
==================================

The ``hdclass`` module contains the core ``HydDown`` class that manages the calculation engine,
time-stepping integration, and thermodynamic state calculations.

.. automodule:: hyddown.hdclass
   :members:
   :undoc-members:
   :show-inheritance:
   :special-members: __init__

Key Classes
-----------

HydDown
~~~~~~~

.. autoclass:: hyddown.hdclass.HydDown
   :members:
   :undoc-members:
   :show-inheritance:
   :special-members: __init__

Key Methods
~~~~~~~~~~~

The ``HydDown`` class implements the following key methods:

- ``run()``: Main integration loop for mass/energy balances
- ``step()``: Single time step calculation
- ``PHproblem()``: Solve thermodynamic state from pressure and enthalpy
- ``UDproblem()``: Solve thermodynamic state from internal energy and density
- ``plot()``: Generate results visualization

Implementation Notes
--------------------

Time Integration
~~~~~~~~~~~~~~~~

The code uses an explicit Euler method for mass balance integration. The time step ``dt``
must be small enough for numerical stability (typically 0.01-1 second).

Calculation Types
~~~~~~~~~~~~~~~~~

The ``calculation.type`` parameter determines the thermodynamic path:

- ``isothermal``: Constant temperature
- ``isenthalpic``: Constant enthalpy, adiabatic expansion without work
- ``isentropic``: Constant entropy, adiabatic expansion with PV work
- ``specified_U``: Constant internal energy
- ``energybalance``: Most general case with heat transfer
