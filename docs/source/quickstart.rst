Quick Start
===========

Basic Usage
-----------

Running HydDown is as simple as:

.. code-block:: bash

   python scripts/hyddown_main.py input.yml

where ``input.yml`` is your input file in YAML syntax.

Simple Example
--------------

Here's a minimal example for a vessel depressurization calculation:

.. code-block:: yaml

   vessel:
     length: 2.0
     diameter: 0.5
     orientation: "vertical"

   initial:
     pressure: 15000000
     temperature: 293.15
     fluid: "Hydrogen"

   calculation:
     type: "isentropic"
     time_step: 0.1
     end_time: 100

   valve:
     flow: "discharge"
     type: "orifice"
     diameter: 0.01
     discharge_coef: 0.84
     back_pressure: 101325

Save this as ``input.yml`` and run:

.. code-block:: bash

   python scripts/hyddown_main.py input.yml

Streamlit Application
---------------------

HydDown includes an interactive web interface built with Streamlit:

.. code-block:: bash

   streamlit run scripts/streamlit_app.py

This provides a user-friendly interface for:

- Setting up calculations without writing YAML
- Visualizing results in real-time
- Comparing with validation data

Available Streamlit Apps
~~~~~~~~~~~~~~~~~~~~~~~~

- ``streamlit_app.py`` - Main application
- ``streamlit_genapp.py`` - General purpose calculator
- ``streamlit_h2app.py`` - Hydrogen-specific calculator
- ``streamlit_sbapp.py`` - Stefan-Boltzmann fire scenario
- ``streamlit_bdv_sbapp.py`` - Blowdown valve with fire

Example Files
-------------

Example input files for various calculation types are located in:

.. code-block:: bash

   src/hyddown/examples/

These include:

- Isothermal filling/discharge
- Isentropic expansion
- Isenthalpic throttling
- Energy balance with heat transfer
- Fire scenarios
- Two-phase calculations
- Relief valve sizing
