Installation
============

Requirements
------------

HydDown requires Python 3.10, 3.11, or 3.12.

Installation via pip
--------------------

The simplest way to install HydDown is via pip:

.. code-block:: bash

   pip install hyddown

Installation from Source
-------------------------

To install HydDown from source for development:

1. Clone the repository:

.. code-block:: bash

   git clone https://github.com/andr1976/HydDown.git
   cd HydDown

2. Install in development mode:

.. code-block:: bash

   pip install -e .

3. Install dependencies:

.. code-block:: bash

   pip install -r requirements.txt

Dependencies
------------

HydDown depends on the following packages:

- numpy >= 1.19.5
- scipy >= 1.6.0
- matplotlib >= 3.3.3
- pandas >= 1.1.4
- CoolProp (thermodynamic backend)
- PyYAML >= 5.4.1
- Cerberus >= 1.3.3 (input validation)
- fluids (vessel geometry calculations)
- ht (heat transfer correlations)
- tqdm (progress bars)

Optional Dependencies
---------------------

For running the Streamlit web application:

.. code-block:: bash

   pip install streamlit

For testing:

.. code-block:: bash

   pip install pytest pytest-cov
