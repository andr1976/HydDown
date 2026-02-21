# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

HydDown is a Python package for calculating hydrogen (or other pure gas phase species) pressure vessel filling and discharge incorporating heat transfer effects. It models vessel response (pressure/temperature) to depressurization, pressurization, and external heat loads (e.g., fire scenarios).

**Key capabilities:**
- Single component two-phase modelling with separate gas/liquid heat transfer
- Multiple thermodynamic calculation methods (isothermal, isenthalpic, isentropic, energy balance)
- Various mass flow equations (orifice, control valve, relief valve, constant mass flow)
- Fire scenario modeling using Stefan-Boltzmann approach
- 1-D transient heat conduction for vessel walls

**Python version support:** 3.10 to 3.12

## Building and Testing

### Installation and Setup

Install in development mode:
```bash
pip install -e .
```

Install dependencies:
```bash
pip install -r requirements.txt
```

### Running the Application

Main script execution:
```bash
python scripts/hyddown_main.py input.yml
```

The default input file is `input.yml` in the root directory. Example input files are located in `src/hyddown/examples/`.

### Streamlit Application

Run the interactive web application:
```bash
streamlit run scripts/streamlit_app.py
```

Additional streamlit apps:
- `streamlit_genapp.py` - General purpose calculator
- `streamlit_h2app.py` - H2-specific calculator
- `streamlit_sbapp.py` - Stefan-Boltzmann fire scenario
- `streamlit_bdv_sbapp.py` - Blowdown valve with Stefan-Boltzmann

### Testing

Run all tests:
```bash
cd src/hyddown
pytest
```

Run specific test:
```bash
cd src/hyddown
pytest test_all.py::test_orifice
```

Run with coverage:
```bash
cd src/hyddown
pytest --cov=. --cov-report=xml
```

### Code Quality

Format code with Black:
```bash
black .
```

Lint with flake8:
```bash
flake8 . --count --select=E9,F63,F7,F82 --show-source --statistics --ignore=F821
flake8 . --count --exit-zero --max-complexity=10 --max-line-length=127 --statistics
```

## Architecture

### Core Modules

**`hdclass.py`** (~1900 lines) - Main calculation engine
- `HydDown` class: Central class managing problem definition, calculations, and results
- Implements the time-stepping integration scheme (explicit Euler)
- Handles multiple calculation types: isothermal, isenthalpic, isentropic, energy balance
- Methods for thermodynamic property calculations (PH, UD problems) for both single component and multicomponent fluids
- Key methods:
  - `run()`: Main integration loop for mass/energy balances
  - `step()`: Single time step calculation
  - `PHproblem()`, `UDproblem()`: Thermodynamic state calculations
  - `plot()`: Results visualization

**`validator.py`** (~860 lines) - Input validation
- Uses Cerberus for schema-based validation
- Defines required/optional parameters for all calculation modes
- Validates vessel geometry, initial conditions, valve parameters, heat transfer settings

**`transport.py`** (~590 lines) - Heat and mass transfer calculations
- Dimensionless numbers: Grashof (Gr), Prandtl (Pr), Nusselt (Nu), Rayleigh (Ra)
- Heat transfer coefficient calculations for natural/forced convection
- Mass flow rate calculations for different valve types:
  - `massflow_rate_mdot()`: Constant mass flow
  - `massflow_rate_orifice()`: Orifice equation
  - `massflow_rate_control_valve()`: Control valve sizing equation
  - `massflow_rate_relief_valve()`: Relief valve (API 520/521)
- Boiling heat transfer (pool boiling, film boiling)

**`fire.py`** (~130 lines) - Fire heat load modeling
- Stefan-Boltzmann radiation + convection heat transfer
- Predefined fire scenarios:
  - `pool_fire_api521()`: 60 kW/m² incident heat flux
  - `pool_fire_scandpower()`: 100 kW/m² incident heat flux
  - `jet_fire_api521()`: 100 kW/m² incident heat flux
  - `jet_fire_scandpower()`: 250 kW/m² incident heat flux

**`thermesh.py`** (~430 lines) - 1-D transient heat conduction
- Adapted from https://github.com/wjbg/thermesh
- Finite element method for vessel wall temperature distribution
- Supports composite materials (multi-layer walls)
- Used for Type III/IV vessels with low thermal conductivity

**`materials.py`** (~340 lines) - Material property database
- Thermal properties for vessel materials (steel, aluminum, composites)
- Temperature-dependent properties where applicable

### Data Flow

1. **Input**: YAML file defines vessel geometry, initial conditions, calculation type, valve parameters, heat transfer settings
2. **Validation**: `validator.py` checks input against schema
3. **Initialization**: `HydDown.__init__()` reads input, initializes arrays, sets up thermodynamic backend
4. **Time Integration**: `HydDown.run()` loops through time steps:
   - Calculate mass flow rate (from `transport.py`)
   - Update mass inventory
   - Calculate heat transfer (convection, radiation, fire)
   - Solve thermodynamic state (P, T from H/U and ρ)
   - Update vessel wall temperature (using `thermesh.py` if enabled)
5. **Output**: Results stored in arrays (time, pressure, temperature, mass flow, etc.)
6. **Plotting**: `HydDown.plot()` generates matplotlib figures

### Thermodynamic Backend

HydDown relies heavily on **CoolProp** for fluid property calculations:
- Single component fluids: Uses HEOS (Helmholtz Equation of State) backend
- Multicomponent mixtures: Supported but slower, requires numerical optimization
- Property pairs: P-T, P-H, D-U, D-H, T-S, etc.
- CoolProp syntax: `PropsSI('Property', 'Input1', value1, 'Input2', value2, 'Fluid')`

**Important**: Single component fluids prefixed with `HEOS::` (e.g., `HEOS::Hydrogen`). Multicomponent mixtures use `&` separator (e.g., `HEOS::Methane[0.9]&Ethane[0.1]`).

### Calculation Types

The `calculation.type` parameter determines the thermodynamic path:

1. **isothermal**: Constant temperature (very slow process with large heat reservoir)
2. **isenthalpic**: Constant enthalpy, adiabatic expansion without work
3. **isentropic**: Constant entropy, adiabatic expansion with PV work
4. **specified_U**: Constant internal energy
5. **energybalance**: Most general case, accounts for heat transfer and work

For `energybalance`, the heat transfer type is specified separately:
- `fixed_U`: Fixed U-value (overall heat transfer coefficient)
- `fixed_Q`: Fixed heat input
- `specified_h`: Specified internal/external heat transfer coefficients
- `detailed`: Detailed heat transfer with wall conduction model
- `fire`: Fire heat load from Stefan-Boltzmann equation

### Valve Flow Types

The `valve.type` parameter determines mass flow calculation:

1. **orifice**: Compressible flow through orifice (requires `diameter`, `discharge_coef`)
2. **control_valve**: Control valve sizing equation (requires `Cv`, `N9`)
3. **relief_valve**: API 520/521 relief valve (requires `diameter`, `set_pressure`)
4. **mdot**: Constant mass flow rate (requires `mass_flow`)

Flow direction set by `valve.flow`: `"discharge"` or `"filling"`

## Input File Structure

YAML files define calculations with required sections:

```yaml
vessel:          # Geometry and material properties
initial:         # Starting pressure, temperature, fluid
calculation:     # Type, time step, end time
valve:           # Flow type, size, coefficients
heat_transfer:   # Heat transfer model (if energybalance)
validation:      # Optional validation data for plotting
```

See `src/hyddown/examples/` for complete examples of different calculation types.

## Common Development Tasks

### Adding a New Fluid Property Calculation

Fluid properties are primarily accessed via CoolProp's `PropsSI()` function. For custom calculations, follow patterns in `transport.py` or `hdclass.py`:

```python
from CoolProp.CoolProp import PropsSI
property = PropsSI('PROPERTY_NAME', 'T', T_value, 'P', P_value, species)
```

### Extending Heat Transfer Models

Heat transfer correlations are in `transport.py`. To add new correlations:
1. Define the correlation function (following existing Nu, h patterns)
2. Update `HydDown.step()` in `hdclass.py` to call the new correlation
3. Add corresponding validation schema in `validator.py`

### Adding New Fire Scenarios

Fire scenarios are defined in `fire.py` as functions returning heat flux [W/m²]:
1. Create new function with Stefan-Boltzmann calculation
2. Update `HydDown.read_input()` to recognize the new scenario
3. Add schema validation for new parameters

### Working with Vessel Geometries

Vessel volumes calculated using `fluids.TANK()` class. Supported types:
- `Flat-end`: Simple cylinder
- `ASME F&D`: Torispherical heads (ASME F&D standard)
- `DIN`: Torispherical heads (DIN standard)
- `Semi-elliptical`: Elliptical heads

Orientation: `horizontal` or `vertical` (affects heat transfer correlations)

## Important Implementation Details

### Time Integration Scheme

- Explicit Euler method for mass balance integration
- Time step `dt` must be small enough for stability (typically 0.01-1 second)
- No automatic time step adjustment - user controls via `calculation.time_step`

### Multicomponent vs Single Component

Single component fluids are significantly faster because:
- CoolProp can directly calculate properties from any pair (P-H, D-U, etc.)
- No iterative optimization needed

Multicomponent fluids require:
- Numerical optimization (scipy.optimize.minimize) to find state
- Only T-P pairs directly supported by CoolProp
- Can be very slow for large systems

### Two-Phase Modeling

HydDown supports single component two-phase systems:
- Tracks liquid and gas temperatures separately
- Different wall temperatures for wetted/unwetted regions
- Uses quality (vapor fraction) to determine phase boundaries
- Implements pool boiling and film boiling correlations

### Wall Heat Conduction

Two modes:
1. **Simple**: Uniform wall temperature (lumped capacitance)
2. **Detailed**: 1-D transient conduction via `thermesh.py` (for Type III/IV vessels)

Detailed mode requires:
- Wall material properties (thermal conductivity, heat capacity, density)
- Mesh discretization parameters
- Significantly slower but more accurate for composite materials

## Validation

HydDown has been extensively validated against:
- Published experimental data (see Manual.md)
- External codes (GeoH2, commercial software)
- API 521 relief valve sizing methods
- Literature correlations for heat transfer

Validation data can be included in input YAML files for comparison plotting:
```yaml
validation:
  pressure:
    time: [0, 10, 20, ...]
    pres: [150, 120, 90, ...]
  temperature:
    gas_high:
      time: [...]
      temp: [...]
```

## References and Documentation

- **Manual**: See `Manual.md` or `docs/MANUAL.pdf` for rigorous explanation of methods
- **Citation**: Andreasen, A., (2021). JOSS, 6(66), 3695, https://doi.org/10.21105/joss.03695
- **Streamlit demo**: https://hyddown-jltaqjxtrsflh2famtkgsj.streamlit.app/
- **CoolProp documentation**: http://www.coolprop.org/
