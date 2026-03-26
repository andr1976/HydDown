# Task: Validate NEM Implementation Against Droste & Schoen Fire Experiments

## Objective

Validate the Non-Equilibrium Model (NEM) implementation in HydDown against experimental fire test data from Droste & Schoen (1988). Compare simulation results for time to PSV activation, time to rupture, pressure evolution, and liquid temperature for three uninsulated propane tank experiments.

## Background

### NEM Implementation Status
The NEM (Non-Equilibrium Model) has been fully implemented in HydDown with the following features:
- Separate gas and liquid phase tracking at different temperatures
- Phase-specific heat transfer (gas-wall, liquid-wall via nucleate boiling)
- Automatic gas-liquid interfacial heat transfer coefficient calculation using natural convection correlations
- Integration with rupture analysis for fire scenarios

**Key implementation files:**
- `src/hyddown/hdclass.py` - NEM solver (lines 2300-2700)
- `src/hyddown/transport.py` - h_gas_liquid_interface() correlation (lines 375-520)
- `src/hyddown/validator.py` - Input validation for NEM parameters

**Core documentation:**
- `NEM_IMPLEMENTATION_METHODOLOGY.md` - Complete methodology
- `NEM_HEAT_TRANSFER_COMPLETE_FIX.md` - Heat transfer implementation details
- `H_GL_CORRELATION_FINAL.md` - Gas-liquid HTC correlation

**Example NEM input file:**
- `src/hyddown/examples/nem_propane_psv_fire.yml`

### Experimental Data Sources

**Primary reference:**
- File: `github/blowdown-paper/background/Droste_Schoen_JHM_20_1988 (1).pdf`
- Original experimental study by Droste & Schoen (1988)

**Secondary reference (simulations by others):**
- File: `github/blowdown-paper/background/237175021.pdf`
- Previous simulation attempts of the same experiments

**Key experimental data to extract:**
- **Figure 2(b)**: Pressure vs time for experiments 1, 2, and 3
- **Figure 2(a)**: Temperature vs time (believed to be liquid temperature)

## Target Experiments

Focus on **Experiments 1, 2, and 3** (uninsulated tanks):
- All use propane as working fluid
- Different tank sizes/geometries
- Same PSV size for all three experiments
- Similar fire conditions

## Validation Metrics

1. **Time to PSV activation** - When pressure reaches PSV set pressure
2. **Time to rupture** - Using `analyse_rupture()` method in HydDown
3. **Pressure evolution** - Compare P(t) curve to Fig 2(b)
4. **Liquid temperature evolution** - Compare T_liquid(t) to Fig 2(a)

## Known Parameters and Uncertainties

### PSV Parameters (UNCERTAIN)
- Reported as "1 inch valve" (this is likely the flange connection size, not orifice)
- **Best guess:** API 526 orifice designation "D"
  - API 526 "D" orifice: diameter = 11.1 mm = 0.0111 m
  - May need to treat as adjustable parameter if validation fails
- **Key:** PSV is the same size for all three experiments

### Fire Parameters (UNCERTAIN)
- **Fire type options:**
  - `fire: "api_pool"` - API 521 pool fire (60 kW/m²)
  - `fire: "scandpower_pool"` - Scandpower pool fire (100 kW/m²)
- **Scaling parameter:**
  - May need `scaling: <1.0` in heat_transfer section to reduce fire intensity
  - This accounts for partial tank engulfment or reduced heat flux
  - Scaling multiplies the incident heat flux
- **Constraint:** Fire type should be similar for all three experiments

### Tank Geometry and Materials
- Extract from experimental paper (Droste & Schoen):
  - Vessel length, diameter, orientation
  - Wall thickness, material (likely carbon steel)
  - Initial fill level (liquid_level as fraction of diameter)
  - Head type (hemispherical, ASME F&D, etc.)

### Initial Conditions
- Extract from experimental paper:
  - Initial pressure
  - Initial temperature
  - Fluid: propane

### PSV Settings
- Extract from experimental paper:
  - Set pressure
  - Blowdown (if reported, otherwise typical 0.2-0.3)
  - Back pressure (likely atmospheric: 101325 Pa)
  - Discharge coefficient (typical 0.975 for PSV)

### Material Properties
- Vessel material: likely CS (carbon steel)
  - Density: 7700 kg/m³
  - Heat capacity: 500 J/kg·K
  - Check if temperature-dependent rupture material specified

## Recommended Approach

### Step 1: Extract Experimental Data
1. Read both PDF files carefully
2. Extract tank geometry for experiments 1, 2, 3
3. Extract initial conditions (P, T, fill level)
4. Digitize pressure curves from Fig 2(b) if possible
5. Digitize temperature curves from Fig 2(a) if possible
6. Note any reported PSV activation times
7. Note any reported rupture times

### Step 2: Create Input Files
Create three YAML files based on `src/hyddown/examples/nem_propane_psv_fire.yml`:
- `nem_droste_exp1.yml`
- `nem_droste_exp2.yml`
- `nem_droste_exp3.yml`

**Template structure:**
```yaml
vessel:
  length: <from paper>
  diameter: <from paper>
  orientation: "horizontal"  # or vertical
  type: "Hemispherical"  # or ASME F&D
  heat_capacity: 500
  density: 7700
  thickness: <from paper>
  liquid_level: <from paper>

initial:
  temperature: <from paper>
  pressure: <from paper>
  fluid: "propane"

calculation:
  type: "energybalance"
  time_step: 0.5
  end_time: <until rupture + margin>
  non_equilibrium: true
  h_gas_liquid: "calc"  # Use automatic correlation

valve:
  flow: "discharge"
  type: "psv"
  diameter: 0.0111  # API 526 "D" - ADJUSTABLE
  discharge_coef: 0.975
  set_pressure: <from paper>
  blowdown: 0.25  # Typical value if not specified
  back_pressure: 101325.0

heat_transfer:
  type: "s-b"  # Stefan-Boltzmann fire
  fire: "api_pool"  # or scandpower_pool - ADJUSTABLE
  scaling: 1.0  # ADJUSTABLE if needed
  h_inner: "calc"

rupture:
  material: "CS_360LT"  # or appropriate material
  fire: "api_pool_peak"  # Should match heat_transfer fire type
```

### Step 3: Run Simulations with Parameter Tuning

**Fixed parameters (same for all 3 experiments):**
- PSV diameter (find value that works for all 3)
- Fire type (find type that works for all 3)
- Fire scaling (find scaling that works for all 3)

**Tuning strategy:**
1. Start with PSV diameter = 0.0111 m (API 526 "D")
2. Start with fire = "api_pool", scaling = 1.0
3. Run all three simulations
4. Compare PSV activation times to experiments
5. Compare pressure evolution to Fig 2(b)
6. Compare liquid temperature to Fig 2(a)
7. Adjust PSV diameter if activation times are systematically off
8. Adjust fire type/scaling if pressure/temperature evolution is off
9. **Important:** Use same PSV and fire parameters for all three

### Step 4: Analyze Rupture

For each experiment, use:
```python
hd = HydDown("nem_droste_exp1.yml")
hd.run()

# Analyze rupture
rupture_result = hd.analyse_rupture()
print(f"Rupture time: {rupture_result['rupture_time']} s")
print(f"Rupture pressure: {rupture_result['rupture_pressure']} bar")
print(f"Rupture temperature: {rupture_result['rupture_temperature']} K")
```

Compare predicted rupture times to experimental values.

### Step 5: Validation Plots

Create comparison plots for each experiment:
```python
import matplotlib.pyplot as plt

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))

# Pressure comparison
ax1.plot(hd.time_array, hd.P/1e5, 'b-', label='NEM Simulation')
ax1.plot(exp_time, exp_pressure, 'ro', label='Droste & Schoen Exp')
ax1.set_ylabel('Pressure [bar]')
ax1.set_xlabel('Time [s]')
ax1.legend()
ax1.grid(True)

# Temperature comparison
ax2.plot(hd.time_array, hd.T_liquid - 273.15, 'b-', label='NEM T_liquid')
ax2.plot(exp_time, exp_temp, 'ro', label='Droste & Schoen Exp')
ax2.set_ylabel('Temperature [°C]')
ax2.set_xlabel('Time [s]')
ax2.legend()
ax2.grid(True)

plt.tight_layout()
plt.savefig(f'validation_exp{i}.pdf')
```

## Key Files to Review

**Before starting, read these:**
1. `github/blowdown-paper/background/Droste_Schoen_JHM_20_1988 (1).pdf` - Extract all experimental details
2. `github/blowdown-paper/background/237175021.pdf` - See how others simulated these experiments
3. `NEM_IMPLEMENTATION_METHODOLOGY.md` - Understand NEM solver
4. `H_GL_CORRELATION_FINAL.md` - Understand automatic h_gl calculation

**Code to understand:**
1. `src/hyddown/hdclass.py` - Main NEM solver, analyse_rupture() method
2. `src/hyddown/transport.py` - h_gas_liquid_interface() correlation
3. `src/hyddown/fire.py` - Fire scenario definitions (api_pool, scandpower_pool, etc.)

## Expected Deliverables

1. **Three validated YAML input files** for experiments 1, 2, 3
2. **Validation plots** comparing simulations to experimental data
3. **Summary table** of validation metrics:
   - Time to PSV activation (sim vs exp)
   - Time to rupture (sim vs exp)
   - Final parameters used (PSV diameter, fire type, scaling)
4. **Discussion** of NEM model performance and any discrepancies

## Notes

- The NEM model with automatic h_gl correlation is ready for production use
- Expect good agreement if experimental parameters are correctly identified
- PSV diameter and fire scaling are the main tunable parameters
- All three experiments should use the same PSV and fire settings
- Rupture analysis is critical for fire scenario validation

## Questions to Answer During Validation

1. Does the automatic h_gl correlation (h_gas_liquid: "calc") provide realistic results?
2. What PSV orifice size best matches all three experiments?
3. What fire type and scaling best matches the experimental conditions?
4. How well does NEM predict liquid stratification compared to equilibrium model?
5. Are rupture times accurately predicted?
