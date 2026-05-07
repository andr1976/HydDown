# Non-Equilibrium Model (NEM) Development Archive

This directory contains the complete development history, validation data, and diagnostic tools for the HydDown Non-Equilibrium Model (NEM) implementation.

## Directory Structure

```
nem_development/
├── README.md                           (this file)
├── diagnostic_scripts/                 (30 Python scripts)
├── intermediate_docs/                  (27 markdown documents)
├── nem_droste_exp1.yml                 (Droste & Schoen Exp 1 config)
├── nem_droste_exp2.yml                 (Droste & Schoen Exp 2 config)
├── nem_droste_exp3.yml                 (Droste & Schoen Exp 3 config)
├── nem_propane_psv_fire_high_hgl.yml   (Test case config)
├── nem_technical_diagram.png           (NEM technical illustration)
├── droste_exp1_validation.png          (Validation plot - Exp 1)
├── droste_exp2_validation.png          (Validation plot - Exp 2)
├── droste_exp3_validation.png          (Validation plot - Exp 3)
└── droste_exp1_diagnostics.png         (Diagnostic plot - Exp 1)
```

## Overview

The Non-Equilibrium Model (NEM) extends HydDown to handle two-temperature systems where the gas and liquid phases are not in thermal equilibrium. This is critical for fire scenarios where:
- External fire rapidly heats the vessel wall
- Wetted (liquid contact) regions remain cooler than unwetted (gas contact) regions
- Gas and liquid phases have different temperatures
- Phase change (boiling/condensation) occurs at different rates than equilibrium models predict

## Validation

The NEM implementation was validated against the Droste & Schoen (1988) full-scale fire experiments:

**Experiment 1 (Oct 1982):**
- Initial: 2.24 bar, 5°C
- Experimental PSV opening: 214s (16.4 bar)
- Experimental rupture: 720s
- Validation results: Excellent agreement for PSV timing and pressure

**Experiment 2 (Nov 1983):**
- Initial: 13.5 bar, 37°C (preheated)
- Experimental PSV opening: 100s (17.3 bar)
- Experimental rupture: 440s
- Validation results: Good agreement after blowdown correction

**Experiment 3 (Dec 1983):**
- Initial: 2.24 bar, 5°C
- Experimental PSV opening: 220s (16.0 bar)
- Experimental rupture: 540s
- Validation results: Excellent agreement

**Overall NEM Validation Accuracy:**
- PSV Opening Time: 8.4% mean absolute error
- Pressure: 9.4% mean absolute error
- Liquid Temperature: 3.5% mean absolute error

## Key Features

The NEM implementation includes:

1. **Separate Energy Balances:** Independent gas and liquid temperature tracking
2. **Phase Change Modeling:** Mass transfer between phases via boiling/condensation
3. **Heat Transfer Mechanisms:**
   - Gas-to-wall convection (natural/forced)
   - Liquid-to-wall convection (natural/forced)
   - Pool boiling heat transfer
   - Gas-liquid interfacial heat transfer (h_gl correlation)
   - Fire heat loads (Stefan-Boltzmann radiation + convection)
4. **Robust Solver:** Fallback chain (brentq → ridder → bisect) for numerical stability
5. **PSV Cycling:** Accurate modeling of pressure relief valve opening/closing

## Fire Modeling

**Pressure Simulation Fire:**
- Model: `scandpower_pool`
- Unscaled heat flux: 88.0 kW/m²
- Scaling factor: 0.5 (50% engulfment)
- Effective heat flux: 44.0 kW/m²
- Flame temperature: 804°C (1077 K)

**Rupture Analysis Fire:**
- Model: `scandpower_pool_peak`
- Peak heat flux: 131.4 kW/m² (unscaled)
- Flame temperature: 939°C (1212.5 K)
- 49% more severe than pressure simulation

## Diagnostic Scripts (30 files)

The `diagnostic_scripts/` directory contains:
- **Validation scripts:** `run_all_droste_validation.py`, `final_results_comparison.py`
- **Heat transfer diagnostics:** Gas-wall, liquid-wall, gas-liquid interfacial
- **Phase transfer analysis:** Boiling/condensation rate studies
- **Solver testing:** Pressure residual, conservation checks
- **Sensitivity studies:** h_gl correlation, relaxation factors
- **Comparison tools:** NEM vs equilibrium model comparisons

## Intermediate Documentation (27 files)

The `intermediate_docs/` directory tracks the development history:

**Key Documentation:**
- `NEM_IMPLEMENTATION_METHODOLOGY.md` - Complete NEM methodology reference
- `NEM_HEAT_TRANSFER_COMPLETE_FIX.md` - Summary of all heat transfer fixes
- `H_GL_CORRELATION_FINAL.md` - Final gas-liquid HTC correlation
- `DROSTE_VALIDATION_REPORT.md` - Detailed validation results

**Development History:**
- Solver improvements and fixes
- Phase transfer modeling evolution
- Heat transfer correlation development
- Energy balance debugging
- PSV cycling logic fixes

## Technical Illustration

`nem_technical_diagram.png` provides a visual overview of the NEM method showing:
- Vessel cross-section with liquid/gas phases
- Fire engulfment and heat transfer mechanisms
- PSV relief and wall heating
- Phase change phenomena (bubbles, droplets)
- Energy balances and validation results

## Usage

To run NEM validation:
```bash
cd docs/nem_development/diagnostic_scripts
python run_all_droste_validation.py
```

To use NEM in your simulations:
```yaml
calculation:
  type: energybalance
  non_equilibrium: true       # Enable NEM
  h_gas_liquid: 500           # Fixed gas-liquid HTC [W/m²K]
  # or
  # h_gas_liquid: "calc"      # Auto-calculate using natural convection correlation

vessel:
  # ... vessel geometry ...
  liquid_level: 0.4668        # Initial liquid level [m] (required for two-phase)

heat_transfer:
  type: "specified_q"         # or "s-b" for Stefan-Boltzmann fire model
  h_inner: "calc"
  q_outer: 80000              # External heat flux [W/m²] (or time-dependent dict)
```
See `src/hyddown/examples/nem_propane_psv_fire.yml` for a complete example.

## References

**Droste & Schoen (1988):**
Full-scale fire tests on liquefied gas road tankers. International Conference on Safety in Road and Rail Tunnels, Basel, Switzerland.

**NEM Theory:**
- Two-temperature energy balance approach
- Separate gas/liquid phases with interfacial heat transfer
- Phase change via mass transfer (boiling/condensation)
- Non-equilibrium thermodynamics for rapid transients

## Maintainer Notes

**Files modified in main codebase:**
- `src/hyddown/hdclass.py` - NEM solver implementation
- `src/hyddown/transport.py` - `h_gas_liquid_interface()` correlation
- `src/hyddown/validator.py` - NEM parameter validation (`non_equilibrium`, `h_gas_liquid` under `calculation`)
- `src/hyddown/examples/nem_propane_psv_fire.yml` - Example NEM case

**Open Issues:**
- Rupture prediction accuracy: Current simplified model doesn't capture localized hot spots
- Fire severity gap: 44 kW/m² (no rupture) vs 131 kW/m² (early rupture)
- Potential improvement: 1-D transient wall conduction (thermesh.py integration)

---

*Last updated: 2026-03-27*
*NEM implementation completed and validated*
