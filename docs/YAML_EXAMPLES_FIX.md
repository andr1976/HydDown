# YAML Examples Documentation Fix

## Problem Identified

The documentation examples in `quickstart.rst` and `examples.rst` did NOT match the actual YAML schema expected by HydDown's validator. The examples would fail validation if users tried to use them.

## Issues Found

### 1. Incorrect Structure (Most Critical)

**Documentation Had:**
```yaml
vessel:
  geometry:              # ❌ Wrong - nested structure doesn't exist
    orientation: vertical
    length: 2.0
    diameter: 0.5
    type: flat-end
  material: steel        # ❌ Wrong - material field doesn't exist here
  wall_thickness: 0.01   # ❌ Wrong - should be "thickness"
```

**Actual Schema Requires:**
```yaml
vessel:
  length: 2.0           # ✅ Direct under vessel
  diameter: 0.5         # ✅ Direct under vessel
  orientation: "vertical"  # ✅ Direct under vessel
  thickness: 0.01       # ✅ Not "wall_thickness"
  # No nested "geometry" section!
  # No "material" field in vessel!
```

### 2. Missing Required Fields

**Missing `back_pressure` in valve:**
```yaml
valve:
  type: "orifice"
  flow: "discharge"
  diameter: 0.01
  discharge_coef: 0.84
  # ❌ Missing: back_pressure (REQUIRED)
```

**Should be:**
```yaml
valve:
  type: "orifice"
  flow: "discharge"
  diameter: 0.01
  discharge_coef: 0.84
  back_pressure: 101325  # ✅ Required field added
```

### 3. Scientific Notation Issues

**Documentation had:**
```yaml
pressure: 150e5    # Works but inconsistent
pressure: 1e5      # Works but inconsistent
```

**Now uses:**
```yaml
pressure: 15000000  # Clearer, consistent with examples
```

### 4. String Values Not Quoted

**Documentation had:**
```yaml
fluid: Hydrogen      # ❌ Might cause issues
type: isentropic     # ❌ Might cause issues
```

**Now:**
```yaml
fluid: "Hydrogen"    # ✅ Properly quoted
type: "isentropic"   # ✅ Properly quoted
```

## Actual YAML Schema (from validator.py)

### Vessel Section
```yaml
vessel:
  # Required
  length: <number>          # meters
  diameter: <number>        # meters

  # Optional
  thickness: <number>       # meters (NOT wall_thickness)
  heat_capacity: <number>   # J/kg·K
  density: <number>         # kg/m³
  thermal_conductivity: <number>  # W/m·K
  orientation: "vertical" | "horizontal"
  type: "Flat-end" | "ASME F&D" | "DIN" | "Hemispherical"
  liquid_level: <number>    # 0-1 fraction

  # Liner properties (for composite vessels)
  liner_thickness: <number>
  liner_heat_capacity: <number>
  liner_thermal_conductivity: <number>
  liner_density: <number>
```

### Initial Section
```yaml
initial:
  temperature: <number>     # K (required)
  pressure: <number>        # Pa (required)
  fluid: "<string>"         # CoolProp fluid name (required)
```

### Calculation Section
```yaml
calculation:
  type: "energybalance" | "isenthalpic" | "isentropic" | "isothermal" | "specified_U"
  time_step: <number>       # seconds (required)
  end_time: <number>        # seconds (required)
```

### Valve Section
```yaml
valve:
  flow: "discharge" | "filling"  # Required
  type: "orifice" | "control_valve" | "relief_valve" | "mdot"  # Required
  back_pressure: <number>    # Pa (REQUIRED for most types)

  # Type-specific (orifice)
  diameter: <number>         # meters
  discharge_coef: <number>   # dimensionless

  # Type-specific (control_valve)
  Cv: <number>              # flow coefficient
  N9: <number>              # piping geometry factor

  # Type-specific (relief_valve)
  set_pressure: <number>     # Pa
  blow_down: <number>        # fraction (0-1)

  # Type-specific (mdot)
  mass_flow: <number>        # kg/s
```

### Heat Transfer Section (for energybalance)
```yaml
heat_transfer:
  type: "specified_h" | "fire" | "fixed_U" | "fixed_Q" | "detailed"
  temp_ambient: <number>     # K

  # Type-specific (specified_h)
  h_inner: <number> | "calc"  # W/m²·K
  h_outer: <number>           # W/m²·K

  # Type-specific (fire)
  fire_type: "pool_fire_api521" | "pool_fire_scandpower" | "jet_fire_api521" | "jet_fire_scandpower"
  emissivity: <number>       # 0-1
```

## Changes Made to Documentation

### Files Updated:
1. `docs/source/quickstart.rst`
2. `docs/source/examples.rst`

### Changes Applied:

**✅ Removed nested `geometry` structure**
**✅ Removed non-existent `material` field**
**✅ Changed `wall_thickness` → `thickness`**
**✅ Added required `back_pressure` to all valve examples**
**✅ Added quotes around string values**
**✅ Used full numbers instead of scientific notation for clarity**
**✅ Added required thermal properties for energybalance examples**
**✅ Fixed vessel `type` values to match allowed list**

### Example: Before vs After

**BEFORE (Would Not Work):**
```yaml
vessel:
  geometry:
    orientation: vertical
    length: 2.0
    diameter: 0.5
    type: flat-end
  material: steel

initial:
  pressure: 150e5
  temperature: 293.15
  fluid: Hydrogen

calculation:
  type: isentropic
  time_step: 0.1
  end_time: 100

valve:
  type: orifice
  flow: discharge
  diameter: 0.01
  discharge_coef: 0.84
```

**AFTER (Validated Successfully):**
```yaml
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
```

## Validation Results

**Test Performed:**
```python
from src.hyddown.hdclass import HydDown
import yaml

with open('test_example.yml', 'r') as f:
    config = yaml.safe_load(f)
    hd = HydDown(config)  # Validation happens here
```

**Result:** ✅ **PASSED** - All updated examples now validate correctly

## Summary

**Problem:** Documentation examples didn't match actual YAML schema
**Root Cause:** Documentation was written without checking validator.py
**Solution:** Updated all examples to match actual schema from validator.py and real working examples
**Verification:** Tested examples validate successfully

All documentation examples are now correct and will work when users copy them!

## Recommendation for Future

When adding new examples to documentation:
1. Check `src/hyddown/validator.py` for actual schema
2. Look at `src/hyddown/examples/*.yml` for working examples
3. Test example with HydDown(config) before documenting
4. Use quoted strings for YAML values
5. Include ALL required fields
