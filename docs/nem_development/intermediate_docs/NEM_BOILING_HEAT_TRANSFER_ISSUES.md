# NEM Nucleate Boiling Heat Transfer Issues

## Summary

The Rohsenow nucleate boiling correlation for wetted wall heat transfer has three critical bugs when used with the Non-Equilibrium Model (NEM):

1. **Wrong phase check**: Uses equilibrium fluid quality instead of liquid phase quality
2. **Wrong saturated properties**: Uses equilibrium fluid object instead of liquid phase object
3. **Wrong temperature**: Uses bulk equilibrium temperature instead of liquid temperature

These bugs cause the wetted heat transfer coefficient to incorrectly drop to zero when NEM develops thermal stratification.

## Issue #1: Wrong Phase Check

**Location:** `hdclass.py` lines 1246, 1267, 1281

**Current Code:**
```python
if self.fluid.Q() >= 0 and self.fluid.Q() <= 1:
    hiw = tp.h_inside_wetted(...)
else:
    hiw = hi  # Uses gas-side coefficient!
```

**Problem:**
- `self.fluid` is the EQUILIBRIUM fluid object
- For NEM with stratification, `self.fluid` might be superheated gas (Q < 0 or Q > 1)
- But the LIQUID PHASE still needs nucleate boiling heat transfer!
- Result: After t~350s, equilibrium fluid becomes superheated → hiw = hi (gas coefficient) → WRONG!

**Evidence from wetted_htc_comparison.py:**
```
t=300s:  h_inside_wetted = 3000 W/m²K  (correct - nucleate boiling)
t=350s:  h_inside_wetted = 0 W/m²K     (WRONG - switched to gas coefficient)
t=400s:  h_inside_wetted = 0 W/m²K     (WRONG - switched to gas coefficient)
```

**Fix:**
```python
# For NEM, check liquid phase quality
if self.non_equilibrium:
    liquid_is_boiling = (self.m_liquid[i-1] > 1e-6 and
                        self.fluid_liquid.Q() >= 0 and
                        self.fluid_liquid.Q() <= 1)
else:
    liquid_is_boiling = self.fluid.Q() >= 0 and self.fluid.Q() <= 1

if liquid_is_boiling:
    hiw = tp.h_inside_wetted(...)
else:
    hiw = hi
```

## Issue #2: Wrong Saturated Properties

**Location:** `transport.py` lines 348-364, called from `hdclass.py` lines 1250-1256, 1282-1288

**Current Code:**
```python
# In hdclass.py
hiw = tp.h_inside_wetted(
    L, T_wall, T_fluid,
    self.transport_fluid_wet,  # fluid parameter
    self.fluid)                # master_fluid parameter - EQUILIBRIUM!

# In transport.py
def h_inside_wetted(L, Tvessel, Tfluid, fluid, master_fluid):
    # Saturated properties from master_fluid
    rhol=master_fluid.saturated_liquid_keyed_output(CP.iDmass),
    rhog=master_fluid.saturated_vapor_keyed_output(CP.iDmass),
    Hvap=(master_fluid.saturated_vapor_keyed_output(CP.iHmass)
          - master_fluid.saturated_liquid_keyed_output(CP.iHmass)),
    sigma=master_fluid.surface_tension(),
```

**Problem:**
- `self.fluid` is the EQUILIBRIUM object
- For NEM, saturated properties should come from `self.fluid_liquid` at liquid phase conditions
- Equilibrium object might be at different state than liquid phase

**Current behavior:**
- Both equilibrium and NEM liquid phase properties happen to be identical at same pressure
- But conceptually wrong - should use liquid phase object

**Fix:**
```python
# For NEM, use liquid phase object for saturated properties
if self.non_equilibrium:
    master_fluid_for_boiling = self.fluid_liquid
else:
    master_fluid_for_boiling = self.fluid

hiw = tp.h_inside_wetted(
    L, T_wall, T_fluid,
    self.transport_fluid_wet,
    master_fluid_for_boiling)  # Correct phase object
```

## Issue #3: Wrong Temperature

**Location:** `hdclass.py` lines 1252-1254, 1284-1286

**Current Code:**
```python
hiw = tp.h_inside_wetted(
    L,
    self.T_inner_wall_wetted[i - 1],  # Wall temperature - OK
    self.T_fluid[i - 1],               # WRONG for NEM - equilibrium T!
    self.transport_fluid_wet,
    self.fluid,
)
```

**Problem:**
- Uses `self.T_fluid[i-1]` which is the equilibrium bulk temperature
- For NEM, should use `self.T_liquid[i-1]` for liquid phase temperature
- This affects:
  - Film temperature for transport properties: `T_film = (T_wall + T_fluid) / 2`
  - Wall superheat in Rohsenow: `Te = T_wall - T_fluid`

**Impact:**
- With equilibrium T_fluid = 348 K (gas temperature)
- Should use T_liquid = 305 K (liquid temperature)
- ΔT = 43 K error in film temperature!
- Wall superheat Te also calculated incorrectly

**Fix:**
```python
if self.non_equilibrium:
    T_liquid_for_boiling = self.T_liquid[i - 1]
else:
    T_liquid_for_boiling = self.T_fluid[i - 1]

hiw = tp.h_inside_wetted(
    L,
    self.T_inner_wall_wetted[i - 1],
    T_liquid_for_boiling,  # Correct liquid temperature
    self.transport_fluid_wet,
    master_fluid_for_boiling,
)
```

## Impact on Results

### Equilibrium vs NEM with h_gl=500 (from wetted_htc_comparison.py)

**Average values (t > 300s):**

| Parameter | Equilibrium | NEM (BUGGY) | Difference |
|-----------|-------------|-------------|------------|
| T_liquid [°C] | 36.1 | 33.5 | -2.6 K |
| T_wall_wetted [°C] | 64.0 | **174.7** | **+110.7 K** |
| ΔT (wall-liquid) [K] | 27.9 | **141.1** | **+113.2 K** |
| h_inside_wetted [W/m²K] | 3000.0 | **241.2** | **-2758.8** |
| q_inner_wetted [kW/m²] | 83.6 | **2.9** | **-80.7** |

**Result:** Wall temperature skyrockets to 175°C because boiling heat transfer shuts off!

### Physical Explanation

1. **Early times (t < 350s):**
   - Equilibrium fluid is two-phase → `self.fluid.Q()` ∈ [0,1] → nucleate boiling active
   - h_inside_wetted = 3000 W/m²K (capped maximum)
   - Heat transfer works correctly

2. **Late times (t > 350s):**
   - NEM develops stratification → gas phase superheats
   - Equilibrium fluid becomes superheated gas → `self.fluid.Q()` < 0 or > 1
   - Code incorrectly thinks boiling stopped → sets `hiw = hi` (gas coefficient ≈ 0)
   - But liquid phase is STILL THERE and STILL BOILING!
   - Result: Wall temperature rises to 254°C with no heat removal

## Recommended Fixes

### Fix 1: NEM-specific logic in hdclass.py

Add NEM detection and use liquid phase properties:

```python
# Around line 1246-1258
if self.heat_method == "specified_h" or self.heat_method == "specified_q":
    T_film = (self.T_fluid[i - 1] + self.T_inner_wall[i - 1]) / 2

    # NEM: Use liquid phase properties for wetted heat transfer
    if self.non_equilibrium:
        # Check if liquid phase exists and is boiling
        liquid_exists = self.m_liquid[i-1] > 1e-6
        if liquid_exists:
            # Update liquid fluid to check phase
            self.fluid_liquid.update(CP.DmassUmass_INPUTS,
                                    self.rho_liquid[i-1],
                                    self.U_liquid[i-1])
            liquid_is_boiling = (self.fluid_liquid.Q() >= 0 and
                                self.fluid_liquid.Q() <= 1)
        else:
            liquid_is_boiling = False

        if liquid_is_boiling:
            # Use liquid temperature and liquid fluid object
            T_film_wet = (self.T_liquid[i-1] + self.T_inner_wall_wetted[i-1]) / 2
            self.transport_fluid_wet.update(CP.PT_INPUTS, self.P[i-1], T_film_wet)

            hiw = tp.h_inside_wetted(
                L,
                self.T_inner_wall_wetted[i-1],
                self.T_liquid[i-1],           # Liquid temperature
                self.transport_fluid_wet,
                self.fluid_liquid,            # Liquid phase object
            )
        else:
            hiw = hi  # No liquid or not boiling
    else:
        # Equilibrium: Use existing logic
        if self.fluid.Q() >= 0 and self.fluid.Q() <= 1:
            self.transport_fluid_wet.update(CP.PT_INPUTS, self.P[i-1], T_film)
            hiw = tp.h_inside_wetted(
                L,
                self.T_inner_wall_wetted[i-1],
                self.T_fluid[i-1],
                self.transport_fluid_wet,
                self.fluid,
            )
        else:
            hiw = hi
```

### Fix 2: Similar changes needed around line 1267-1290

The same NEM-specific logic needs to be added in the other branch where `hiw` is calculated.

## Testing

After fixes, expect:
- h_inside_wetted to remain active (non-zero) in NEM when liquid is present
- Wall temperature to stay reasonable (not skyrocket to 254°C)
- Heat transfer from wall to liquid to continue throughout simulation
- Results to show proper nucleate boiling behavior

## Files to Modify

1. **`src/hyddown/hdclass.py`**
   - Lines 1246-1258: Add NEM logic for first hiw calculation
   - Lines 1267-1290: Add NEM logic for second hiw calculation

2. **Optional: `src/hyddown/transport.py`**
   - Could add parameter documentation clarifying that `master_fluid` should be the liquid phase object for NEM

## Related Analysis Files

- `nem_boiling_properties_diagnostic.py`: Diagnostic showing property differences
- `wetted_htc_comparison.py`: Demonstrates the bug impact on results
- `wetted_htc_equilibrium_vs_nem_h500.pdf`: Visualization of incorrect behavior
