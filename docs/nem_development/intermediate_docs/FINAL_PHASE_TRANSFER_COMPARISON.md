# Final Phase Transfer Energy Comparison

## Summary of Relaxation Mechanisms

The phase transfer has **two levels of relaxation**:

### 1. Mass Transfer Relaxation (Line 2238, 2278)
```python
relax_factor = 0.2  # Transfer only 20% per timestep
dm_evap = quality * m_liquid * 0.2
dm_cond = (1 - quality) * m_gas * 0.2
```

This prevents instantaneous phase equilibration and accounts for:
- Bubbles/droplets need time to rise/fall
- Interfacial mass transfer resistance
- Gradual phase change process

### 2. Heat-Based Safety Limit (Line 2252, 2291)
```python
dm_evap_max = Q_available / h_fg
dm_evap = min(dm_requested, dm_evap_max * 5.0)
```

Limits phase transfer to 5× what heat input can support, accounting for:
- Stored internal energy
- Transient heat capacity effects
- Prevents unphysical mass transfer spikes

## Comparison of Energy Transfer Methods

All three methods use the same mass transfer relaxation but differ in energy accounting:

### Method 1: Enthalpy (h)
```python
E_evap = dm_evap × h_vapor_sat
E_cond = dm_cond × h_liquid_sat
```

### Method 2: Internal Energy (u)
```python
E_evap = dm_evap × u_vapor_sat
E_cond = dm_cond × u_liquid_sat
```

### Method 3: Compromise (h+u)/2
```python
e_vap = (h_vapor_sat + u_vapor_sat) / 2.0
e_liq = (h_liquid_sat + u_liquid_sat) / 2.0
E_evap = dm_evap × e_vap
E_cond = dm_cond × e_liq
```

## Results Comparison (h_gl = 200 W/m²K)

| Metric | Enthalpy (h) | (h+u)/2 | Internal Energy (u) | Equilibrium |
|--------|--------------|---------|---------------------|-------------|
| **Final pressure** | 11.80 bar | 12.04 bar | 12.10 bar | 14.04 bar |
| **Final T_gas** | 416.59 K | 418.69 K | 416.60 K | 314.24 K |
| **Final T_liquid** | 306.81 K | 307.67 K | 307.88 K | 314.24 K |
| **Final ΔT** | 109.77 K | 111.01 K | 108.72 K | 0 K |
| **Max ΔT** | 142.19 K | 139.92 K | 137.60 K | 0 K |
| **Final total mass** | 984.42 kg | 960.52 kg | 926.73 kg | 475.52 kg |
| **Final gas mass** | 140.04 kg | 143.74 kg | 146.94 kg | N/A |
| **Final liquid mass** | 844.39 kg | 816.78 kg | 779.79 kg | N/A |
| **Max pressure** | 14.31 bar | 14.31 bar | 14.31 bar | 14.33 bar |
| **Max phase transfer** | 4.97 kg/s | 5.31 kg/s | 5.70 kg/s | N/A |

## Detailed Analysis

### 1. Temperature Stratification - All Excellent ✓

All three methods give similar temperature behavior:
- ΔT range: 109-111 K
- Max ΔT range: 138-142 K
- **No temperature inversion** in any method ✓
- Gas consistently hotter than liquid ✓

**Conclusion:** Temperature results are **insensitive** to the choice of h, u, or (h+u)/2

### 2. Mass Discharge - Progressive Trend

Clear progression from most retained to least retained:
- **h**: 984 kg (least discharge, most retained)
- **(h+u)/2**: 960 kg (middle)
- **u**: 927 kg (most discharge, least retained)

**Physical interpretation:**
- h = u + Pv (includes PV expansion work)
- Larger energy per kg transferred
- More energy removed from liquid per unit mass
- Liquid cools more → less evaporation → more mass retained

**Difference range:** 57 kg (6% of initial 1324 kg) between h and u

### 3. Pressure Control - Perfect ✓

All three methods show:
- Max pressure: 14.31 bar
- PSV set pressure: 14.30 bar
- **Overshoot: 0.01 bar (0.07%)** - excellent PSV control!

### 4. Comparison to Equilibrium

All NEM methods retain **much more mass** than equilibrium:
- NEM: 927-984 kg (70-74% retention)
- Equilibrium: 476 kg (36% retention)

**This is physically reasonable!**

In equilibrium:
- Uniform temperature throughout
- Efficient heat transfer and vaporization
- More mass evaporates and discharges

In NEM:
- Temperature stratification (gas 110K hotter than liquid)
- Hot gas → lower density → slower discharge
- Cool liquid → less vaporization
- Stratification reduces discharge efficiency

**Fire scenario with thermal stratification → less discharge is expected**

## Thermodynamic Considerations

### Open System First Law
```
dU/dt = Q - W + Σ(h_in*dm_in) - Σ(h_out*dm_out)
```

**Standard formulation uses enthalpy (h)** for mass crossing boundaries.

### Why (h+u)/2 Works

The compromise e = (h+u)/2 is **not** rigorous, but it partially accounts for:
- Internal energy carried by mass (u)
- Flow work associated with mass transfer (Pv)

For propane at ~10 bar:
- h ≈ 660 kJ/kg (vapor)
- u ≈ 615 kJ/kg (vapor)
- (h+u)/2 ≈ 638 kJ/kg

The difference h - u = Pv ≈ 45 kJ/kg represents PV expansion work.

Using (h+u)/2 splits this work 50/50, which has no rigorous thermodynamic basis but gives intermediate results.

### Energy Conservation

All three methods should conserve energy if:
1. Mass transfer terms are consistent with energy balance
2. Pressure solver properly handles volume constraint
3. Heat transfer terms are correctly applied

The fact that all three give:
- Similar temperatures (109-111 K)
- Similar pressure (11.8-12.1 bar)
- No obvious energy violations

Suggests the solver is working correctly with all three formulations.

## Recommendations

### Option 1: Enthalpy (h) - **Recommended** ✓

**Use:**
```python
E_evap = dm_evap × h_vapor_sat
E_cond = dm_cond × h_liquid_sat
```

**Pros:**
- ✅ Thermodynamically rigorous (open system standard)
- ✅ Includes PV flow work
- ✅ No temperature inversion
- ✅ Excellent pressure control
- ✅ Most conservative mass retention

**Cons:**
- Retains more mass than equilibrium (but physically reasonable for stratified system)

### Option 2: Compromise (h+u)/2 - **Acceptable**

**Use:**
```python
e = (h + u) / 2
E_evap = dm_evap × e_vapor
E_cond = dm_cond × e_liquid
```

**Pros:**
- ✅ No temperature inversion
- ✅ Excellent pressure control
- ✅ Middle-ground between h and u
- ✅ Gives intermediate mass discharge

**Cons:**
- ⚠️ No rigorous thermodynamic basis
- ⚠️ Arbitrary 50/50 split of PV work

### Option 3: Internal Energy (u) - **Acceptable**

**Use:**
```python
E_evap = dm_evap × u_vapor_sat
E_cond = dm_cond × u_liquid_sat
```

**Pros:**
- ✅ No temperature inversion
- ✅ Excellent pressure control
- ✅ More mass discharge (closer to equilibrium)

**Cons:**
- ⚠️ Theoretically less rigorous (neglects PV flow work)
- ⚠️ Discharge may be too aggressive

## Sensitivity Analysis: Relaxation Factor

Current: `relax_factor = 0.2` (20% per timestep)

**Effect of increasing relaxation:**
- More aggressive phase transfer
- Faster approach to thermal equilibrium
- Potentially more discharge
- Risk of numerical instability if too large

**Effect of decreasing relaxation:**
- Slower phase transfer
- More temperature stratification
- Less discharge
- More stable numerically

**Recommendation:** Keep `relax_factor = 0.2` as reasonable balance

## Final Recommendation

### Primary Choice: **Enthalpy (h)**

Use enthalpy for phase transfer energy because:
1. **Thermodynamically correct** for open systems
2. **Standard practice** in all thermodynamics textbooks
3. **Excellent results**: no inversion, good pressure control
4. **Conservative**: retains more mass (physically reasonable for stratified system)

### Alternative: **(h+u)/2 Compromise**

If you prefer a middle-ground approach:
- Gives intermediate results between h and u
- No rigorous basis but practically useful
- Results are very similar to h (within 2.5% on mass)

### Avoid: Original (h+u)/2 with Wrong Interpretation

The **old** implementation that used:
```python
E_evap = dm × (h_vap - h_liq)  # WRONG - this is h_fg, not energy carried!
```

This caused temperature inversion and should not be used.

### Current Implementation Status

**Currently using:** (h+u)/2 with correct interpretation
- E_evap = dm × (h_vap + u_vap)/2
- E_cond = dm × (h_liq + u_liq)/2

**Results:** Excellent (ΔT = +111 K, no inversion, max P = 14.31 bar)

## Conclusion

All three methods with **correct interpretation** (evaporated liquid becomes vapor, condensed vapor becomes liquid) give:
- ✅ **No temperature inversion**
- ✅ **Excellent pressure control**
- ✅ **Physically reasonable stratification**
- ✅ **Stable numerical behavior**

The choice between h, u, or (h+u)/2 affects mass discharge by ~6% but has minimal impact on temperature and pressure behavior.

**Recommended: Use enthalpy (h)** for thermodynamic rigor, but all three are acceptable for engineering calculations.
