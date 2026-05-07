# Phase Transfer Energy: Enthalpy (h) vs Internal Energy (u)

## Corrected Formulations

Both formulations now use the correct interpretation:

### Using Enthalpy (h)
```python
# Evaporated liquid becomes saturated vapor
E_evap = dm_evap × h_vapor_sat  (from liquid parent phase)

# Condensed vapor becomes saturated liquid
E_cond = dm_cond × h_liquid_sat  (from gas parent phase)
```

### Using Internal Energy (u)
```python
# Evaporated liquid becomes saturated vapor
E_evap = dm_evap × u_vapor_sat  (from liquid parent phase)

# Condensed vapor becomes saturated liquid
E_cond = dm_cond × u_liquid_sat  (from gas parent phase)
```

## Energy Balance Application

Both are applied identically:
```python
# Gas phase:
U_gas_new = U_gas_old + ... + E_evap - E_cond

# Liquid phase:
U_liquid_new = U_liquid_old + ... - E_evap + E_cond
```

## Results Comparison (h_gl = 200 W/m²K)

| Metric | Enthalpy (h) | Internal Energy (u) | Equilibrium | Notes |
|--------|--------------|---------------------|-------------|-------|
| **Final pressure** | 11.80 bar | 12.10 bar | 14.04 bar | u closer to equilibrium |
| **Final T_gas** | 416.59 K | 416.60 K | 314.24 K | Nearly identical |
| **Final T_liquid** | 306.81 K | 307.88 K | 314.24 K | Nearly identical |
| **Final ΔT** | 109.77 K | 108.72 K | 0 K | Nearly identical |
| **Max ΔT** | 142.19 K | 137.60 K | 0 K | Nearly identical |
| **Final total mass** | 984.42 kg | 926.73 kg | 475.52 kg | h retains more mass |
| **Final gas mass** | 140.04 kg | 146.94 kg | N/A | Similar |
| **Final liquid mass** | 844.39 kg | 779.79 kg | N/A | h retains more liquid |
| **Max pressure** | 14.31 bar | 14.31 bar | 14.33 bar | Both match PSV perfectly ✓ |

## Key Observations

### 1. Temperature Behavior - Nearly Identical ✓
Both h and u give almost the same temperature stratification:
- ΔT ≈ 109 K (no inversion!)
- Max ΔT ≈ 138-142 K
- Gas stays hotter than liquid throughout simulation

### 2. Mass Discharge - Significant Difference
- **With h**: Final mass = 984 kg (retains 74% of initial mass)
- **With u**: Final mass = 927 kg (retains 70% of initial mass)
- **Equilibrium**: Final mass = 476 kg (retains 36% of initial mass)

Both NEM cases discharge much less than equilibrium!

### 3. Pressure Control - Excellent ✓
Both h and u have:
- Max pressure = 14.31 bar
- PSV set pressure = 14.30 bar
- **Perfect control** (only 0.01 bar overshoot!)

### 4. Physical Interpretation

**Why does NEM retain more mass than equilibrium?**

In equilibrium model:
- Single uniform temperature throughout
- All fluid at same T ≈ 314 K
- Efficient discharge at consistent properties

In NEM:
- Gas superheat: T_gas ≈ 417 K >> T_liquid ≈ 308 K
- Hot gas → lower density → slower mass discharge
- Liquid stays cool → less vaporization → more mass retained
- Temperature stratification reduces effective discharge rate

**This is physically reasonable for a fire scenario with thermal stratification!**

## Difference Between h and u

### Theoretical Difference

At saturation:
```
h = u + Pv = u + P/ρ
```

For propane at ~10 bar:
- h_fg ≈ 360 kJ/kg (latent enthalpy)
- u_fg ≈ 335 kJ/kg (latent internal energy)
- Difference ≈ 25 kJ/kg (7% of h_fg)

### Why h gives more retained mass

**With enthalpy:**
- E_evap = dm × h_vapor (larger energy transfer)
- More energy removed from liquid per kg evaporated
- Liquid cools more → less further evaporation
- Result: More liquid remains

**With internal energy:**
- E_evap = dm × u_vapor (smaller energy transfer)
- Less energy removed from liquid per kg evaporated
- Liquid stays warmer → more evaporation
- Result: More liquid evaporates (less remains)

The difference: h includes PV expansion work, u does not.

## Which is Correct?

### Thermodynamic First Law for Open Systems

**Standard formulation:**
```
dU/dt = Q - W + Σ(h_in × dm_in/dt) - Σ(h_out × dm_out/dt)
```

Where:
- U = internal energy (what we're solving for)
- h = enthalpy (what mass carries when crossing boundaries)

**Mass carries enthalpy, not internal energy!**

When fluid mass crosses a control volume boundary (like phase interface), it carries:
- Its internal energy (u)
- Plus flow work (Pv)
- **Total = enthalpy (h = u + Pv)**

### For Closed System (No Mass Transfer)

For a closed system with internal mass redistribution:
```
dU/dt = Q - W
```

But NEM has **mass transfer between control volumes** (liquid phase and gas phase are separate control volumes).

### Recommendation: Use Enthalpy (h) ✓

**Reasons:**
1. **Thermodynamically rigorous**: Open system first law requires enthalpy for mass transfer
2. **Standard practice**: All thermodynamics textbooks use h for flow across boundaries
3. **Includes flow work**: h = u + Pv accounts for expansion work during phase change
4. **Results are stable**: No temperature inversion, good pressure control

**Using internal energy (u) is an approximation** that:
- Neglects PV flow work
- Happens to give similar temperature results (109 K vs 108 K)
- But gives different mass discharge behavior

## Impact on Results

### Temperature Stratification: Minimal Difference
- Both h and u give ΔT ≈ 109 K
- Temperature evolution nearly identical
- Both avoid unphysical inversion ✓

### Mass Discharge: ~6% Difference
- h: 984 kg final (26% discharged)
- u: 927 kg final (30% discharged)
- Difference: 57 kg (6% of initial mass)

For engineering purposes, this is **acceptable accuracy**.

### Energy Conservation

Both methods should conserve energy if implemented correctly. The difference is where the PV work appears:
- **With h**: PV work implicit in enthalpy
- **With u**: PV work must be explicit in pressure-volume term

Since we're solving for U (internal energy), using h requires that the PV work is properly accounted for in the pressure solver.

## Conclusion

### Recommended: Use Enthalpy (h)

**Implementation:**
```python
# Evaporation
u_vap_sat = self.fluid_liquid.saturated_vapor_keyed_output(CP.iUmass)
h_vap_sat = self.fluid_liquid.saturated_vapor_keyed_output(CP.iHmass)
E_evap = dm_evap × h_vap_sat  # Use h

# Condensation
u_liq_sat = self.fluid_gas.saturated_liquid_keyed_output(CP.iUmass)
h_liq_sat = self.fluid_gas.saturated_liquid_keyed_output(CP.iHmass)
E_cond = dm_cond × h_liq_sat  # Use h
```

**Justification:**
- ✅ Thermodynamically correct for open system
- ✅ Standard practice in thermodynamics
- ✅ Includes PV flow work
- ✅ Excellent temperature behavior
- ✅ Perfect pressure control
- ✅ Physically reasonable mass discharge

### Alternative: Internal Energy (u) is Acceptable

If using u for simplicity:
- Temperature results nearly identical
- Slightly more mass discharge (~6%)
- Theoretically less rigorous but practically acceptable

**Both methods are vastly better than the original (h+u)/2 compromise!**

## Summary Table

| Aspect | (h+u)/2 Original | Internal Energy (u) | Enthalpy (h) ✓ |
|--------|------------------|---------------------|----------------|
| **Temperature inversion** | ❌ Yes (-11K) | ✅ No (+109K) | ✅ No (+110K) |
| **Thermodynamic rigor** | ❌ Inconsistent | ⚠️ Approximate | ✅ Correct |
| **Final mass** | 162 kg | 927 kg | 984 kg |
| **Pressure control** | Good | Excellent | Excellent |
| **Recommendation** | ❌ Don't use | ⚠️ Acceptable | ✅ **Recommended** |
