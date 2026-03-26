# NEM Implementation Status

## Changes Made

### 1. Removed Fallback Code ✅
Per user request: "fail hard so we can make the solver robust and stable"

**Removed:**
- Fallback at line ~2350: h_valve calculation fallback (now raises exception if fails)
- Fallback at line ~2457-2546: Pressure equilibrium violation fallback (now raises RuntimeError)

**Result:** Solver now strictly enforces pressure equilibrium. Any solver failure will stop execution with detailed error message.

### 2. Added Phase Transfer Safety Limits ✅
**Problem:** Massive phase transfer spikes (343 kg/s) causing unphysical pressure spikes

**Solution:** Limit phase transfer rate based on available heat:
```python
dm_evap_max = Q_liquid / h_fg
dm_evap = min(dm_evap_requested, dm_evap_max * 5.0)
```

**Result:** Phase transfer rate reduced from 343 kg/s to 20 kg/s (17× reduction)

### 3. Fixed Timing References ✅
Use `Q_inner[i-1]` and `Q_inner_wetted[i-1]` (previous timestep) instead of `[i]` (current timestep not yet calculated)

## Current Results

### Comparison: Equilibrium vs NEM (Propane PSV Fire)

| Metric | Equilibrium | NEM | Issue |
|--------|-------------|-----|-------|
| **Final P** | 14.04 bar | 11.89 bar | ⚠️ Lower |
| **Final T** | 314 K | Gas: 593K, Liquid: 347K | ❌ Extreme superheat |
| **Final mass** | 476 kg | 217 kg | ❌ 2.2× more discharge |
| **Max P** | 14.33 bar | 16.19 bar | ⚠️ 1.9 bar overshoot |
| **Phase transfer** | N/A | Max: 20 kg/s, Avg: 1.7 kg/s | ✅ Reasonable |

### Specific Issues

#### Issue 1: Extreme Gas Superheat
- Final gas temperature: **593 K** (320°C!)
- This is physically unrealistic for propane in fire scenario
- Suggests insufficient evaporative cooling

#### Issue 2: Temperature Inversion
- At t=464s (max pressure): T_liquid = 337K **>** T_gas = 321K
- This is **unphysical** for a stratified system
- Hot gas should float above cool liquid, not vice versa

#### Issue 3: Excessive Mass Discharge
- NEM discharges: 1323 - 217 = **1106 kg**
- Equilibrium discharges: 1323 - 476 = **847 kg**
- NEM discharges 30% more mass

#### Issue 4: Liquid Level Geometry
- Final liquid mass: 215 kg
- Final liquid level: **0.000 m** (contradiction!)
- Suggests geometry calculation issue

## Possible Root Causes

### 1. Gas-Liquid Heat Transfer (Q_gas_liquid)
The gas-liquid heat transfer term may not be correctly coupling the two phases. If this is too weak:
- Gas heats up without transferring heat to liquid
- Liquid stays cool or heats independently
- Results in temperature inversion

### 2. Phase Transfer Logic
Currently checking if **previous timestep** was in two-phase region:
```python
self.fluid_liquid.update(CP.DmassUmass_INPUTS, self.rho_liquid[i-1], self.U_liquid[i-1])
quality_liquid = self.fluid_liquid.Q()
```

**Problem:** This is reactive, not predictive. By the time we detect two-phase, the phase has already moved far from saturation.

**Alternative:** Check if phase is deviating from saturation temperature:
```python
T_sat = get_T_sat(P_current)
if T_liquid > T_sat + tolerance:
    # Force evaporation based on temperature deviation
```

### 3. Energy Transfer with Phase Change
Currently using `(h + u)/2` compromise:
```python
e_liq = (h_liq_sat + u_liq_sat) / 2.0
e_vap = (h_vap_sat + u_vap_sat) / 2.0
E_evap = dm_evap * (e_vap - e_liq)
```

This may not be thermodynamically consistent. For evaporation:
- Mass leaves liquid at h_liq → takes enthalpy
- Mass enters gas at h_vap → brings enthalpy
- Net effect: subtract h_fg from liquid, add h_fg to gas

But internal energy changes are:
- Liquid: ΔU = Δ(mU) - PΔV work
- Gas: ΔU = Δ(mU) + PΔV work

The (h+u)/2 compromise doesn't properly account for PV work.

### 4. Pressure Solver Bounds
The pressure solver uses:
```python
P_min = 0.5 * self.P[i-1]
P_max = 2.0 * self.P[i-1]
```

During rapid depressurization or pressure spikes, this may be too restrictive or too loose.

## Recommendations

### Option 1: Investigate Gas-Liquid Heat Transfer
**Priority: HIGH**

Check the Q_gas_liquid calculation:
1. Is the heat transfer coefficient correct?
2. Is the interface area properly calculated?
3. Is the temperature difference (T_gas - T_liquid) used correctly?

**Expected:** Q_gas_liquid should flow from hot phase to cold phase and equilibrate temperatures over time.

### Option 2: Saturation Temperature Enforcement
**Priority: MEDIUM**

Instead of checking quality at previous state, enforce saturation limits:
```python
P_current = self.P[i-1]  # or guess
T_sat = get_T_sat_from_P(P_current)

# Liquid cannot be much hotter than saturation
if T_liquid > T_sat + 2K:
    dm_evap = rate_function(T_liquid - T_sat)

# Gas cannot be much cooler than saturation
if T_gas < T_sat - 2K:
    dm_cond = rate_function(T_sat - T_gas)
```

### Option 3: Fix Energy Transfer Thermodynamics
**Priority: MEDIUM**

Use proper enthalpy for mass transfer:
```python
# Mass leaving liquid takes enthalpy with it
E_out_liquid = dm_evap * h_liquid_current

# Mass entering gas brings enthalpy
E_in_gas = dm_evap * h_vapor_sat

# Apply to energy balance:
U_liquid_new = U_liquid - E_out_liquid
U_gas_new = U_gas + E_in_gas
```

But need to be careful about PV work terms.

### Option 4: Reduce Timestep
**Priority: LOW (expensive)**

Current timestep: 0.5s may be too large for NEM coupling dynamics.
Try: 0.1s or 0.05s

Trade-off: 5-10× longer simulation time

## Next Steps

1. **Investigate Q_gas_liquid calculation** - ensure it's physically correct
2. **Add diagnostic output** - print T_gas, T_liquid, Q_gas_liquid at each timestep
3. **Test with smaller timestep** - verify if dt=0.1s improves results
4. **Consider saturation temperature approach** - may be more robust than quality-based

## Test Cases

### Quick Test
```bash
python analyze_nem_pressure.py
```
Shows: PSV opening time, max pressure, phase transfer rates

### Full Comparison
```bash
python compare_nem.py
```
Shows: Equilibrium vs NEM side-by-side with plots

### Energy Conservation
```bash
python debug_energy_detailed.py
```
Shows: Timestep-by-timestep energy balance
