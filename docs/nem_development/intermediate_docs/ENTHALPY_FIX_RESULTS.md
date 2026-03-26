# Phase Transfer Energy Fix: Using Enthalpy Instead of (h+u)/2

## Change Made

**Previous implementation** (lines 2258-2266):
```python
# Energy transferred with phase change (using h+u)/2 compromise
e_liq = (h_liq_sat + u_liq_sat) / 2.0
e_vap = (h_vap_sat + u_vap_sat) / 2.0
E_evap = dm_evap_mass * (e_vap - e_liq)
```

**New implementation**:
```python
# Energy transferred with phase change (using enthalpy)
# When mass evaporates:
#   - Liquid loses: dm * h_liquid (enthalpy leaves with mass)
#   - Gas gains: dm * h_vapor (enthalpy arrives with mass)
h_liq_sat = self.fluid_liquid.saturated_liquid_keyed_output(CP.iHmass)
h_vap_sat = self.fluid_liquid.saturated_vapor_keyed_output(CP.iHmass)
E_evap = dm_evap_mass * (h_vap_sat - h_liq_sat)  # Net energy = dm * h_fg
```

## Results Comparison (h_gl = 200 W/m²K)

| Metric | (h+u)/2 | Enthalpy (h) | Status |
|--------|---------|--------------|---------|
| **Final ΔT** | -11.27 K | **+172.43 K** | ✅ Fixed inversion |
| **Temperature order** | Liquid > Gas | **Gas > Liquid** | ✅ Physical |
| **Max ΔT** | 88.29 K | 172.43 K | ✅ More stratification |
| **Final mass** | 162.48 kg | 211.06 kg | ✓ More realistic |
| **Max pressure** | 16.44 bar | 16.95 bar | Similar |
| **Phase transfer rate** | 20.08 kg/s | 24.27 kg/s | Slightly higher |

## Physical Interpretation

### Why Enthalpy is Correct

When mass transfers between phases, it carries **enthalpy** with it, not internal energy:

**For evaporation (liquid → gas):**
1. Mass dm leaves liquid phase at saturation: carries away h_liq × dm
2. Mass dm enters gas phase at saturation: brings in h_vap × dm
3. Net energy transfer: E = dm × (h_vap - h_liq) = dm × h_fg

**For condensation (gas → liquid):**
1. Mass dm leaves gas phase: carries away h_vap × dm
2. Mass dm enters liquid phase: brings in h_liq × dm
3. Net energy transfer: E = dm × (h_liq - h_vap) = -dm × h_fg

The enthalpy h = u + P/ρ naturally includes the PV work done during phase change.

### Why (h+u)/2 Was Wrong

The compromise e = (h+u)/2 attempted to account for the fact that:
- Mass transfer involves enthalpy (h)
- Energy balance uses internal energy (U)

But this is thermodynamically inconsistent! It led to:
- **Temperature inversion**: Liquid became hotter than gas (unphysical)
- **Poor energy balance**: Phases didn't couple correctly
- **Incorrect energy flows**: Neither enthalpy nor internal energy conservation

### Why It Works Now

Using enthalpy for mass transfer is correct because:

1. **Mass carries enthalpy**: When fluid mass moves, it carries h (not u)
2. **Energy balance closes**: The difference h - u = P/ρ is the PV work, which is internally consistent in a closed vessel
3. **Thermodynamically rigorous**: Standard approach in thermodynamics textbooks

## Remaining Issues

While the temperature inversion is fixed, NEM still differs from equilibrium:

| Metric | NEM (h_gl=200) | Equilibrium | Difference |
|--------|----------------|-------------|------------|
| Final mass | 211 kg | 476 kg | **2.3× more discharge** |
| Final pressure | 12.52 bar | 14.04 bar | 1.5 bar lower |
| Max pressure | 16.95 bar | 14.33 bar | **2.6 bar higher** |

**Possible causes:**
1. **Gas superheat** (533K) leads to lower density → higher discharge velocity
2. **Pressure spikes** above PSV set pressure allow more mass to escape
3. **Phase transfer rate limits** may still be too permissive

## Next Steps

1. ✅ **Fixed**: Temperature inversion eliminated by using enthalpy
2. ⚠️ **Investigate**: Why NEM discharges 2× more mass than equilibrium
3. ⚠️ **Tune**: h_gl value (currently 200 W/m²K may still be too high/low)
4. ⚠️ **Check**: Phase transfer rate safety limits (5× multiplier may need adjustment)

## Conclusion

**Major improvement!** Using enthalpy instead of (h+u)/2 for phase transfer energy:
- ✅ Eliminates unphysical temperature inversion
- ✅ Maintains correct temperature stratification (gas > liquid)
- ✅ More thermodynamically rigorous
- ✅ Simulation is stable

The NEM model now has correct physics for phase transfer energy accounting.
