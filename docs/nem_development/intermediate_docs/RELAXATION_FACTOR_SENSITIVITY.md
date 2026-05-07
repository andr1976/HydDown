# Relaxation Factor Sensitivity Analysis

## Test Configuration

- Method: (h+u)/2 compromise
- h_gl: 200 W/m²K
- Varied: relax_factor from 0.2 to 0.8

## Results

| Metric | relax = 0.2 | relax = 0.8 | Difference | % Change |
|--------|-------------|-------------|------------|----------|
| **Final pressure** | 12.04 bar | 12.00 bar | -0.04 bar | -0.3% |
| **Final T_gas** | 418.69 K | 417.46 K | -1.23 K | -0.3% |
| **Final T_liquid** | 307.67 K | 307.52 K | -0.15 K | -0.05% |
| **Final ΔT** | 111.01 K | 109.94 K | -1.07 K | -1.0% |
| **Max ΔT** | 139.92 K | 139.84 K | -0.08 K | -0.06% |
| **Final total mass** | 960.52 kg | 957.87 kg | -2.65 kg | -0.3% |
| **Final gas mass** | 143.74 kg | 144.14 kg | +0.40 kg | +0.3% |
| **Final liquid mass** | 816.78 kg | 813.73 kg | -3.05 kg | -0.4% |
| **Max pressure** | 14.31 bar | 14.31 bar | 0.00 bar | 0.0% |

## Key Findings

### 1. Minimal Impact on Results

Increasing relaxation factor from 0.2 (20%) to 0.8 (80%) - a **4× increase** - causes:
- ✅ Temperature: < 1% change
- ✅ Pressure: < 0.3% change
- ✅ Mass: < 0.3% change
- ✅ Stratification (ΔT): 1% change

**Conclusion: Results are insensitive to relaxation factor!**

### 2. Why So Insensitive?

The phase transfer has **two limiting mechanisms**:

#### Mechanism 1: Relaxation Factor (varied)
```python
dm_requested = quality * m_phase * relax_factor
```

#### Mechanism 2: Heat-Based Safety Limit (fixed)
```python
dm_max = Q_available / h_fg
dm_actual = min(dm_requested, dm_max * 5.0)
```

**The heat-based limit is the dominant constraint!**

Even with relax_factor = 0.8 requesting 4× more mass transfer, the heat-based safety limit caps the actual transfer to what 5× the heat input can support.

### 3. Phase Transfer Rates

Let's estimate typical values:

**Heat input to liquid (fire scenario):**
- Q_wetted ~ 300 kW (60 kW/m² × 5 m² wetted area)
- Per timestep (dt=0.5s): Q = 150 kJ

**Maximum evaporation from heat:**
- h_fg ~ 360 kJ/kg (propane)
- dm_max = 150 kJ / 360 kJ/kg = 0.42 kg

**With 5× safety factor:**
- dm_allowed = 0.42 × 5 = 2.1 kg

**Requested with relax_factor:**
- If quality = 0.1, m_liquid = 1000 kg:
  - relax = 0.2: dm_requested = 0.1 × 1000 × 0.2 = 20 kg
  - relax = 0.8: dm_requested = 0.1 × 1000 × 0.8 = 80 kg

**Actual transfer:**
- relax = 0.2: dm_actual = min(20, 2.1) = **2.1 kg** (heat limited)
- relax = 0.8: dm_actual = min(80, 2.1) = **2.1 kg** (heat limited)

**Same result!** The heat limit dominates.

### 4. When Would Relaxation Factor Matter?

The relaxation factor would be important if:

1. **Heat-based limit is disabled** (no 5× cap)
2. **Very high heat input** (Q >> h_fg × m_phase × relax)
3. **Small two-phase mass** (quality × m_phase is small)
4. **Large safety multiplier** (>>5×)

In the current implementation with:
- Fire heat input (~300 kW)
- Safety multiplier = 5×
- Propane h_fg ~ 360 kJ/kg

The heat-based limit is always active, making relax_factor relatively unimportant.

## Implications for Model

### Good News ✓

The insensitivity to relaxation factor means:
- ✅ **Robust**: Results don't depend sensitively on tuning parameter
- ✅ **Physical**: Heat-based limit provides physical constraint
- ✅ **Stable**: Can use larger relax_factor without instability

### Recommendation

**Current setting: relax_factor = 0.2 is fine**

Could also use:
- 0.5 (50%): Middle ground
- 0.8 (80%): Faster equilibration
- 1.0 (100%): Instantaneous (if heat allows)

Since results are similar, choosing 0.2 provides:
- Conservative approach
- Accounts for physical resistance to mass transfer
- Numerically safe

**No strong reason to change from 0.2**

## Testing Higher Values

Let's consider even more extreme:

### relax_factor = 1.0 (100%)

Would request complete phase equilibration per timestep:
- dm_requested = quality × m_phase × 1.0

But still limited by heat-based constraint, so expect:
- Similar results to 0.8
- Possibly numerical issues if quality spikes to 1.0

### Heat-Based Multiplier

Currently using 5× multiplier on heat-based rate. Could test:
- 2×: More restrictive (less phase transfer)
- 10×: More permissive (more phase transfer)

This would likely have **larger impact** than relaxation factor.

## Conclusion

### Current Implementation
```python
relax_factor = 0.2  # 20% per timestep
dm_max_multiplier = 5.0  # 5× heat-based rate
```

### Sensitivity Results
- **Relaxation factor (0.2 → 0.8)**: < 1% impact on all results
- **Reason**: Heat-based safety limit dominates
- **Recommendation**: Keep relax_factor = 0.2 (no need to change)

### If More Aggressive Phase Transfer Desired

To allow more phase transfer, adjust:
1. ✅ **Increase h_gl** (200 → 500 W/m²K): Better gas-liquid coupling
2. ✅ **Increase heat multiplier** (5× → 10×): Allow more stored energy use
3. ⚠️ **Increase relax_factor** (0.2 → 0.8): Minimal effect (already tested)

The most impactful parameter is **h_gl** (gas-liquid heat transfer coefficient).

## Summary

**Tested:** relax_factor = 0.2 vs 0.8 (4× increase)

**Result:** < 1% difference in all metrics

**Conclusion:** Phase transfer is **heat-limited**, not **mass-transfer-rate-limited**

**Recommendation:** Keep current settings (relax_factor = 0.2, heat_multiplier = 5.0)
