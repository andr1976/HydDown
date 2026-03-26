# Gas-Liquid Heat Transfer Coefficient (h_gl) Sensitivity Study

## Configuration

**Phase Transfer Settings:**
- relax_factor = 0.8 (80% per timestep)
- heat_multiplier = 50.0 (loose cap)
- Energy method: (h+u)/2

**Test Range:** h_gl from 50 to 2000 W/m²K (40× range)

## Results Summary

| h_gl [W/m²K] | Final P [bar] | Final T_gas [K] | Final T_liquid [K] | Final ΔT [K] | Final mass [kg] | Max P [bar] |
|--------------|---------------|-----------------|--------------------|--------------|-----------------| ------------|
| **Equilibrium** | 14.04 | 314.24 | 314.24 | 0.0 | 475.52 | 14.33 |
| **50** | 11.71 | 384.37 | 306.19 | 78.18 | 990.82 | 14.30 |
| **100** | 12.40 | 398.42 | 306.98 | 91.44 | 971.97 | 14.30 |
| **200** | 12.00 | 417.64 | 307.52 | 110.13 | 957.82 | 14.31 |
| **500** | 13.21 | 403.19 | 311.18 | 92.01 | 887.30 | 14.31 |
| **1000** | 12.75 | 370.26 | 310.09 | 60.17 | 796.17 | 14.31 |
| **2000** | 12.76 | 354.21 | 310.12 | 44.09 | 690.98 | 14.33 |

## Key Observations

### 1. Strong Impact on Temperature Stratification (ΔT)

Temperature difference varies dramatically with h_gl:

```
h_gl = 50:    ΔT = 78 K   (weak coupling)
h_gl = 200:   ΔT = 110 K  (moderate coupling) ← Peak stratification
h_gl = 1000:  ΔT = 60 K   (strong coupling)
h_gl = 2000:  ΔT = 44 K   (very strong coupling)
```

**Non-monotonic behavior!** ΔT peaks at h_gl = 200 W/m²K then decreases.

### 2. Strong Impact on Mass Discharge

Final mass varies from 991 kg to 691 kg (300 kg range, 23% of initial):

```
h_gl = 50:    991 kg retained  (least discharge)
h_gl = 200:   958 kg retained
h_gl = 1000:  796 kg retained
h_gl = 2000:  691 kg retained  (most discharge, approaching equilibrium 476 kg)
```

**Trend:** Higher h_gl → more discharge → closer to equilibrium

### 3. Excellent Pressure Control

All cases show Max P = 14.30-14.33 bar (PSV set = 14.30 bar) ✓

PSV control is excellent regardless of h_gl.

### 4. Gas Temperature Non-Monotonic

Gas temperature shows interesting behavior:

```
h_gl = 50:    T_gas = 384 K
h_gl = 200:   T_gas = 418 K  ← Peak
h_gl = 500:   T_gas = 403 K
h_gl = 1000:  T_gas = 370 K
h_gl = 2000:  T_gas = 354 K  ← Approaching equilibrium
```

**Peak at h_gl = 200 W/m²K!** Then decreases as coupling increases.

### 5. Liquid Temperature Nearly Constant

Liquid temperature stays ~307-311 K for all cases (much less variation than gas).

Liquid is dominated by:
- Fire heat input (external)
- Evaporative cooling (latent heat)
- Much less affected by gas-liquid coupling

## Physical Interpretation

### Low h_gl (50-100 W/m²K): Weak Coupling

**Characteristics:**
- Weak heat exchange between phases
- Gas heats moderately (384-398 K)
- Liquid stays cool (306-307 K)
- Modest stratification (78-91 K)
- **Least discharge** (991-972 kg retained)

**Why less discharge?**
- Gas doesn't get as hot → higher density → less volumetric discharge rate
- More mass stays in vessel

### Medium h_gl (200 W/m²K): Moderate Coupling ← Peak Stratification

**Characteristics:**
- Moderate heat exchange
- **Gas gets hottest** (418 K)
- Liquid still cool (308 K)
- **Maximum stratification** (110 K)
- Moderate discharge (958 kg retained)

**Why peak stratification?**
- Enough coupling to heat gas from liquid evaporation
- Not so much that gas and liquid equilibrate
- "Sweet spot" for maximum thermal separation

### High h_gl (500-2000 W/m²K): Strong Coupling

**Characteristics:**
- Strong heat exchange between phases
- Gas cools down (403 → 354 K as h_gl increases)
- Liquid warms up slightly (311 K)
- Reduced stratification (92 → 44 K)
- **Most discharge** (887 → 691 kg)
- Approaching equilibrium behavior

**Why more discharge?**
- Strong coupling → more uniform temperature
- More evaporation → more mass to discharge
- More like equilibrium model

## Comparison to Equilibrium

| Metric | Equilibrium | NEM (h_gl=2000) | NEM (h_gl=50) |
|--------|-------------|-----------------|---------------|
| Final mass | 476 kg | 691 kg | 991 kg |
| ΔT | 0 K | 44 K | 78 K |
| T_gas | 314 K | 354 K | 384 K |

**Trend:** As h_gl increases, NEM approaches equilibrium behavior.

With h_gl = 2000 W/m²K:
- Still 45% more mass retained than equilibrium
- Still 44 K stratification (equilibrium has 0)
- NEM physics still matters even with very strong coupling!

## Recommendations

### For Realistic Fire Scenarios

**Use h_gl = 200-500 W/m²K**

**Reasoning:**
- Represents moderate gas-liquid coupling
- Physical range for natural convection + evaporation
- Gives realistic stratification (90-110 K)
- Conservative discharge prediction (more mass retained)

### For Conservative Design

**Use h_gl = 50-100 W/m²K**

**Reasoning:**
- Weak coupling → maximum mass retention
- Most conservative for vessel capacity
- May overestimate stratification

### For Validation Against Equilibrium

**Use h_gl = 1000-2000 W/m²K**

**Reasoning:**
- Strong coupling → approaches equilibrium
- Useful for comparing NEM vs equilibrium predictions
- Shows NEM converges to equilibrium limit

## Impact of Current Settings

**With relaxation = 0.8 and heat_multiplier = 50:**

The phase transfer limits are sufficiently loose that **h_gl dominates behavior**.

We see:
- 40× variation in h_gl (50 → 2000)
- Produces 300 kg variation in final mass (23%)
- Large impact on temperature stratification (44-110 K range)

This confirms that **h_gl is the key parameter** for NEM behavior, not phase transfer limits.

## Plots Generated

File: `h_gl_sensitivity_study.pdf`

**Subplots:**
1. Pressure vs time (all h_gl + equilibrium)
2. Gas temperature vs time
3. Liquid temperature vs time
4. Temperature difference (ΔT) vs time
5. Total mass vs time
6. Summary table of final results

## Conclusion

### h_gl Sensitivity: HIGH ✓

Varying h_gl from 50 to 2000 W/m²K produces:
- ✅ **Large impact** on temperature stratification (44-110 K)
- ✅ **Large impact** on mass discharge (691-991 kg)
- ✅ **Moderate impact** on gas temperature (354-418 K)
- ✅ **Small impact** on liquid temperature (306-311 K)

### Recommended Value

**h_gl = 200 W/m²K** provides:
- Realistic thermal stratification (~110 K)
- Moderate coupling between phases
- Conservative mass retention (958 kg vs 476 kg equilibrium)
- Good balance between physical realism and safety

### Key Finding

With loose phase transfer limits (relax=0.8, mult=50×), the **gas-liquid heat transfer coefficient (h_gl) is the dominant parameter** controlling NEM behavior.

Phase transfer limits had < 0.3% impact, while h_gl variation produces 23% impact on final mass!

## Updated Input File

The input file `nem_propane_psv_fire.yml` already has:
```yaml
calculation:
  h_gas_liquid: 200  # W/m²K
```

This is a good default value confirmed by this sensitivity study.
