# Phase Transfer Comparison: Pre vs Post vs Pre+Post

## Summary Table

| Approach | Final T_gas | Final T_liquid | ΔT | Final m_liquid | Energy Residuals | Phase Transfer Rate |
|----------|-------------|----------------|-----|----------------|------------------|---------------------|
| **Pre-solver only** | 385 K | 346 K | 46 K | 22 kg | 5-13% | ~0.02 kg/step |
| **Post-solver only** | 662 K | 306 K | 356 K | 1161 kg | 8-10% | ~0.001 kg/step |
| **Pre + Post** | 614 K | 306 K | 308 K | 1095 kg | ~17% | ~0.02 kg/step |

## Analysis

### Pre-Solver Only (Original)
✅ **Pros:**
- Moderate temperature stratification (46K)
- Good energy conservation (5-13%)
- Reasonable phase transfer rate
- Most liquid evaporates (22kg remaining)

❌ **Cons:**
- Checks phase change at old density (inconsistent)
- May over-predict evaporation

### Post-Solver Only
✅ **Pros:**
- Thermodynamically consistent state check
- Good energy conservation (8-10%)

❌ **Cons:**
- Extreme superheat (662K gas!)
- Minimal evaporation (1161kg liquid remains)
- Unrealistic stratification (356K)
- Once at equilibrium, phases don't transfer

### Pre + Post Combined
✅ **Pros:**
- Captures immediate heating response (pre)
- Corrects inconsistencies (post)
- Moderate phase transfer rate

❌ **Cons:**
- Still high superheat (614K)
- Poor energy conservation (17%)
- Post-solver barely triggers (pre dominates)

## Physical Interpretation

The core issue: **All methods struggle with superheat accumulation**

### Why Gas Gets Hot

1. Gas receives heat from wall: Q_gas = h*A*(T_wall - T_gas)
2. As gas heats up, it wants to increase pressure
3. Pressure solver enforces P_gas = P_liquid
4. To maintain same P at higher U, gas density must decrease OR temperature must increase
5. With limited evaporation, gas temperature rises dramatically

### The Missing Physics

Real vessels have:
1. **Continuous equilibration**: Vapor pressure continuously adjusts
2. **Interfacial mass transfer**: Not just at saturation, but driven by T gradients
3. **Convection**: Gas circulation mixes hot/cold regions
4. **Smaller effective timesteps**: Nature doesn't use 0.5s steps!

Our explicit Euler with 0.5s timesteps can't capture this.

## Recommendations

### Option 1: Smaller Timesteps ⭐ RECOMMENDED
Current: dt = 0.5s
Try: dt = 0.1s or 0.05s

**Rationale:**
- More frequent pressure equilibration
- More frequent phase transfer checks
- Better captures continuous coupling
- Post-solver becomes more effective

**Trade-off:** 5-10x longer simulation time

### Option 2: Saturation Temperature Enforcement
After each heat addition:
```python
# Get saturation temperature at current pressure
T_sat = get_Tsat(P)

# Enforce physical limits
if T_liquid > T_sat + 2K:
    # Force evaporation
    dm_evap = (T_liquid - T_sat) / h_fg * some_rate

if T_gas < T_sat - 2K:
    # Force condensation
    dm_cond = (T_sat - T_gas) / h_fg * some_rate
```

**Rationale:** Liquid physically cannot exist much above saturation temperature

### Option 3: Interface-Driven Mass Transfer
Model evaporation rate explicitly:
```python
# Mass transfer at interface
A_interface = calculate_interface_area(liquid_level)
h_mass_transfer = 100  # W/m²K (tunable)
dm_evap = h_mass_transfer * A_interface * (T_liquid - T_sat) / h_fg * dt
```

**Rationale:** More physics-based than discrete quality checks

### Option 4: Keep Pre-Solver Only
Revert to original pre-solver approach:
- Best energy conservation (5-13%)
- Most reasonable results (ΔT = 46K)
- Fastest simulation
- Liquid almost fully evaporates (physically reasonable for fire scenario)

## My Recommendation

**Use pre-solver only** with possibly:
- Slightly reduced timestep (0.25s instead of 0.5s) if needed
- Add saturation temperature check as safety limit

The pre-solver approach gave the most physically reasonable results:
- Moderate stratification (46K)
- Most liquid evaporates during fire (expected)
- Good energy conservation
- Simulation completes quickly

The extreme superheat in other approaches (614-662K) suggests they're missing the continuous evaporation that would cool the gas.

## Test Case Context

**Scenario:** Propane vessel in pool fire
- Fire heat flux: 60 kW/m² (API 521)
- PSV opens at 14.3 bar
- Duration: 600 seconds

**Expected behavior:**
- Liquid should evaporate rapidly (fire heat + pressure relief)
- Gas temperature rises but limited by evaporative cooling
- Pressure limited by PSV

**Pre-solver results align best with expectations.**
