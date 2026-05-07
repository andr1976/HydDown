# Phase Transfer Analysis: Pre-Solver vs Post-Solver

## Implementations Tested

### 1. Pre-Solver Phase Transfer (Original)
**Method:**
- After heat transfer, check if tentative U puts phase in two-phase region
- Transfer mass/energy BEFORE solving for pressure
- Solve P with modified masses and energies

**Results:**
- Final T_gas: 385 K
- Final T_liquid: 346 K
- ΔT: 46 K
- Phase transfer: ~0.02-0.03 kg/timestep (moderate)
- Energy residuals: 5-13%

**Issue:**
- Checking phase change at old density with new U is artificial
- May over-predict phase transfer

### 2. Post-Solver Phase Transfer Only (Current)
**Method:**
- Solve for P with heat-only energies (no phase transfer)
- AFTER P-solver, check if final (P,U) state is in two-phase
- Transfer mass/energy and re-solve if needed

**Results:**
- Final T_gas: 662 K (!)
- Final T_liquid: 306 K
- ΔT: 356 K (extreme superheat!)
- Phase transfer: ~0.001-0.002 kg/timestep (very low)
- Energy residuals: 8-10%

**Issue:**
- Gas becomes extremely superheated because liquid barely evaporates
- Physically unrealistic stratification
- Once at equilibrium (P,U), phases are stable and don't transfer

## Physical Interpretation

### Why Pre-Solver Gives More Transfer

When we check `update(DmassUmass, ρ_old, U_new)`:
- If U increased significantly, this may show quality > 0
- Triggers evaporation
- But this is checking at OLD density, which is inconsistent

### Why Post-Solver Gives Less Transfer

When we check `update(PUmass, P_solved, U)`:
- This is the thermodynamically consistent state
- If both (P,U) are from energy balance, phases are in equilibrium
- Very little phase change detected

The issue: **Without mass transfer during heating, gas and liquid evolve independently and get farther apart**

## The Real Physics

In a real vessel with two phases:
1. Liquid receives heat → temperature rises
2. As liquid heats, vapor pressure increases
3. Some liquid evaporates to maintain saturation
4. Evaporation cools the liquid (latent heat)
5. This happens **continuously**, not in discrete jumps

Our explicit Euler method with discrete timesteps can't capture this continuous coupling.

## Possible Solutions

### Option A: Hybrid Approach
- Keep pre-solver phase transfer (with relaxation)
- ALSO check post-solver for consistency
- If post-solver detects opposite phase change, reduce pre-solver transfer

### Option B: Saturation Temperature Coupling
- Calculateexpected T_sat from current P
- If liquid T > T_sat: force evaporation
- If gas T < T_sat: force condensation
- Use temperature deviation as driver, not quality

### Option C: Interface Mass Transfer Rate
- Model evaporation rate explicitly: dm/dt = h_gl * A * (T_liq - T_sat) / h_fg
- Independent of thermodynamic state checking
- More physics-based but requires tuning h_gl

### Option D: Smaller Timesteps with Post-Solver Only
- Current timestep: 0.5s might be too large
- Reduce to 0.1s or 0.05s
- Post-solver transfer becomes more effective
- But 10x slower simulation

## Recommendation

I suggest **Option B** - saturation temperature coupling:

```python
# After heat transfer, before P-solver
P_guess = self.P[i-1]

# Get saturation temperature at current pressure
fluid_sat = CP.AbstractState("HEOS", species)
fluid_sat.update(CP.PQ_INPUTS, P_guess, 0.0)
T_sat = fluid_sat.T()

# Check liquid superheat
if T_liquid_tentative > T_sat + tolerance:
    # Liquid is superheated → must evaporate
    dm_evap = calculate_evaporation_rate(...)

# Check gas subcooling
if T_gas_tentative < T_sat - tolerance:
    # Gas is subcooled → must condense
    dm_cond = calculate_condensation_rate(...)
```

This enforces physical constraint: **liquid cannot exist above saturation temperature at given pressure**.

## Energy Conservation Note

Both approaches show 5-13% energy residuals, which suggests the issue is not phase transfer method but:
1. Explicit Euler truncation error
2. Discharge enthalpy estimation
3. Gas-liquid heat transfer approximation
4. Relaxation factor damping

These are acceptable for engineering calculations.

## Current Status

The post-solver-only approach is implemented but gives unrealistic superheat. We need to add saturation-based phase transfer or revert to hybrid approach.
