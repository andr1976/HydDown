# NEM Solver Comparison: brentq vs Nelder-Mead

## Approaches Tested

### 1. brentq-based solver (FAST)
**Method:**
- Outer loop: `brentq` solves for pressure P
- Objective: volume_residual(P) = (V_gas + V_liquid - V_vessel) / V_vessel
- Inner loops: For given P, use `brentq` to find ρ_gas and ρ_liquid from (U,P) constraints
- Each density iteration: solve ρ such that update(DmassUmass, ρ, U) gives correct P

**Performance:**
- Simulation time: ~4 seconds (1200 timesteps)
- Energy residuals: 5-13%
- Stable and reliable

**Code structure:**
```python
def volume_residual(P):
    # Find ρ_gas: brentq(rho_from_UP_gas, bounds)
    # Find ρ_liquid: brentq(rho_from_UP_liquid, bounds)
    # Return (V_total - V_vessel) / V_vessel

P_solution = brentq(volume_residual, P_min, P_max)
```

### 2. Nelder-Mead solver (SLOW)
**Method:**
- Single optimization: `minimize(..., method='Nelder-Mead')`
- Objective: volume_objective(P) = ((V_total - V_vessel) / V_vessel)²
- Inner loops: For given P, use `fsolve` to find ρ_gas and ρ_liquid from (U,P) constraints
- Each density iteration: solve ρ such that update(DmassUmass, ρ, U) gives correct P

**Performance:**
- Simulation time: ~180 seconds (1200 timesteps) - **40x slower!**
- Energy residuals: 5-13% (similar to brentq)
- Warnings about slow convergence in fsolve

**Code structure:**
```python
def volume_objective(P_guess):
    # Find ρ_gas: fsolve(rho_residual_gas, init)
    # Find ρ_liquid: fsolve(rho_residual_liquid, init)
    # Return ((V_total - V_vessel) / V_vessel)²

result = minimize(volume_objective, x0=[P_guess], method='Nelder-Mead', bounds=...)
P_solution = result.x[0]
```

## Why is Nelder-Mead slower?

1. **Triple nesting:** Nelder-Mead calls objective function many times, each call solves two fsolve problems
2. **No gradient info:** Simplex methods like Nelder-Mead explore the function landscape without derivative information
3. **Many evaluations:** Nelder-Mead typically requires 100-200 function evaluations per optimization
4. **fsolve overhead:** Using fsolve instead of brentq adds Newton iteration overhead

## Recommendation

**Use the brentq-based solver** because:
- 40x faster (4s vs 180s)
- Same energy conservation accuracy (5-13%)
- More robust (bracketing method guarantees convergence if solution exists)
- Cleaner separation: outer brentq for P, inner brentq for densities

The user's suggestion to try Nelder-Mead was implemented and tested, but the performance penalty is too high for practical use.

## Results Summary

Both solvers produce similar physics:
- Final pressure: ~20.7 bar
- Final gas temperature: ~387 K
- Final liquid temperature: ~341 K
- Temperature difference: ~46 K
- Energy residuals: 5-13%

The brentq approach is clearly superior for this application.
