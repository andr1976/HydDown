# NEM Solver Stability Analysis

## Current Implementation (hdclass.py:2647-2651)

### Solver Setup
```python
# Bounds for pressure
P_min = max(self.P[i-1] * 0.5, 1e5)  # At least 1 bar
P_max = self.P[i-1] * 2.0

# Solve for pressure that gives correct volume
P_solution = brentq(volume_residual, P_min, P_max, xtol=1e-6)
```

### Key Observations

#### 1. **Starting Guess Strategy**
- **Current**: No explicit initial guess provided to `brentq`
- **Behavior**: `brentq` evaluates function at both bounds `[P_min, P_max]`
- **Issue**: If previous pressure is a good guess, it's not being utilized
- **Better approach**: Could use `P[i-1]` as initial guess with hybrid solver

#### 2. **Bounds Strategy**
- **Current**: Symmetric bounds around previous pressure (0.5x to 2x)
- **Lower bound**: `max(0.5 * P[i-1], 1 bar)` - minimum of 1 bar
- **Upper bound**: `2.0 * P[i-1]` - double previous pressure
- **Issues**:
  - Fixed multipliers (0.5, 2.0) may be too restrictive
  - When PSV is opening/closing, pressure can change rapidly
  - During phase transitions, pressure can jump discontinuously
  - No adaptive adjustment based on residual signs

#### 3. **Fallback Solver**
- **Current**: NONE - fails hard if `brentq` doesn't converge
- **Exception**: Returns `1.0` if CoolProp update fails (line 2628)
- **Issue**: No attempt to try alternative solvers or bound adjustments
- **Better approach**: Could try:
  - `scipy.optimize.ridder` - similar to brentq but more robust
  - `scipy.optimize.bisect` - slower but guaranteed convergence if bracketed
  - Adaptive bound expansion if root not bracketed

#### 4. **Root Bracketing**
- **Error message**: `"f(a) and f(b) must have different signs"`
- **Meaning**: The root is NOT in the interval `[P_min, P_max]`
- **Causes**:
  - Actual pressure outside [0.5*P_prev, 2*P_prev] range
  - Rapid pressure changes during PSV cycling
  - Large timestep causing too much change per step
  - CoolProp thermodynamic limits being exceeded

#### 5. **Failure Scenario Analysis**

From error at t=348.9s with PSV diameter 0.015m:
```
P_guess range: [1058804.51, 4235218.05] Pa  (10.6 to 42.4 bar)
m_gas: 89.359 kg, m_liquid: 824.454 kg
```

- Previous pressure ~21 bar (midpoint of guess range)
- PSV set pressure: 15.6 bar
- Likely scenario: PSV just opened, causing rapid depressurization
- Actual pressure may have dropped below 10.6 bar (P_min)
- Or thermodynamic state became inconsistent with energy balance

## Potential Improvements (NOT YET IMPLEMENTED)

### Option A: Use Previous Pressure as Initial Guess
```python
from scipy.optimize import brentq, ridder

try:
    # First try brentq with current bounds
    P_solution = brentq(volume_residual, P_min, P_max, xtol=1e-6)
except ValueError:
    # If not bracketed, try to expand bounds adaptively
    # Check residual at P_prev
    res_prev = volume_residual(self.P[i-1])

    if abs(res_prev) < 0.01:  # Close to solution
        P_solution = self.P[i-1]
    else:
        # Expand bounds based on residual sign
        if res_prev > 0:
            P_max = min(P_max * 2.0, 100e5)  # Need higher P
        else:
            P_min = max(P_min * 0.5, 1e5)  # Need lower P

        P_solution = brentq(volume_residual, P_min, P_max, xtol=1e-6)
```

### Option B: Adaptive Bounds with Wider Initial Range
```python
# More generous bounds during PSV operation
if abs(self.mass_rate[i-1]) > 1e-6:  # PSV is active
    P_min = max(self.P[i-1] * 0.3, 1e5)
    P_max = self.P[i-1] * 3.0
else:
    P_min = max(self.P[i-1] * 0.5, 1e5)
    P_max = self.P[i-1] * 2.0
```

### Option C: Fallback Solver Chain
```python
from scipy.optimize import brentq, ridder, bisect

try:
    # Try brentq first (fastest)
    P_solution = brentq(volume_residual, P_min, P_max, xtol=1e-6, maxiter=100)
except ValueError:
    try:
        # Try ridder (more robust)
        P_solution = ridder(volume_residual, P_min, P_max, xtol=1e-6, maxiter=100)
    except ValueError:
        # Check if we can bracket the root by expanding bounds
        P_min_expanded = max(P_min * 0.1, 1e5)
        P_max_expanded = min(P_max * 10.0, 100e5)

        try:
            P_solution = bisect(volume_residual, P_min_expanded, P_max_expanded, xtol=1e-6)
        except ValueError:
            # Last resort: use previous pressure if residual is small
            res_prev = volume_residual(self.P[i-1])
            if abs(res_prev) < 0.05:
                P_solution = self.P[i-1]
            else:
                raise  # Give up
```

### Option D: Reduce Timestep Automatically
```python
# Detect rapid changes and reduce timestep
if abs(self.P[i-1] - self.P[i-2]) / self.P[i-2] > 0.1:  # >10% change
    print(f"Warning: Large pressure change detected, consider reducing timestep")
    # Could automatically subdivide this timestep
```

## Comparison of Scipy Root-Finding Methods

| Method | Speed | Robustness | Notes |
|--------|-------|------------|-------|
| `brentq` | Fast | Good if bracketed | Current implementation |
| `ridder` | Medium | Better than brentq | Good for noisy functions |
| `bisect` | Slow | Very robust | Guaranteed if bracketed |
| `newton` | Very fast | Poor | Needs derivative, can diverge |
| `fsolve` | Fast | Medium | Can handle unbrackt but less reliable |

## Recommendations (Priority Order)

1. **Add diagnostic output** to understand failure modes better
2. **Implement adaptive bounds** that expand when PSV is active
3. **Add fallback solver** (ridder → bisect chain)
4. **Use previous pressure** to check if near solution before full solve
5. **Warn user** if timestep appears too large based on pressure changes
6. **Consider** using `scipy.optimize.root_scalar` with hybrid method

## Current Workarounds (User-Side)

1. **Reduce timestep**: `time_step: 0.05` or smaller
2. **Larger PSV**: Avoids rapid pressure oscillations
3. **Lower fire heat flux**: Slower pressure changes
4. **Monitor** pressure rate of change - if >20% per step, reduce dt

## References

- scipy.optimize.brentq: https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.brentq.html
- scipy.optimize.ridder: https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.ridder.html
- Root-finding comparison: https://docs.scipy.org/doc/scipy/reference/optimize.html#root-finding
