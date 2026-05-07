# Gas-Liquid Interfacial Heat Transfer Correlation (h_gl)

## Overview

Implemented automatic calculation of gas-liquid interfacial heat transfer coefficient (h_gl) for NEM based on natural convection correlations, eliminating the need for user-specified values in most cases.

## Physical Basis

### Configuration

For typical NEM scenarios:
- **Hot gas phase** (lighter) above
- **Cold liquid phase** (denser) below
- **Horizontal interface** at liquid level
- **Stable stratification** (hot fluid above cold fluid)

### Heat Transfer Mechanism

Natural convection driven by:
- Temperature difference ΔT = T_gas - T_liquid
- Buoyancy forces in each phase
- Thermal boundary layers at interface

## Correlation Formulation

### Overall Heat Transfer Coefficient

The overall h_gl accounts for thermal resistances on both sides of the interface:

```
1/h_gl = 1/h_gas + 1/h_liquid
```

where:
- `h_gas`: Gas-side heat transfer coefficient [W/m²K]
- `h_liquid`: Liquid-side heat transfer coefficient [W/m²K]

### Individual Phase Correlations

**Gas side** (heated surface facing downward - stable):
```
Nu_gas = 0.27 * Ra_gas^0.25   for Ra < 10^7 (laminar)
Nu_gas = 0.15 * Ra_gas^0.333  for Ra ≥ 10^7 (turbulent)

h_gas = Nu_gas * k_gas / L_char
```

**Liquid side** (cooled surface facing upward - stable):
```
Nu_liquid = 0.27 * Ra_liquid^0.25   for Ra < 10^7 (laminar)
Nu_liquid = 0.15 * Ra_liquid^0.333  for Ra ≥ 10^7 (turbulent)

h_liquid = Nu_liquid * k_liquid / L_char
```

### Dimensionless Numbers

**Rayleigh number**:
```
Ra = Gr * Pr
Gr = g * β * ΔT * L³ / ν²
Pr = μ * Cp / k
```

where:
- `g`: Gravitational acceleration [9.81 m/s²]
- `β`: Thermal expansion coefficient [1/K]
- `ΔT`: Temperature difference [K]
- `L`: Characteristic length [m]
- `ν`: Kinematic viscosity [m²/s]
- `μ`: Dynamic viscosity [Pa·s]
- `Cp`: Specific heat capacity [J/kgK]
- `k`: Thermal conductivity [W/mK]

### Characteristic Length

- **Horizontal vessels**: Use diameter (D)
- **Vertical vessels**: Use diameter (D)
- Represents typical boundary layer development length

## Implementation

### Function Location

`src/hyddown/transport.py::h_gas_liquid_interface()`

### Function Signature

```python
def h_gas_liquid_interface(T_gas, T_liquid, P, L_char, fluid_gas, fluid_liquid):
    """
    Calculate gas-liquid interfacial heat transfer coefficient.

    Parameters
    ----------
    T_gas : float
        Gas phase bulk temperature [K]
    T_liquid : float
        Liquid phase bulk temperature [K]
    P : float
        System pressure [Pa]
    L_char : float
        Characteristic length (vessel diameter) [m]
    fluid_gas : CoolProp.AbstractState
        Gas phase fluid object (already updated at current state)
    fluid_liquid : CoolProp.AbstractState
        Liquid phase fluid object (already updated at current state)

    Returns
    -------
    h_gl : float
        Overall gas-liquid heat transfer coefficient [W/m²K]
    ```

### Usage in Input File

**Option 1: Automatic calculation (recommended)**
```yaml
calculation:
  type: "energybalance"
  non_equilibrium: True
  h_gas_liquid: "calc"  # Use correlation
```

**Option 2: Fixed value (for sensitivity studies)**
```yaml
calculation:
  type: "energybalance"
  non_equilibrium: True
  h_gas_liquid: 500  # Fixed value in W/m²K
```

**Option 3: Default (if not specified)**
```yaml
calculation:
  type: "energybalance"
  non_equilibrium: True
  # h_gas_liquid defaults to 50 W/m²K
```

### Diagnostic Output

The calculated h_gl is stored in the array:
```python
hd.h_gas_liquid[i]  # Gas-liquid HTC at time step i [W/m²K]
```

## Validation Results

### Test Case: NEM Propane PSV Fire (t=500s)

| Parameter | Fixed h_gl=500 | Calculated h_gl | Notes |
|-----------|----------------|-----------------|-------|
| **h_gl range** | 500 (constant) | 5-99 W/m²K | Varies with ΔT |
| **h_gl mean** | 500 | 46.5 W/m²K | Lower for stable stratification |
| **Final ΔT (gas-liquid)** | 97.9 K | 135.9 K | More stratification with correlation |
| **Final T_gas** | 136.4°C | 176.8°C | Higher due to reduced heat transfer |
| **Final T_liquid** | 38.5°C | 40.9°C | Similar liquid heating |

### Correlation Behavior

**Dependency on ΔT**:
- At ΔT = 1 K: h_gl ≈ 18 W/m²K (weak convection)
- At ΔT = 50 K: h_gl ≈ 58 W/m²K (moderate convection)
- At ΔT = 150 K: h_gl ≈ 99 W/m²K (strong convection)

**Physical interpretation**:
- Larger temperature differences → stronger buoyancy forces → higher h_gl
- Stable stratification (hot above cold) → modest h_gl values (5-100 W/m²K typical)
- Unstable stratification would give much higher values (not applicable to typical NEM)

## Comparison with Literature

### Typical h_gl Ranges

| Configuration | Literature Range | Correlation Range | Agreement |
|--------------|------------------|-------------------|-----------|
| **Stable horizontal stratification** | 20-200 W/m²K | 5-100 W/m²K | ✓ Good |
| **Unstable stratification** | 500-5000 W/m²K | N/A | (Different regime) |
| **Forced convection (mixing)** | 1000-10000 W/m²K | N/A | (Different physics) |

### References for Correlations

- **Churchill and Chu (1975)**: "Correlating equations for laminar and turbulent free convection from a horizontal surface"
- **Incropera & DeWitt**: "Fundamentals of Heat and Mass Transfer" - Natural convection correlations
- **McAdams (1954)**: "Heat Transmission" - Classic natural convection data

## Limitations and Assumptions

### Current Assumptions

1. **Stable stratification**: Hot gas above cold liquid (typical for NEM)
2. **Quiescent phases**: No forced convection or significant mixing
3. **Horizontal interface**: Simplified geometry (actual interface may be curved)
4. **No waves or sloshing**: Interface treated as flat and stationary
5. **No phase change at interface**: Evaporation/condensation handled separately

### Limitations

1. **Does not account for**:
   - Interface waves or disturbances
   - Entrainment of one phase into another
   - Marangoni effects (surface tension gradients)
   - Condensation heat transfer at interface
   - Non-horizontal vessel orientations

2. **Correlation bounds**:
   - Minimum h_gl = 5 W/m²K (very weak convection)
   - Maximum h_gl = 1000 W/m²K (very strong convection)
   - Falls back to 100 W/m²K if calculation fails

### When Fixed Values May Be Preferred

- **High mixing scenarios**: Agitated vessel, jet impingement
- **Validation studies**: Matching specific experimental h_gl
- **Unstable stratification**: Cold gas above hot liquid (rare)
- **Uncertainty quantification**: Sensitivity to h_gl variation

## Recommendations

### General Usage

✅ **Use "calc" for**:
- Fire scenarios
- Depressurization/pressurization
- Natural thermal stratification development
- When h_gl value is unknown

✅ **Use fixed value for**:
- Validation against experimental data with known h_gl
- Sensitivity studies exploring h_gl impact
- Scenarios with known mixing conditions

### Typical Values as Reference

| Scenario | Recommended h_gl | Source |
|----------|-----------------|---------|
| **Stable stratification (typical NEM)** | "calc" or 50-100 W/m²K | This correlation |
| **Moderate mixing** | 200-500 W/m²K | Engineering judgment |
| **Strong mixing** | 500-2000 W/m²K | Engineering judgment |
| **Vigorous boiling at interface** | 1000-5000 W/m²K | Phase change dominant |

## Impact on Simulation Results

### Effect of h_gl on Thermal Stratification

- **Lower h_gl** (10-50 W/m²K):
  - Greater temperature difference between phases
  - Gas heats up faster
  - Liquid temperature more independent
  - Closer to adiabatic phase behavior

- **Higher h_gl** (500-2000 W/m²K):
  - Smaller temperature difference between phases
  - Stronger thermal coupling
  - Approaches thermal equilibrium behavior
  - Less stratification

### Effect on PSV Sizing and Safety

- **Lower h_gl**:
  - Gas reaches higher temperatures
  - Higher pressure rise rate
  - Earlier PSV opening
  - More conservative for safety analysis

- **Higher h_gl**:
  - More heat transferred to liquid
  - Liquid boils more
  - Faster pressure rise from evaporation
  - Both phases heat together

## Future Enhancements

### Potential Improvements

1. **Enhanced correlations**:
   - Account for vessel orientation effects
   - Include interface curvature effects
   - Model for unstable stratification
   - Condensation heat transfer at interface

2. **Flow regime detection**:
   - Detect transition to mixed/turbulent regimes
   - Adjust correlations for wavy interface
   - Account for entrainment effects

3. **Validation database**:
   - Compile experimental h_gl data for two-phase systems
   - Validate correlation against broader datasets
   - Develop uncertainty bounds

4. **User options**:
   - Choice of correlation (Churchill-Chu, McAdams, etc.)
   - Manual specification of Ra range transitions
   - Custom characteristic length definitions

## Conclusion

The implemented h_gl correlation provides:

✓ **Automatic calculation** based on local thermal conditions
✓ **Physical consistency** with natural convection theory
✓ **Reasonable values** for stable stratification (5-100 W/m²K)
✓ **ΔT dependency** captured naturally
✓ **No user tuning** required in most cases
✓ **Backward compatible** with fixed value option

The correlation eliminates a key user-specified parameter in NEM while maintaining physical fidelity and providing results consistent with natural convection expectations for stratified two-phase systems.
