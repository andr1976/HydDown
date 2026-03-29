# Non-Equilibrium Model (NEM) Implementation Methodology

## Overview

The Non-Equilibrium Model (NEM), also called Partial Phase Equilibrium (PPE) model, extends the standard equilibrium thermodynamic model to handle two-phase systems where the gas and liquid phases can have different temperatures while maintaining pressure equilibrium.

## Physical Basis

### Key Assumptions

1. **Pressure equilibrium**: Both phases are at the same pressure at all times (P_gas = P_liquid)
2. **Thermal non-equilibrium**: Gas and liquid can have different temperatures (T_gas ≠ T_liquid)
3. **Separate phase tracking**: Each phase has its own mass, temperature, density, and internal energy
4. **No slip velocity**: Phases share the same volume (no relative motion)
5. **Mass conservation**: Total mass = m_gas + m_liquid

### Physical Phenomena Modeled

- **Thermal stratification**: Temperature difference between phases due to limited heat transfer
- **Phase-specific wall heat transfer**: Different heat transfer mechanisms for gas (natural convection) and liquid (nucleate boiling)
- **Gas-liquid interfacial heat transfer**: Parameterized with coefficient h_gl [W/m²K]
- **Phase change**: Evaporation/condensation when phases enter two-phase region

## Implementation Architecture

### State Variables

For each phase, track independent thermodynamic states:

**Gas phase:**
- `m_gas[i]`: Mass [kg]
- `T_gas[i]`: Temperature [K]
- `U_gas[i]`: Specific internal energy [J/kg]
- `rho_gas[i]`: Density [kg/m³]

**Liquid phase:**
- `m_liquid[i]`: Mass [kg]
- `T_liquid[i]`: Temperature [K]
- `U_liquid[i]`: Specific internal energy [J/kg]
- `rho_liquid[i]`: Density [kg/m³]

**Shared variables:**
- `P[i]`: Pressure [Pa] (same for both phases)
- `liquid_level[i]`: Liquid height [m]
- `T_vessel[i]`: Unwetted wall temperature [K]
- `T_vessel_wetted[i]`: Wetted wall temperature [K]

### CoolProp Fluid Objects

Three separate fluid objects are maintained:

1. **`self.fluid`**: Equilibrium mixture state (for compatibility, set to gas state for NEM)
2. **`self.fluid_gas`**: Gas phase properties
3. **`self.fluid_liquid`**: Liquid phase properties

## Time Integration Scheme

### Step 1: Mass Balance (Valve Flow)

Update phase masses based on valve flow:

```
if discharge:
    m_gas[i] = m_gas[i-1] - mdot * dt
    m_liquid[i] = m_liquid[i-1]
else:  # filling
    m_gas[i] = m_gas[i-1] - mdot * dt  (mdot < 0 for filling)
    m_liquid[i] = m_liquid[i-1]
```

### Step 2: Phase Transfer (Evaporation/Condensation)

Check if either phase is in two-phase region and transfer mass:

**Liquid phase check:**
```
fluid_liquid.update(DmassUmass_INPUTS, rho_liquid[i-1], U_liquid[i-1])
Q_liquid = fluid_liquid.Q()

if 0 < Q_liquid < 1:
    # Some liquid should evaporate
    dm_evap = Q_liquid * m_liquid[i-1] * relax_factor
    m_liquid[i] -= dm_evap
    m_gas[i] += dm_evap
    # Transfer energy with latent heat
```

**Gas phase check (if no evaporation):**
```
fluid_gas.update(DmassUmass_INPUTS, rho_gas[i-1], U_gas[i-1])
Q_gas = fluid_gas.Q()

if 0 < Q_gas < 1:
    # Some gas should condense
    dm_cond = (1 - Q_gas) * m_gas[i-1] * relax_factor
    m_gas[i] -= dm_cond
    m_liquid[i] += dm_cond
```

### Step 3: Heat Transfer Calculations

#### 3.1 Gas-Side Wall Heat Transfer

**Heat transfer coefficient:**
```
if non_equilibrium:
    T_for_gas_side = T_gas[i-1]
else:
    T_for_gas_side = T_fluid[i-1]

h_inside = h_inside_calc(L, T_for_gas_side, T_wall_inner, ...)
```

**Heat transfer rate:**
```
A_unwetted = surf_area_inner - wetted_area
Q_inner[i] = A_unwetted * h_inside * (T_wall_inner - T_gas[i-1])
```

#### 3.2 Liquid-Side Wall Heat Transfer (Wetted Region)

**Key principle**: For NEM, nucleate boiling occurs if liquid exists, regardless of quality check.

**Heat transfer coefficient:**
```
if non_equilibrium:
    liquid_exists = m_liquid[i-1] > 1e-6
    if liquid_exists:
        # Update fluid_liquid to saturation for saturated properties
        fluid_liquid.update(PQ_INPUTS, P[i-1], 0.0)

        # Update transport_fluid_wet for transport properties
        T_film_wet = (T_liquid[i-1] + T_wall_wetted[i-1]) / 2
        transport_fluid_wet.update(PT_INPUTS, P[i-1], T_film_wet)

        # Calculate nucleate boiling HTC
        h_inside_wetted = h_inside_wetted_calc(
            L, T_wall_wetted, T_liquid[i-1],
            transport_fluid_wet,  # For k, μ, Cp at film T
            fluid_liquid          # For σ, ρ_sat, h_fg at saturation
        )
    else:
        h_inside_wetted = h_inside  # No liquid, use gas-side
```

**Heat transfer rate:**
```
if non_equilibrium:
    T_fluid_wet = T_liquid[i-1]
else:
    T_fluid_wet = T_fluid[i-1]

Q_inner_wetted[i] = wetted_area * h_inside_wetted * (T_wall_wetted - T_fluid_wet)
```

**Critical implementation detail**: Use `T_liquid` for ΔT calculation, not `T_fluid` (which equals `T_gas` in NEM).

#### 3.3 Gas-Liquid Interfacial Heat Transfer

```
A_interface = inner_vol.A_cross_sectional(liquid_level)
Q_gas_liquid = h_gl * A_interface * (T_gas[i-1] - T_liquid[i-1])
```

where `h_gl` is either specified by user or auto-calculated via
`h_gas_liquid_interface()` natural convection correlation (typical range: 100-2000 W/m²K).
Default if omitted: 50 W/m²K.

### Step 4: Energy Balance

#### Phase Transfer Energy

Energy carried by phase-changing mass uses an average of enthalpy and internal
energy at saturation as a compromise between open-system (enthalpy) and
closed-system (internal energy) formulations:

```
# Evaporation: energy carried by mass leaving liquid as vapor
e_vap_sat = (h_vap_sat + u_vap_sat) / 2
E_evap = dm_evap * e_vap_sat

# Condensation: energy carried by mass leaving gas as liquid
e_liq_sat = (h_liq_sat + u_liq_sat) / 2
E_cond = dm_cond * e_liq_sat
```

#### 4.1 Gas Phase Energy Balance

```
U_gas_start = U_gas[i-1] * m_gas[i-1]

# Enthalpy flow (discharge/filling through gas phase)
if discharge:
    h_out = h_gas_discharge  # Vapor enthalpy at current state
    H_flow = -mdot[i-1] * h_out * dt
else:
    h_in = h_reservoir
    H_flow = -mdot[i-1] * h_in * dt  (mdot < 0)

# Heat transfers
Q_wall_to_gas = Q_inner[i] * dt
Q_interface = Q_gas_liquid * dt  (positive = gas to liquid)
E_phase_change = E_evap - E_cond  (energy from evaporation/condensation)

U_gas_end = U_gas_start + H_flow + Q_wall_to_gas - Q_interface + E_phase_change
U_gas_specific = U_gas_end / m_gas[i]
```

#### 4.2 Liquid Phase Energy Balance

```
U_liquid_start = U_liquid[i-1] * m_liquid[i-1]

# Heat transfers
Q_wall_to_liquid = Q_inner_wetted[i] * dt
Q_interface = Q_gas_liquid * dt  (positive = liquid receives heat from gas)
E_phase_change = -E_evap + E_cond  (opposite sign from gas)

U_liquid_end = U_liquid_start + Q_wall_to_liquid + Q_interface + E_phase_change
U_liquid_specific = U_liquid_end / m_liquid[i]
```

### Step 5: Thermodynamic State Update (Pressure Solver)

**Problem**: Given specific internal energies (U_gas_specific, U_liquid_specific) and masses, find the common pressure P where both phases fit in the vessel volume.

**Constraint**: Volume conservation
```
V_total = m_gas / rho_gas + m_liquid / rho_liquid
```

**Iterative solver** (Brent's method):
```
def volume_residual(P):
    # Update gas phase at (P, U_gas_specific)
    fluid_gas.update(PUmass_INPUTS, P, U_gas_specific)
    rho_gas = fluid_gas.rhomass()
    T_gas[i] = fluid_gas.T()

    # Update liquid phase at (P, U_liquid_specific)
    fluid_liquid.update(PUmass_INPUTS, P, U_liquid_specific)
    rho_liquid = fluid_liquid.rhomass()
    T_liquid[i] = fluid_liquid.T()

    # Calculate volume residual
    V_calculated = m_gas[i] / rho_gas + m_liquid[i] / rho_liquid
    return V_calculated - V_total

# Solve for P
P[i] = brentq(volume_residual, P_min, P_max)
```

**Pressure bounds**: Based on previous timestep pressure
```
P_min = max(P[i-1] * 0.5, 1e5)   # At least 1 bar
P_max = P[i-1] * 2.0
```

**Fallback solver chain**: brentq → ridder → bisect (with expanded bounds 0.2x to 5x)

### Step 6: Geometric Calculations

All geometric calculations use the `fluids.TANK` object (`self.inner_vol`).

**Liquid level** (from liquid density and mass):
```
V_liquid = m_liquid[i] / rho_liquid[i]
liquid_level[i] = inner_vol.h_from_V(V_liquid)
```

**Wetted area** (liquid contact with wall):
```
wetted_area = inner_vol.SA_from_h(liquid_level[i])
```

**Gas-liquid interface area** (cross-sectional area at liquid surface):
```
A_interface = inner_vol.A_cross_sectional(liquid_level[i])
```

These methods handle all supported head types (hemispherical, ASME F&D, DIN,
semi-elliptical, flat-end) and both horizontal/vertical orientations.

### Step 7: Wall Temperature Update

**Unwetted wall** (exposed to gas):
```
T_vessel[i] = T_vessel[i-1] + (Q_outer - Q_inner) * dt / (m_wall * cp_wall * fraction_unwetted)
```

**Wetted wall** (exposed to liquid):
```
if wetted_area > 1e-10:
    T_vessel_wetted[i] = T_vessel_wetted[i-1] + (Q_outer_wetted - Q_inner_wetted) * dt / (m_wall * cp_wall * fraction_wetted)
else:
    T_vessel_wetted[i] = T_vessel_wetted[i-1]  # Safety check
```

## Key Implementation Details

### 1. Boiling Heat Transfer for NEM

**Old (incorrect) approach:**
```
if fluid.Q() >= 0 and fluid.Q() <= 1:  # Check equilibrium quality
    hiw = h_inside_wetted(...)
```

**New (correct) approach:**
```
if non_equilibrium:
    if m_liquid > 1e-6:  # Liquid exists
        hiw = h_inside_wetted(...)  # Always use boiling correlation
```

**Rationale**: In NEM, phases are tracked separately. If liquid exists, it IS liquid and can boil when wall is hot, regardless of whether equilibrium mixture would be two-phase.

### 2. Temperature Selection for Heat Transfer

| Surface | Equilibrium Mode | NEM Mode |
|---------|-----------------|----------|
| **Unwetted wall → gas** | T_fluid | T_gas (explicit) |
| **Wetted wall → liquid** | T_fluid | T_liquid (explicit) |
| **Gas ↔ liquid interface** | N/A | T_gas - T_liquid |

### 3. Fluid Object State Management

**For boiling heat transfer:**
```
# Step A: Update to saturation for saturated properties
fluid_liquid.update(PQ_INPUTS, P, 0.0)

# Step B: Use in h_inside_wetted (needs σ, ρ_sat, h_fg)
hiw = h_inside_wetted(..., fluid_liquid)

# Step C: Later restored to actual state for phase transfer check
fluid_liquid.update(DmassUmass_INPUTS, rho_liquid, U_liquid)
```

**Safety**: Between steps A and C, `fluid_liquid` is only passed to `h_inside_wetted()`, which requires saturation state. No property queries (rho, U, T, h) are made.

### 4. Backward Compatibility

All NEM-specific logic is wrapped in conditionals:
```
if self.non_equilibrium:
    # NEM logic
else:
    # Equilibrium logic (unchanged)
```

This ensures equilibrium calculations remain identical to pre-NEM implementation.

## Input File Configuration

### Required Parameters

```yaml
calculation:
  type: "energybalance"
  non_equilibrium: true     # Enable NEM
  h_gas_liquid: 500         # Gas-liquid HTC [W/m²K] - fixed value
  # or
  h_gas_liquid: "calc"      # Auto-calculate using natural convection correlation

vessel:
  liquid_level: 0.4668      # Initial liquid level [m] (required for two-phase)
  # ... other vessel geometry parameters ...

initial:
  temperature: 278          # Bulk/saturation temperature [K]
  pressure: 550000          # Initial pressure [Pa]
  fluid: "propane"          # Single component only
```

**Note on initial temperatures:** Gas is automatically initialized at `T_sat + 5 K`
(slight superheat) and liquid at `T_sat` at the given pressure.
There are no separate `temperature_gas`/`temperature_liquid` input parameters.

### Optional Parameters

```yaml
calculation:
  h_gas_liquid: "calc"      # Auto-calculate h_gl (default if omitted: 50 W/m²K)
```

**Note:** The phase transfer relaxation factor is hardcoded at 0.8 and is not
configurable via the input file. The wetted area and interface area are computed
automatically from the vessel geometry using `fluids.TANK`.

## Validation and Diagnostics

### Key Checks During Simulation

1. **Mass conservation**: `m_total = m_gas + m_liquid` (constant if no valve flow)
2. **Pressure equilibrium**: `P_gas == P_liquid` (enforced by solver)
3. **Volume constraint**: `V = m_gas/rho_gas + m_liquid/rho_liquid` (enforced by solver)
4. **Energy balance closure**: Track energy in/out to verify conservation

### Common Physical Behaviors

- **Early heating phase**: T_gas rises faster than T_liquid (stratification develops)
- **PSV cycling**: Depressurization reduces T_gas, narrows gas-liquid temperature gap
- **Fire scenario**: Unwetted wall reaches high temperature, wetted wall stays cooler
- **Liquid boiling**: Wetted wall temperature intermediate between T_gas and T_liquid

### Diagnostic Outputs

Key arrays for analysis:
- `hd.T_gas`, `hd.T_liquid`: Phase temperatures [K]
- `hd.m_gas`, `hd.m_liquid`: Phase masses [kg]
- `hd.U_gas`, `hd.U_liquid`: Phase specific internal energies [J/kg]
- `hd.rho_gas`, `hd.rho_liquid`: Phase densities [kg/m³]
- `hd.Q_inner`, `hd.Q_inner_wetted`: Wall heat transfer rates [W]
- `hd.h_gas_liquid`: Gas-liquid interfacial HTC [W/m²K]
- `hd.h_inside`, `hd.h_inside_wetted`: Wall heat transfer coefficients [W/m²K]
- `hd.mdot_phase_transfer`: Phase transfer rate [kg/s] (positive = condensation)
- `hd.liquid_level`: Liquid level height [m]

## Limitations and Future Work

### Current Limitations

1. **Single component only**: Multi-component mixtures not supported with NEM
2. **No slip assumption**: Phases share same volume, no sloshing/entrainment
3. **Simplified phase transfer**: Evaporation/condensation uses relaxation factor (0.8), not rate-limited
4. **Natural convection only for h_gl**: Auto-calculated h_gl uses natural convection correlation; forced convection effects (jetting, sloshing) are not included

### Potential Enhancements

1. **Multi-component NEM**: Extend to mixtures with different phase compositions
2. **Superheating/subcooling limits**: Physical bounds on temperature excursions
3. **Enhanced phase transfer**: Rate-limited models based on interfacial area and mass transfer coefficient
4. **Forced/mixed convection for h_gl**: Account for PSV discharge jetting and sloshing effects

## Summary

The NEM implementation provides a rigorous framework for modeling thermal non-equilibrium in two-phase systems by:

1. **Tracking separate phase states** (T, m, U, ρ for each phase)
2. **Using phase-specific properties** for all heat transfer calculations
3. **Maintaining pressure equilibrium** through iterative volume constraint solver
4. **Explicitly modeling all energy pathways**: wall-gas, wall-liquid, gas-liquid, valve flow, phase change
5. **Ensuring backward compatibility** with equilibrium model for validation

The methodology is physically consistent, numerically stable, and enables accurate simulation of fire scenarios, PSV cycling, and other transients where thermal stratification significantly affects vessel response.
