# Phase Transfer Energy Analysis

## Current Definitions (using internal energy u)

### Line 2265 - Evaporation Energy
```python
E_evap = dm_evap_mass * (u_vap_sat - u_liq_sat)
```

### Line 2304 - Condensation Energy
```python
E_cond = dm_cond_mass * (u_liq_sat - u_vap_sat)
```

## How They're Used

### Gas Phase Energy Balance (lines 2365, 2378)
```python
U_gas_new = U_gas_old + ... + E_evap - E_cond
```

### Liquid Phase Energy Balance (line 2386)
```python
U_liquid_new = U_liquid_old + ... - E_evap + E_cond
```

## Sign Analysis

### Evaporation (dm > 0, liquid → gas)

**Energy value:**
- E_evap = dm × (u_vap - u_liq)
- Since u_vap > u_liq: E_evap > 0 ✓ (positive)

**Applied to phases:**
- Gas: +E_evap → Gas gains energy ✓
- Liquid: -E_evap → Liquid loses energy ✓

**Physical interpretation:**
- Evaporating mass leaves liquid carrying u_liq
- Evaporating mass enters gas at u_vap
- Net: Gas gains dm×u_vap, Liquid loses dm×u_liq
- But we're tracking the DIFFERENCE: dm×(u_vap - u_liq)

**Is this correct?** Let's think carefully:
- Liquid before: U_liq_old, m_liq_old
- Liquid after: U_liq_new, m_liq_new = m_liq_old - dm

If dm mass leaves at specific internal energy u_liq:
```
U_liq_new = U_liq_old - dm * u_liq
```

But we're doing:
```
U_liq_new = U_liq_old - E_evap = U_liq_old - dm*(u_vap - u_liq)
```

This is WRONG! We should subtract dm×u_liq, not dm×(u_vap - u_liq).

### Condensation (dm > 0, gas → liquid)

**Energy value:**
- E_cond = dm × (u_liq - u_vap)
- Since u_liq < u_vap: E_cond < 0 ✗ (negative!)

**Applied to phases:**
- Gas: -E_cond = -(negative) = **positive** (Gas GAINS energy) ✗ WRONG
- Liquid: +E_cond = +(negative) = **negative** (Liquid LOSES energy) ✗ WRONG

**Should be:**
- Gas should LOSE energy (mass leaving + releases latent heat)
- Liquid should GAIN energy (mass entering + absorbs latent heat)

## The Fundamental Issue

The current approach treats E_evap and E_cond as the "net energy difference" between phases. But this is not how mass transfer works!

**What should happen:**

### For Evaporation (liquid → gas, dm > 0):

Liquid loses:
```
ΔU_liquid = -dm * u_liq  (mass leaves carrying u_liq)
```

Gas gains:
```
ΔU_gas = +dm * u_vap  (mass enters at u_vap)
```

Total system:
```
ΔU_total = -dm*u_liq + dm*u_vap = dm*(u_vap - u_liq) = dm*u_fg
```

But this is NOT zero! Energy is conserved only if we account for the heat required for phase change:
```
Q_phase_change = dm * u_fg  (must be supplied externally)
```

### For Condensation (gas → liquid, dm > 0):

Gas loses:
```
ΔU_gas = -dm * u_vap  (mass leaves carrying u_vap)
```

Liquid gains:
```
ΔU_liquid = +dm * u_liq  (mass enters at u_liq)
```

Total system:
```
ΔU_total = -dm*u_vap + dm*u_liq = -dm*(u_vap - u_liq) = -dm*u_fg
```

This releases energy (negative ΔU means energy released as heat).

## Correct Formulation

**Define phase transfer energy as the energy carried by the transferring mass:**

### For Evaporation:
```python
# Liquid loses energy equal to what the departing mass carried
E_liquid_loss = dm_evap * u_liq_sat

# Gas gains energy equal to what the arriving mass brings
E_gas_gain = dm_evap * u_vap_sat

# In energy balance:
U_liquid_new = U_liquid_old - E_liquid_loss = U_liquid_old - dm_evap * u_liq_sat
U_gas_new = U_gas_old + E_gas_gain = U_gas_old + dm_evap * u_vap_sat
```

### For Condensation:
```python
# Gas loses energy equal to what the departing mass carried
E_gas_loss = dm_cond * u_vap_sat

# Liquid gains energy equal to what the arriving mass brings
E_liquid_gain = dm_cond * u_liq_sat

# In energy balance:
U_gas_new = U_gas_old - E_gas_loss = U_gas_old - dm_cond * u_vap_sat
U_liquid_new = U_liquid_old + E_liquid_gain = U_liquid_old + dm_cond * u_liq_sat
```

## Proposed Fix

**Lines 2258-2265:** Change evaporation energy definition
```python
# OLD (wrong):
E_evap = dm_evap_mass * (u_vap_sat - u_liq_sat)

# NEW (correct):
E_evap_from_liquid = dm_evap_mass * u_liq_sat  # Energy leaving liquid
E_evap_to_gas = dm_evap_mass * u_vap_sat       # Energy entering gas
```

**Lines 2297-2304:** Change condensation energy definition
```python
# OLD (wrong):
E_cond = dm_cond_mass * (u_liq_sat - u_vap_sat)

# NEW (correct):
E_cond_from_gas = dm_cond_mass * u_vap_sat      # Energy leaving gas
E_cond_to_liquid = dm_cond_mass * u_liq_sat     # Energy entering liquid
```

**Lines 2365, 2378, 2386:** Update energy balance
```python
# Gas phase:
U_gas_new = (U_gas_old
            + ...
            + E_evap_to_gas       # Mass arriving from evaporation
            - E_cond_from_gas)    # Mass leaving to condensation

# Liquid phase:
U_liquid_new = (U_liquid_old
               + ...
               - E_evap_from_liquid   # Mass leaving to evaporation
               + E_cond_to_liquid)    # Mass arriving from condensation
```

## Wait - What About Energy Conservation?

With separate accounting:
```
ΔU_gas = +dm_evap*u_vap - dm_cond*u_vap
ΔU_liquid = -dm_evap*u_liq + dm_cond*u_liq
ΔU_total = dm_evap*(u_vap - u_liq) - dm_cond*(u_vap - u_liq)
         = (dm_evap - dm_cond) * u_fg
```

This is NOT zero unless dm_evap = dm_cond!

**Where does the energy come from/go?**

The answer: **Latent heat must be supplied/removed externally!**

- For evaporation: u_fg energy must come from heat input (Q_inner_wetted)
- For condensation: u_fg energy released as heat (increases Q_gas_liquid)

But in our model, this heat is already accounted for in Q_inner and Q_inner_wetted!

**The issue:** We're double-counting!

If Q_inner_wetted already provides the heat for evaporation, we shouldn't ALSO add u_vap to the gas.

## Alternative: Use Enthalpy (h)

When using enthalpy (h = u + P/ρ):
- Mass leaving carries enthalpy (flow work included)
- Mass entering brings enthalpy
- This is the standard formulation in thermodynamics!

**For evaporation:**
```python
E_evap_from_liquid = dm_evap * h_liq_sat
E_evap_to_gas = dm_evap * h_vap_sat
```

**But then:**
```
ΔU_liquid = -dm*h_liq?  NO! U and H are different!
```

The issue: Our energy balance is in terms of U (internal energy), but mass carries h (enthalpy).

**Conversion:**
```
h = u + Pv = u + P/ρ
H = U + PV

When mass dm leaves at (P, v_liq):
  - Enthalpy leaves: dm * h_liq
  - Internal energy leaves: dm * u_liq
  - PV work done: dm * P*v_liq = dm * (h_liq - u_liq)
```

For an open system (mass can leave), the first law is:
```
dU = δQ - δW + Σ h_in*dm_in - Σ h_out*dm_out
```

NOT:
```
dU = δQ - δW + Σ u_in*dm_in - Σ u_out*dm_out  (WRONG for open systems!)
```

## Conclusion

**For NEM with phase transfer, we MUST use enthalpy (h) for mass transfer terms!**

The energy balance for an open system (with mass transfer) is:
```
dU/dt = Q - W + Σ (h_in * dm_in/dt) - Σ (h_out * dm_out/dt)
```

Where:
- Q = heat transfer
- W = work (for us, W = 0 since constant volume)
- h_in/out = specific enthalpy of incoming/outgoing mass

Using internal energy (u) for mass transfer is fundamentally incorrect for open systems!
