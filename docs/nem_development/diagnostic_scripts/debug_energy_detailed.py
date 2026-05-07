#!/usr/bin/env python3
"""
Detailed energy conservation analysis - check phase transfer energy
"""

import yaml
import sys
import numpy as np
sys.path.insert(0, 'src')

from hyddown.hdclass import HydDown

# Load input
with open('src/hyddown/examples/nem_propane_psv_fire.yml') as f:
    input_data = yaml.load(f, Loader=yaml.FullLoader)

# Short simulation
input_data['calculation']['end_time'] = 10.0

print("="*80)
print("DETAILED ENERGY CONSERVATION ANALYSIS")
print("="*80)

# Run NEM
hdown = HydDown(input_data)
hdown.run(disable_pbar=True)

dt = hdown.tstep

# Analyze each timestep
print("\n" + "="*80)
print("TIMESTEP-BY-TIMESTEP ENERGY ACCOUNTING")
print("="*80)

for i in range(1, min(21, len(hdown.time_array))):
    # Previous state
    U_gas_prev = hdown.m_gas[i-1] * hdown.U_gas[i-1]
    U_liquid_prev = hdown.m_liquid[i-1] * hdown.U_liquid[i-1]
    U_total_prev = U_gas_prev + U_liquid_prev

    # Current state
    U_gas_curr = hdown.m_gas[i] * hdown.U_gas[i]
    U_liquid_curr = hdown.m_liquid[i] * hdown.U_liquid[i]
    U_total_curr = U_gas_curr + U_liquid_curr

    # Energy changes
    dU_gas = U_gas_curr - U_gas_prev
    dU_liquid = U_liquid_curr - U_liquid_prev
    dU_total = U_total_curr - U_total_prev

    # Heat input
    Q_gas = hdown.Q_inner[i] * dt
    Q_liquid = hdown.Q_inner_wetted[i] * dt
    Q_total = Q_gas + Q_liquid

    # Phase transfer
    dm_phase = hdown.m_gas[i] - hdown.m_gas[i-1]  # Positive = evaporation

    # Mass discharge
    dm_discharge = hdown.mass_fluid[i-1] - hdown.mass_fluid[i]

    # Mass flow energy (enthalpy out)
    if dm_discharge > 0:
        # Estimate discharge enthalpy
        h_discharge = hdown.U_gas[i-1] + hdown.P[i-1] / hdown.rho_gas[i-1] if hdown.rho_gas[i-1] > 0 else 0
        E_discharge = dm_discharge * h_discharge
    else:
        E_discharge = 0.0

    # Expected energy change
    dU_expected = Q_total - E_discharge

    # Residual
    residual = dU_total - dU_expected

    if i % 1 == 0:  # Print every step for first 20
        print(f"\nt = {hdown.time_array[i]:.1f} s:")
        print(f"  Heat input:")
        print(f"    Q_gas:     {Q_gas/1e3:8.3f} kJ")
        print(f"    Q_liquid:  {Q_liquid/1e3:8.3f} kJ")
        print(f"    Q_total:   {Q_total/1e3:8.3f} kJ")

        print(f"  Internal energy change:")
        print(f"    dU_gas:    {dU_gas/1e3:8.3f} kJ")
        print(f"    dU_liquid: {dU_liquid/1e3:8.3f} kJ")
        print(f"    dU_total:  {dU_total/1e3:8.3f} kJ")

        print(f"  Phase transfer:")
        print(f"    dm_phase:  {dm_phase:8.4f} kg {'(evap)' if dm_phase > 0 else '(cond)'}")

        print(f"  Discharge:")
        print(f"    dm_out:    {dm_discharge:8.4f} kg")
        print(f"    E_out:     {E_discharge/1e3:8.3f} kJ")

        print(f"  Energy balance:")
        print(f"    Expected:  {dU_expected/1e3:8.3f} kJ")
        print(f"    Actual:    {dU_total/1e3:8.3f} kJ")
        print(f"    Residual:  {residual/1e3:8.3f} kJ ({abs(residual)/abs(dU_expected)*100:5.1f}%)")

        # Check individual phase energy balances
        print(f"  Phase-specific:")
        print(f"    U_gas:     {U_gas_curr/1e6:8.4f} MJ (rho={hdown.rho_gas[i]:6.2f} kg/m³, T={hdown.T_gas[i]:6.2f} K)")
        print(f"    U_liquid:  {U_liquid_curr/1e6:8.4f} MJ (rho={hdown.rho_liquid[i]:6.2f} kg/m³, T={hdown.T_liquid[i]:6.2f} K)")
        print(f"    P:         {hdown.P[i]/1e5:8.3f} bar")

print("\n" + "="*80)
