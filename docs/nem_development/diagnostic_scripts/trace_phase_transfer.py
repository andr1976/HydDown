#!/usr/bin/env python3
"""
Trace phase transfer in NEM model.
"""

import yaml
import sys
import numpy as np
import CoolProp.CoolProp as CP
sys.path.insert(0, 'src')

from hyddown.hdclass import HydDown

# Load and run NEM
with open('src/hyddown/examples/nem_propane_psv_fire.yml') as f:
    input_data = yaml.load(f, Loader=yaml.FullLoader)

input_data['calculation']['non_equilibrium'] = True
input_data['calculation']['end_time'] = 60.0  # Short run

print("="*80)
print("PHASE TRANSFER TRACING")
print("="*80)

hdown = HydDown(input_data)

print(f"\nInitial State:")
print(f"  m_gas:    {hdown.m_gas0:.3f} kg")
print(f"  m_liquid: {hdown.m_liquid0:.3f} kg")
print(f"  T_gas:    {hdown.T_gas0:.2f} K")
print(f"  T_liquid: {hdown.T_liquid0:.2f} K")
print(f"  P:        {hdown.p0/1e5:.3f} bar")

hdown.run(disable_pbar=True)

print(f"\nFinal State:")
print(f"  m_gas:    {hdown.m_gas[-1]:.3f} kg")
print(f"  m_liquid: {hdown.m_liquid[-1]:.3f} kg")
print(f"  T_gas:    {hdown.T_gas[-1]:.2f} K")
print(f"  T_liquid: {hdown.T_liquid[-1]:.2f} K")
print(f"  P:        {hdown.P[-1]/1e5:.3f} bar")

print(f"\nPhase Mass Changes:")
print(f"  Δm_gas:    {hdown.m_gas[-1] - hdown.m_gas[0]:.3f} kg")
print(f"  Δm_liquid: {hdown.m_liquid[-1] - hdown.m_liquid[0]:.3f} kg")
print(f"  Δm_total:  {(hdown.m_gas[-1] + hdown.m_liquid[-1]) - (hdown.m_gas[0] + hdown.m_liquid[0]):.3f} kg")

# Check phase transfer details
total_evap = np.sum(hdown.mdot_phase_transfer[hdown.mdot_phase_transfer < 0] * hdown.tstep)
total_cond = np.sum(hdown.mdot_phase_transfer[hdown.mdot_phase_transfer > 0] * hdown.tstep)
net_transfer = np.sum(hdown.mdot_phase_transfer * hdown.tstep)

print(f"\nPhase Transfer Summary:")
print(f"  Total evaporation:   {abs(total_evap):.3f} kg (liquid → gas)")
print(f"  Total condensation:  {total_cond:.3f} kg (gas → liquid)")
print(f"  Net transfer:        {net_transfer:.3f} kg")
print(f"  Convention: dm_transfer > 0 means condensation (gas → liquid)")
print(f"              dm_transfer < 0 means evaporation (liquid → gas)")

# Sample some timesteps
print(f"\nSample Timesteps:")
print(f"{'Time [s]':<10} {'T_gas [K]':<12} {'T_liq [K]':<12} {'T_sat [K]':<12} {'dm/dt [kg/s]':<15} {'Action':<15}")
print("-" * 90)

for idx in [0, 10, 20, 40, 60, 80, 100, 119]:
    if idx >= len(hdown.time_array):
        break
    t = hdown.time_array[idx]
    T_g = hdown.T_gas[idx]
    T_l = hdown.T_liquid[idx]
    P = hdown.P[idx]

    # Get saturation temperature
    hdown.fluid.update(CP.PQ_INPUTS, P, 0.5)
    T_sat = hdown.fluid.T()

    dm_dt = hdown.mdot_phase_transfer[idx]

    if abs(dm_dt) < 1e-6:
        action = "No transfer"
    elif dm_dt > 0:
        action = f"Cond: {dm_dt:.3e}"
    else:
        action = f"Evap: {abs(dm_dt):.3e}"

    print(f"{t:<10.1f} {T_g:<12.2f} {T_l:<12.2f} {T_sat:<12.2f} {dm_dt:<15.3e} {action:<15}")

print("\n" + "="*80)
