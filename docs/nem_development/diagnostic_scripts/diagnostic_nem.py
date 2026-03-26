#!/usr/bin/env python3
"""
Diagnostic script to understand NEM energy balance.
"""

import yaml
import sys
import numpy as np
sys.path.insert(0, 'src')

from hyddown.hdclass import HydDown

# Load input
with open('src/hyddown/examples/nem_propane_psv_fire.yml') as f:
    input_data = yaml.load(f, Loader=yaml.FullLoader)

# Run both models
print("="*80)
print("DIAGNOSTIC: Energy Balance Analysis")
print("="*80)

# Equilibrium
input_eq = input_data.copy()
input_eq['calculation'] = input_data['calculation'].copy()
input_eq['calculation']['non_equilibrium'] = False

print("\n[1] Running Equilibrium Model...")
hdown_eq = HydDown(input_eq)
hdown_eq.run(disable_pbar=True)

# NEM
input_nem = input_data.copy()
input_nem['calculation'] = input_data['calculation'].copy()
input_nem['calculation']['non_equilibrium'] = True

print("[2] Running Non-Equilibrium Model...")
hdown_nem = HydDown(input_nem)
hdown_nem.run(disable_pbar=True)

print("\n" + "="*80)
print("ENERGY BALANCE ANALYSIS")
print("="*80)

# Calculate total heat input
dt = hdown_eq.tstep

# Equilibrium
Q_total_eq = np.sum((hdown_eq.Q_inner + hdown_eq.Q_inner_wetted) * dt)
mass_discharged_eq = hdown_eq.mass_fluid[0] - hdown_eq.mass_fluid[-1]

# NEM
Q_total_nem = np.sum((hdown_nem.Q_inner + hdown_nem.Q_inner_wetted) * dt)
mass_discharged_nem = hdown_nem.mass_fluid[0] - hdown_nem.mass_fluid[-1]

print(f"\n{'='*40}")
print(f"Total Heat Input:")
print(f"{'='*40}")
print(f"Equilibrium:      {Q_total_eq/1e6:.2f} MJ")
print(f"Non-Equilibrium:  {Q_total_nem/1e6:.2f} MJ")
print(f"Difference:       {(Q_total_eq - Q_total_nem)/1e6:.2f} MJ ({abs(Q_total_eq - Q_total_nem)/Q_total_eq*100:.1f}%)")

print(f"\n{'='*40}")
print(f"Mass Discharge:")
print(f"{'='*40}")
print(f"Equilibrium:      {mass_discharged_eq:.2f} kg ({mass_discharged_eq/hdown_eq.mass_fluid[0]*100:.1f}%)")
print(f"Non-Equilibrium:  {mass_discharged_nem:.2f} kg ({mass_discharged_nem/hdown_nem.mass_fluid[0]*100:.1f}%)")

print(f"\n{'='*40}")
print(f"Internal Energy Change:")
print(f"{'='*40}")
# For equilibrium
U_initial_eq = hdown_eq.U_mass[0] * hdown_eq.mass_fluid[0]
U_final_eq = hdown_eq.U_mass[-1] * hdown_eq.mass_fluid[-1]
dU_eq = U_final_eq - U_initial_eq
print(f"Equilibrium:      {dU_eq/1e6:.2f} MJ")

# For NEM
U_initial_nem = hdown_nem.U_gas[0] * hdown_nem.m_gas[0] + hdown_nem.U_liquid[0] * hdown_nem.m_liquid[0]
U_final_nem = hdown_nem.U_gas[-1] * hdown_nem.m_gas[-1] + hdown_nem.U_liquid[-1] * hdown_nem.m_liquid[-1]
dU_nem = U_final_nem - U_initial_nem
print(f"Non-Equilibrium:  {dU_nem/1e6:.2f} MJ")

print(f"\n{'='*40}")
print(f"Energy Balance Check (Q_in - Enthalpy_out = dU):")
print(f"{'='*40}")
# For equilibrium: Q_in - H_out*dm_out should equal dU
# Approximate H_out as average enthalpy
H_avg_eq = np.mean(hdown_eq.H_mass)
Energy_out_eq = H_avg_eq * mass_discharged_eq
Balance_eq = Q_total_eq - Energy_out_eq - dU_eq
print(f"Equilibrium:")
print(f"  Q_in:           {Q_total_eq/1e6:.2f} MJ")
print(f"  H_out*dm:       {Energy_out_eq/1e6:.2f} MJ")
print(f"  dU:             {dU_eq/1e6:.2f} MJ")
print(f"  Balance error:  {Balance_eq/1e6:.2f} MJ ({abs(Balance_eq)/Q_total_eq*100:.1f}%)")

# For NEM
H_avg_nem_gas = np.mean([hdown_nem.fluid_gas.hmass() if hdown_nem.m_gas[i] > 0 else 0 for i in range(len(hdown_nem.m_gas))])
Energy_out_nem = H_avg_nem_gas * mass_discharged_nem
Balance_nem = Q_total_nem - Energy_out_nem - dU_nem
print(f"\nNon-Equilibrium:")
print(f"  Q_in:           {Q_total_nem/1e6:.2f} MJ")
print(f"  H_out*dm:       {Energy_out_nem/1e6:.2f} MJ")
print(f"  dU:             {dU_nem/1e6:.2f} MJ")
print(f"  Balance error:  {Balance_nem/1e6:.2f} MJ ({abs(Balance_nem)/Q_total_nem*100:.1f}%)")

print(f"\n{'='*40}")
print(f"Heat Transfer Breakdown:")
print(f"{'='*40}")
# Gas phase heat
Q_gas_eq = np.sum(hdown_eq.Q_inner * dt)
Q_liquid_eq = np.sum(hdown_eq.Q_inner_wetted * dt)

Q_gas_nem = np.sum(hdown_nem.Q_inner * dt)
Q_liquid_nem = np.sum(hdown_nem.Q_inner_wetted * dt)

print(f"Equilibrium:")
print(f"  Q_gas:    {Q_gas_eq/1e6:.2f} MJ ({Q_gas_eq/Q_total_eq*100:.1f}%)")
print(f"  Q_liquid: {Q_liquid_eq/1e6:.2f} MJ ({Q_liquid_eq/Q_total_eq*100:.1f}%)")

print(f"\nNon-Equilibrium:")
print(f"  Q_gas:    {Q_gas_nem/1e6:.2f} MJ ({Q_gas_nem/Q_total_nem*100:.1f}%)")
print(f"  Q_liquid: {Q_liquid_nem/1e6:.2f} MJ ({Q_liquid_nem/Q_total_nem*100:.1f}%)")

print(f"\n{'='*40}")
print(f"Phase Mass Evolution (NEM):")
print(f"{'='*40}")
print(f"Initial:  Gas={hdown_nem.m_gas[0]:.2f} kg, Liquid={hdown_nem.m_liquid[0]:.2f} kg")
print(f"Final:    Gas={hdown_nem.m_gas[-1]:.2f} kg, Liquid={hdown_nem.m_liquid[-1]:.2f} kg")
print(f"Change:   Gas={hdown_nem.m_gas[-1]-hdown_nem.m_gas[0]:.2f} kg, Liquid={hdown_nem.m_liquid[-1]-hdown_nem.m_liquid[0]:.2f} kg")
total_phase_transfer = np.sum(hdown_nem.mdot_phase_transfer * dt)
print(f"Total phase transfer: {total_phase_transfer:.2f} kg")

print("\n" + "="*80)
