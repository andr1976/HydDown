#!/usr/bin/env python3
"""
Check conservation of mass, volume, and energy in NEM simulation
"""

import yaml
import sys
import numpy as np
sys.path.insert(0, 'src')

from hyddown.hdclass import HydDown

# Load input
with open('src/hyddown/examples/nem_propane_psv_fire.yml') as f:
    input_data = yaml.load(f, Loader=yaml.FullLoader)

# Short simulation for debugging
input_data['calculation']['end_time'] = 50.0  # 100 steps

print("="*80)
print("CONSERVATION CHECK: NEM Simulation")
print("="*80)

# Run NEM
hdown = HydDown(input_data)
hdown.run(disable_pbar=True)

print("\n" + "="*80)
print("CHECKING CONSERVATION LAWS")
print("="*80)

# Check every 10 timesteps
check_indices = range(0, len(hdown.time_array), 10)

print("\nMASS CONSERVATION:")
print(f"{'Time [s]':<10} {'m_gas [kg]':<12} {'m_liquid [kg]':<12} {'m_total [kg]':<12} {'V_gas [m³]':<12} {'V_liquid [m³]':<12} {'V_total [m³]':<12}")
print("-" * 100)

for i in check_indices:
    m_total = hdown.m_gas[i] + hdown.m_liquid[i]

    # Calculate volumes
    if hdown.rho_gas[i] > 0:
        V_gas = hdown.m_gas[i] / hdown.rho_gas[i]
    else:
        V_gas = 0.0

    if hdown.rho_liquid[i] > 0:
        V_liquid = hdown.m_liquid[i] / hdown.rho_liquid[i]
    else:
        V_liquid = 0.0

    V_total = V_gas + V_liquid

    print(f"{hdown.time_array[i]:<10.1f} {hdown.m_gas[i]:<12.2f} {hdown.m_liquid[i]:<12.2f} {m_total:<12.2f} {V_gas:<12.4f} {V_liquid:<12.4f} {V_total:<12.4f}")

print(f"\nVessel volume: {hdown.vol:.4f} m³")
print(f"Initial mass:  {hdown.mass_fluid[0]:.2f} kg")
print(f"Final mass:    {hdown.mass_fluid[-1]:.2f} kg")
print(f"Mass discharged: {hdown.mass_fluid[0] - hdown.mass_fluid[-1]:.2f} kg")

# Energy conservation
print("\n" + "="*80)
print("ENERGY CONSERVATION:")
print("="*80)

dt = hdown.tstep

for i in check_indices[1:]:  # Skip initial
    # Energy at this timestep
    U_gas = hdown.m_gas[i] * hdown.U_gas[i]
    U_liquid = hdown.m_liquid[i] * hdown.U_liquid[i]
    U_total = U_gas + U_liquid

    # Energy at previous timestep
    U_gas_prev = hdown.m_gas[i-1] * hdown.U_gas[i-1]
    U_liquid_prev = hdown.m_liquid[i-1] * hdown.U_liquid[i-1]
    U_total_prev = U_gas_prev + U_liquid_prev

    # Energy change
    dU = U_total - U_total_prev

    # Heat input
    Q_in = (hdown.Q_inner[i] + hdown.Q_inner_wetted[i]) * dt

    # Mass flow energy (enthalpy carried out)
    # Note: For discharge, we lose enthalpy
    dm = hdown.mass_fluid[i-1] - hdown.mass_fluid[i]
    if dm > 0:
        # Estimate discharge enthalpy (use gas properties as rough estimate)
        h_discharge = hdown.U_gas[i-1] + hdown.P[i-1] / hdown.rho_gas[i-1] if hdown.rho_gas[i-1] > 0 else 0
        E_out = dm * h_discharge
    else:
        E_out = 0.0

    # Energy balance: dU = Q_in - E_out
    expected_dU = Q_in - E_out
    residual = dU - expected_dU

    if i % 10 == 0:
        print(f"\nt = {hdown.time_array[i]:.1f} s:")
        print(f"  dU (actual):   {dU/1e3:.3f} kJ")
        print(f"  Q_in:          {Q_in/1e3:.3f} kJ")
        print(f"  E_out:         {E_out/1e3:.3f} kJ")
        print(f"  Expected dU:   {expected_dU/1e3:.3f} kJ")
        print(f"  Residual:      {residual/1e3:.3f} kJ ({abs(residual)/abs(expected_dU)*100:.1f}%)")
        print(f"  P:             {hdown.P[i]/1e5:.2f} bar")
        print(f"  T_gas:         {hdown.T_gas[i]:.1f} K")
        print(f"  T_liquid:      {hdown.T_liquid[i]:.1f} K")

print("\n" + "="*80)
print("VOLUME CONSTRAINT CHECK:")
print("="*80)

for i in check_indices:
    if hdown.rho_gas[i] > 0:
        V_gas = hdown.m_gas[i] / hdown.rho_gas[i]
    else:
        V_gas = 0.0

    if hdown.rho_liquid[i] > 0:
        V_liquid = hdown.m_liquid[i] / hdown.rho_liquid[i]
    else:
        V_liquid = 0.0

    V_total = V_gas + V_liquid
    V_error = V_total - hdown.vol

    print(f"t = {hdown.time_array[i]:6.1f} s: V_total = {V_total:.6f} m³, V_vessel = {hdown.vol:.6f} m³, Error = {V_error:.6f} m³ ({abs(V_error)/hdown.vol*100:.3f}%)")

print("\n" + "="*80)
