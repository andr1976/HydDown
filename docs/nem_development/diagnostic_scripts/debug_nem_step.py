#!/usr/bin/env python3
"""
Debug NEM - print first 5 timesteps to understand what's happening
"""

import yaml
import sys
sys.path.insert(0, 'src')

from hyddown.hdclass import HydDown

# Load input
with open('src/hyddown/examples/nem_propane_psv_fire.yml') as f:
    input_data = yaml.load(f, Loader=yaml.FullLoader)

# Short simulation
input_data['calculation']['end_time'] = 2.5  # Just 5 steps at 0.5s each

# Run NEM
print("Running NEM for 5 timesteps...")
hdown = HydDown(input_data)
hdown.run(disable_pbar=True)

print("\n" + "="*80)
print("TIMESTEP-BY-TIMESTEP ANALYSIS")
print("="*80)

for i in range(min(6, len(hdown.time_array))):
    print(f"\nTimestep {i}: t = {hdown.time_array[i]:.2f} s")
    print(f"  P:          {hdown.P[i]/1e5:.3f} bar")
    print(f"  T_gas:      {hdown.T_gas[i]:.2f} K")
    print(f"  T_liquid:   {hdown.T_liquid[i]:.2f} K")
    print(f"  m_gas:      {hdown.m_gas[i]:.2f} kg")
    print(f"  m_liquid:   {hdown.m_liquid[i]:.2f} kg")
    print(f"  U_gas:      {hdown.U_gas[i]/1e3:.2f} kJ/kg")
    print(f"  U_liquid:   {hdown.U_liquid[i]/1e3:.2f} kJ/kg")
    print(f"  Q_inner:    {hdown.Q_inner[i]/1e3:.2f} kW (dry wall)")
    print(f"  Q_wetted:   {hdown.Q_inner_wetted[i]/1e3:.2f} kW (wetted wall)")
    print(f"  mdot:       {hdown.mass_rate[i]:.4f} kg/s")
    print(f"  mdot_phase: {hdown.mdot_phase_transfer[i]:.4e} kg/s")

print("\n" + "="*80)
