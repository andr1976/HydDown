#!/usr/bin/env python3
"""
Debug NEM energy balance - check if U is actually changing
"""

import yaml
import sys
import numpy as np
sys.path.insert(0, 'src')

from hyddown.hdclass import HydDown

# Monkey-patch to add debug output
original_step = HydDown.step

def debug_step(self):
    i = self.itr
    if i >= 1 and i <= 3 and self.non_equilibrium:  # Only debug first few steps
        print(f"\n=== Timestep {i}: t={self.time_array[i]:.2f}s ===")

        # Get values BEFORE step
        U_gas_before = self.m_gas[i-1] * self.U_gas[i-1]
        U_liquid_before = self.m_liquid[i-1] * self.U_liquid[i-1]

        print(f"BEFORE step:")
        print(f"  U_gas_total:   {U_gas_before/1e6:.3f} MJ")
        print(f"  U_liquid_total: {U_liquid_before/1e6:.3f} MJ")
        print(f"  Q_inner[i]:    {self.Q_inner[i]/1e3:.3f} kW")
        print(f"  Q_wetted[i]:   {self.Q_inner_wetted[i]/1e3:.3f} kW")
        print(f"  dt:            {self.tstep:.2f} s")

        # Expected energy addition
        dU_gas_expected = self.Q_inner[i] * self.tstep
        dU_liquid_expected = self.Q_inner_wetted[i] * self.tstep

        print(f"EXPECTED changes:")
        print(f"  dU_gas:   {dU_gas_expected/1e3:.3f} kJ")
        print(f"  dU_liquid: {dU_liquid_expected/1e3:.3f} kJ")

    # Call original
    result = original_step(self)

    if i >= 1 and i <= 3 and self.non_equilibrium:
        U_gas_after = self.m_gas[i] * self.U_gas[i]
        U_liquid_after = self.m_liquid[i] * self.U_liquid[i]

        print(f"AFTER step:")
        print(f"  U_gas_total:   {U_gas_after/1e6:.3f} MJ")
        print(f"  U_liquid_total: {U_liquid_after/1e6:.3f} MJ")
        print(f"  ACTUAL dU_gas:   {(U_gas_after - U_gas_before)/1e3:.3f} kJ")
        print(f"  ACTUAL dU_liquid: {(U_liquid_after - U_liquid_before)/1e3:.3f} kJ")

    return result

HydDown.step = debug_step

# Load input
with open('src/hyddown/examples/nem_propane_psv_fire.yml') as f:
    input_data = yaml.load(f, Loader=yaml.FullLoader)

# Short simulation
input_data['calculation']['end_time'] = 2.0  # Just 4 steps

# Run NEM
print("Running NEM with energy debug...")
hdown = HydDown(input_data)
hdown.run(disable_pbar=True)

print("\n" + "="*80)
print("SUMMARY")
print("="*80)
print(f"Initial U_gas: {hdown.U_gas[0]:.2f} kJ/kg")
print(f"Final U_gas:   {hdown.U_gas[-1]:.2f} kJ/kg")
print(f"Initial U_liquid: {hdown.U_liquid[0]:.2f} kJ/kg")
print(f"Final U_liquid:   {hdown.U_liquid[-1]:.2f} kJ/kg")
