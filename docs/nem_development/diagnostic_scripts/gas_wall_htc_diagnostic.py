#!/usr/bin/env python
"""Diagnostic for gas/wall heat transfer coefficient in NEM."""

import sys
sys.path.insert(0, 'src')
from hyddown.hdclass import HydDown
import yaml
import numpy as np

print("\n" + "="*80)
print("GAS/WALL HEAT TRANSFER DIAGNOSTIC FOR NEM")
print("="*80)

# Run NEM simulation
with open('src/hyddown/examples/nem_propane_psv_fire.yml') as f:
    input_nem = yaml.load(f, Loader=yaml.FullLoader)

input_nem['calculation']['h_gas_liquid'] = 500
input_nem['calculation']['end_time'] = 400.0

print("\nRunning NEM simulation (h_gl=500)...")
hd = HydDown(input_nem)
hd.run()

print(f"\n" + "="*80)
print("ANALYSIS OF GAS-SIDE (UNWETTED) HEAT TRANSFER")
print("="*80)

times = [200, 300, 400]

for t in times:
    idx = np.argmin(np.abs(hd.time_array - t))

    print(f"\nAt t = {hd.time_array[idx]:.1f} s:")
    print(f"  Phase Temperatures:")
    print(f"    T_gas:    {hd.T_gas[idx]:.2f} K ({hd.T_gas[idx]-273.15:.1f}°C)")
    print(f"    T_liquid: {hd.T_liquid[idx]:.2f} K ({hd.T_liquid[idx]-273.15:.1f}°C)")
    print(f"    T_fluid (equilibrium): {hd.T_fluid[idx]:.2f} K ({hd.T_fluid[idx]-273.15:.1f}°C)")
    print(f"    ΔT (gas-liquid): {hd.T_gas[idx] - hd.T_liquid[idx]:.2f} K")

    print(f"\n  Wall Temperatures:")
    print(f"    T_wall (unwetted): {hd.T_vessel[idx]:.2f} K ({hd.T_vessel[idx]-273.15:.1f}°C)")
    print(f"    T_wall (wetted):   {hd.T_vessel_wetted[idx]:.2f} K ({hd.T_vessel_wetted[idx]-273.15:.1f}°C)")

    print(f"\n  Heat Transfer Coefficients:")
    print(f"    h_inside (unwetted/gas-side): {hd.h_inside[idx]:.1f} W/m²K")
    print(f"    h_inside_wetted (liquid-side): {hd.h_inside_wetted[idx]:.1f} W/m²K")

    print(f"\n  Temperature Driving Forces:")
    print(f"    ΔT_gas (wall-gas):    {hd.T_vessel[idx] - hd.T_gas[idx]:.2f} K (SHOULD use this)")
    print(f"    ΔT_fluid (wall-fluid): {hd.T_vessel[idx] - hd.T_fluid[idx]:.2f} K (CURRENTLY using this)")
    print(f"    Error in ΔT: {(hd.T_vessel[idx] - hd.T_fluid[idx]) - (hd.T_vessel[idx] - hd.T_gas[idx]):.2f} K")

    print(f"\n  Heat Fluxes:")
    print(f"    q_inner (unwetted): {hd.q_inner[idx]/1000:.1f} kW/m²")
    print(f"    q_inner_wetted:     {hd.q_inner_wetted[idx]/1000:.1f} kW/m²")

    # Calculate what it SHOULD be
    q_inner_correct = hd.h_inside[idx] * (hd.T_vessel[idx] - hd.T_gas[idx])
    q_inner_current = hd.q_inner[idx]

    print(f"\n  Heat Flux Analysis:")
    print(f"    Current (using T_fluid):  {q_inner_current/1000:.2f} kW/m²")
    print(f"    Correct (using T_gas):    {q_inner_correct/1000:.2f} kW/m²")
    print(f"    Error: {(q_inner_current - q_inner_correct)/1000:.2f} kW/m² ({100*(q_inner_current - q_inner_correct)/q_inner_correct:.1f}%)")

print(f"\n" + "="*80)
print("ISSUE IDENTIFIED")
print("="*80)
print("""
For NEM, the UNWETTED wall is in contact with the GAS PHASE, but the code uses:
- T_fluid (equilibrium temperature) instead of T_gas
- equilibrium fluid properties instead of gas phase properties

This causes INCORRECT heat transfer calculation on the unwetted (gas-side) wall!

LOCATIONS TO FIX in hdclass.py:
1. Line ~1230: T_film calculation for h_inside_mixed
2. Line ~1304: T_film calculation for h_inside
3. Line ~1805: T_fluid passed to tp.h_inner()
4. Line ~1410, 1868: Q_inner calculation using T_fluid
5. Line ~1416, 1874: q_inner calculation using T_fluid

Should use T_gas for NEM, T_fluid for equilibrium.
""")
print("="*80)
