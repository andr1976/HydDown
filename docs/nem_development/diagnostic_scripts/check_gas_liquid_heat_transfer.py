#!/usr/bin/env python
"""Check gas-liquid heat transfer magnitude."""

import yaml
import sys
sys.path.insert(0, 'src')
from hyddown.hdclass import HydDown
import numpy as np

# Load input and enable NEM
with open('src/hyddown/examples/nem_propane_psv_fire.yml') as f:
    input_dict = yaml.load(f, Loader=yaml.FullLoader)
input_dict['calculation']['non_equilibrium'] = True

# Run NEM simulation
hd = HydDown(input_dict)
hd.run()

print("\n" + "="*80)
print("GAS-LIQUID HEAT TRANSFER ANALYSIS")
print("="*80)

# Find max temperature difference
idx_max_dT = np.argmax(np.abs(hd.T_gas - hd.T_liquid))
max_dT = hd.T_gas[idx_max_dT] - hd.T_liquid[idx_max_dT]

print(f"\nMaximum temperature difference:")
print(f"  t = {hd.time_array[idx_max_dT]:.2f} s")
print(f"  T_gas = {hd.T_gas[idx_max_dT]:.2f} K")
print(f"  T_liquid = {hd.T_liquid[idx_max_dT]:.2f} K")
print(f"  ΔT = {max_dT:.2f} K")

# Calculate interface area at this time
if hd.m_liquid[idx_max_dT] > 1e-6 and hd.rho_liquid[idx_max_dT] > 1e-6:
    V_liquid = hd.m_liquid[idx_max_dT] / hd.rho_liquid[idx_max_dT]
    liquid_level = hd.inner_vol.h_from_V(V_liquid)
    # Calculate interface area (cross-sectional area at liquid level)
    # For horizontal cylinder, this is complex - use approximation
    import math
    # Use vessel diameter from inner_vol object
    D = hd.inner_vol.D
    A_interface = math.pi * D**2 / 4  # Cross-sectional area
else:
    A_interface = 0.0
    liquid_level = 0.0

print(f"  Liquid level = {liquid_level:.3f} m")
print(f"  Interface area = {A_interface:.3f} m²")

# Calculate heat transfer rate with different h_gl values
h_gl_values = [50, 100, 200, 500, 1000]

print(f"\n  Heat transfer rate Q_gas_liquid (W):")
for h_gl in h_gl_values:
    Q_gl = h_gl * A_interface * max_dT
    print(f"    h_gl = {h_gl:4d} W/m²K: Q = {Q_gl/1e3:8.2f} kW")

# Compare to fire heat input
print(f"\n  For comparison:")
print(f"    Fire heat flux: ~60 kW/m²")
print(f"    Typical wetted area: ~10 m²")
print(f"    Typical fire heat input: ~600 kW")

# Check actual h_gl used
if "h_gas_liquid" in input_dict['calculation']:
    h_gl_actual = input_dict['calculation']['h_gas_liquid']
else:
    h_gl_actual = 50.0

print(f"\n  Actual h_gl used: {h_gl_actual} W/m²K")
Q_gl_actual = h_gl_actual * A_interface * max_dT
print(f"  Actual Q_gas_liquid: {Q_gl_actual/1e3:.2f} kW")

# Calculate time constant for thermal equilibration
# τ = (m_gas * cp_gas + m_liquid * cp_liquid) / (h_gl * A)
cp_gas = 2000  # J/kgK (approximate for propane gas)
cp_liquid = 2500  # J/kgK (approximate for propane liquid)
thermal_capacity = (hd.m_gas[idx_max_dT] * cp_gas +
                    hd.m_liquid[idx_max_dT] * cp_liquid)

if A_interface > 0:
    time_constant = thermal_capacity / (h_gl_actual * A_interface)
    print(f"\n  Thermal equilibration time constant:")
    print(f"    τ = (m_gas*cp_gas + m_liq*cp_liq) / (h_gl*A)")
    print(f"    τ = {time_constant:.1f} seconds")
    print(f"    At t={hd.time_array[idx_max_dT]:.1f}s, we're {hd.time_array[idx_max_dT]/time_constant:.2f}×τ into simulation")

print("\n" + "="*80)
print("\nCONCLUSION:")
print("  If h_gl = 50 W/m²K is too low, phases won't equilibrate fast enough.")
print("  Large ΔT can persist throughout simulation.")
print("  Consider increasing h_gl to 200-500 W/m²K for faster equilibration.")
print("="*80 + "\n")
