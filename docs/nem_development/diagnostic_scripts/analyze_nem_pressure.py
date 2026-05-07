#!/usr/bin/env python
"""Analyze NEM pressure spike."""

import matplotlib.pyplot as plt
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
hd_nem = HydDown(input_dict)
hd_nem.run()

# Find when PSV opens
psv_set = 14.3e5  # Pa
psv_open_idx = np.where(hd_nem.P > psv_set)[0]

if len(psv_open_idx) > 0:
    idx_open = psv_open_idx[0]
    t_open = hd_nem.time_array[idx_open]

    print(f"\n{'='*80}")
    print("PSV OPENING ANALYSIS")
    print(f"{'='*80}")
    print(f"PSV opens at t = {t_open:.2f} s (index {idx_open})")
    print(f"  P = {hd_nem.P[idx_open]/1e5:.2f} bar")
    print(f"  T_gas = {hd_nem.T_gas[idx_open]:.2f} K")
    print(f"  T_liquid = {hd_nem.T_liquid[idx_open]:.2f} K")
    print(f"  m_gas = {hd_nem.m_gas[idx_open]:.2f} kg")
    print(f"  m_liquid = {hd_nem.m_liquid[idx_open]:.2f} kg")
    print(f"  mass_rate = {hd_nem.mass_rate[idx_open]:.3f} kg/s")

    # Check a few timesteps after opening
    print(f"\nTimesteps after PSV opens:")
    for j in range(5):
        idx = idx_open + j
        if idx < len(hd_nem.time_array):
            print(f"  t={hd_nem.time_array[idx]:6.2f}s: P={hd_nem.P[idx]/1e5:6.2f} bar, "
                  f"mdot={hd_nem.mass_rate[idx]:7.3f} kg/s, "
                  f"T_gas={hd_nem.T_gas[idx]:6.1f}K, "
                  f"T_liq={hd_nem.T_liquid[idx]:6.1f}K")

# Find maximum pressure
idx_max_p = np.argmax(hd_nem.P)
print(f"\n{'='*80}")
print("MAXIMUM PRESSURE")
print(f"{'='*80}")
print(f"Max pressure at t = {hd_nem.time_array[idx_max_p]:.2f} s")
print(f"  P = {hd_nem.P[idx_max_p]/1e5:.2f} bar")
print(f"  T_gas = {hd_nem.T_gas[idx_max_p]:.2f} K")
print(f"  T_liquid = {hd_nem.T_liquid[idx_max_p]:.2f} K")
print(f"  ΔT = {hd_nem.T_gas[idx_max_p] - hd_nem.T_liquid[idx_max_p]:.2f} K")
print(f"  m_gas = {hd_nem.m_gas[idx_max_p]:.2f} kg")
print(f"  m_liquid = {hd_nem.m_liquid[idx_max_p]:.2f} kg")
print(f"  mass_rate = {hd_nem.mass_rate[idx_max_p]:.3f} kg/s")

# Check phase transfer rate
print(f"\n{'='*80}")
print("PHASE TRANSFER RATES")
print(f"{'='*80}")
max_phase_transfer = np.max(np.abs(hd_nem.mdot_phase_transfer))
avg_phase_transfer = np.mean(np.abs(hd_nem.mdot_phase_transfer[1:]))
print(f"Max phase transfer rate: {max_phase_transfer:.4f} kg/s")
print(f"Avg phase transfer rate: {avg_phase_transfer:.4f} kg/s")

# Find when phase transfer is active
active_transfer = np.abs(hd_nem.mdot_phase_transfer) > 1e-6
n_active = np.sum(active_transfer)
print(f"Phase transfer active in {n_active}/{len(hd_nem.time_array)} timesteps ({n_active/len(hd_nem.time_array)*100:.1f}%)")

# Plot pressure and mass flow rate around PSV opening
if len(psv_open_idx) > 0:
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(10, 8), sharex=True)

    # Pressure
    ax1.plot(hd_nem.time_array, hd_nem.P/1e5, 'b-', label='Pressure')
    ax1.axhline(psv_set/1e5, color='r', linestyle='--', label='PSV set')
    ax1.axvline(t_open, color='gray', linestyle=':', alpha=0.5)
    ax1.set_ylabel('Pressure [bar]')
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    # Mass flow rate
    ax2.plot(hd_nem.time_array, hd_nem.mass_rate, 'g-', label='Discharge rate')
    ax2.axvline(t_open, color='gray', linestyle=':', alpha=0.5)
    ax2.set_ylabel('Mass flow [kg/s]')
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    # Temperature
    ax3.plot(hd_nem.time_array, hd_nem.T_gas, 'r-', label='T_gas')
    ax3.plot(hd_nem.time_array, hd_nem.T_liquid, 'b-', label='T_liquid')
    ax3.axvline(t_open, color='gray', linestyle=':', alpha=0.5)
    ax3.set_ylabel('Temperature [K]')
    ax3.set_xlabel('Time [s]')
    ax3.legend()
    ax3.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('nem_pressure_analysis.pdf')
    print(f"\n✓ Analysis plot saved to 'nem_pressure_analysis.pdf'")

print(f"\n{'='*80}")
