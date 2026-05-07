#!/usr/bin/env python
"""Diagnostic for PSV opening behavior - wetted heat flux and gas temperature."""

import sys
sys.path.insert(0, 'src')
from hyddown.hdclass import HydDown
import yaml
import numpy as np
import matplotlib.pyplot as plt

print("\n" + "="*80)
print("PSV OPENING DIAGNOSTIC")
print("="*80)

# Run NEM simulation
with open('src/hyddown/examples/nem_propane_psv_fire.yml') as f:
    input_nem = yaml.load(f, Loader=yaml.FullLoader)

print("\nRunning NEM simulation...")
hd = HydDown(input_nem)
hd.run()

# Find when PSV opens
mass_rate_threshold = 0.1  # kg/s
psv_open = np.abs(hd.mass_rate) > mass_rate_threshold
psv_open_times = hd.time_array[psv_open]

if len(psv_open_times) > 0:
    first_opening = psv_open_times[0]
    print(f"\n✓ PSV first opens at t = {first_opening:.1f} s")
else:
    print("\n! PSV never opens")
    first_opening = hd.time_array[-1]

# Create detailed time series plot
fig, axes = plt.subplots(4, 2, figsize=(16, 14))

# Subplot 1: Pressure and PSV activity
ax1 = axes[0, 0]
ax1.plot(hd.time_array, hd.P/1e5, 'b-', linewidth=2, label='Pressure')
ax1.axhline(14.3, color='r', linestyle='--', alpha=0.5, label='PSV set')
ax1_twin = ax1.twinx()
ax1_twin.plot(hd.time_array, np.abs(hd.mass_rate), 'r-', linewidth=1, alpha=0.7, label='Mass rate')
ax1.set_xlabel('Time [s]')
ax1.set_ylabel('Pressure [bar]', color='b')
ax1_twin.set_ylabel('Mass rate [kg/s]', color='r')
ax1.set_title('Pressure and PSV Activity')
ax1.legend(loc='upper left')
ax1_twin.legend(loc='upper right')
ax1.grid(True, alpha=0.3)

# Subplot 2: Temperatures
ax2 = axes[0, 1]
ax2.plot(hd.time_array, hd.T_gas - 273.15, 'r-', linewidth=2, label='T_gas')
ax2.plot(hd.time_array, hd.T_liquid - 273.15, 'b-', linewidth=2, label='T_liquid')
ax2.plot(hd.time_array, hd.T_vessel_wetted - 273.15, 'g--', linewidth=2, label='T_wall_wetted')
ax2.plot(hd.time_array, hd.T_vessel - 273.15, 'orange', linestyle='-.', linewidth=1, label='T_wall_unwetted')
ax2.set_xlabel('Time [s]')
ax2.set_ylabel('Temperature [°C]')
ax2.set_title('Phase and Wall Temperatures')
ax2.legend()
ax2.grid(True, alpha=0.3)

# Subplot 3: Wetted heat flux and temperature difference
ax3 = axes[1, 0]
ax3.plot(hd.time_array, hd.q_inner_wetted/1000, 'b-', linewidth=2, label='q_wetted')
ax3_twin = ax3.twinx()
dT_wetted = hd.T_vessel_wetted - hd.T_liquid
ax3_twin.plot(hd.time_array, dT_wetted, 'r-', linewidth=1, alpha=0.7, label='ΔT (wall-liq)')
ax3.set_xlabel('Time [s]')
ax3.set_ylabel('q_wetted [kW/m²]', color='b')
ax3_twin.set_ylabel('ΔT [K]', color='r')
ax3.set_title('Wetted Heat Flux and Temperature Difference')
ax3.legend(loc='upper left')
ax3_twin.legend(loc='upper right')
ax3.grid(True, alpha=0.3)

# Subplot 4: Wetted heat transfer coefficient
ax4 = axes[1, 1]
ax4.plot(hd.time_array, hd.h_inside_wetted, 'b-', linewidth=2)
ax4.set_xlabel('Time [s]')
ax4.set_ylabel('h_inside_wetted [W/m²K]')
ax4.set_title('Wetted Heat Transfer Coefficient')
ax4.grid(True, alpha=0.3)

# Subplot 5: Liquid mass and level
ax5 = axes[2, 0]
ax5.plot(hd.time_array, hd.m_liquid, 'b-', linewidth=2, label='Liquid mass')
ax5_twin = ax5.twinx()
ax5_twin.plot(hd.time_array, hd.liquid_level, 'r-', linewidth=1, alpha=0.7, label='Liquid level')
ax5.set_xlabel('Time [s]')
ax5.set_ylabel('Liquid mass [kg]', color='b')
ax5_twin.set_ylabel('Liquid level [m]', color='r')
ax5.set_title('Liquid Inventory')
ax5.legend(loc='upper left')
ax5_twin.legend(loc='upper right')
ax5.grid(True, alpha=0.3)

# Subplot 6: Gas-liquid heat transfer
ax6 = axes[2, 1]
if hasattr(hd, 'Q_gas_liquid'):
    ax6.plot(hd.time_array, hd.Q_gas_liquid/1000, 'purple', linewidth=2)
    ax6.set_xlabel('Time [s]')
    ax6.set_ylabel('Q_gas_liquid [kW]')
    ax6.set_title('Gas-Liquid Heat Transfer')
    ax6.grid(True, alpha=0.3)
else:
    ax6.text(0.5, 0.5, 'Q_gas_liquid not available', ha='center', va='center')
    ax6.axis('off')

# Subplot 7: Total heat transfers
ax7 = axes[3, 0]
ax7.plot(hd.time_array, hd.Q_inner/1000, 'r-', linewidth=2, label='Q_inner (gas-side)')
ax7.plot(hd.time_array, hd.Q_inner_wetted/1000, 'b-', linewidth=2, label='Q_inner_wetted (liq-side)')
ax7.plot(hd.time_array, hd.Q_outer/1000, 'orange', linewidth=1, alpha=0.7, label='Q_outer')
ax7.set_xlabel('Time [s]')
ax7.set_ylabel('Heat Transfer [kW]')
ax7.set_title('Total Heat Transfers')
ax7.legend()
ax7.grid(True, alpha=0.3)

# Subplot 8: Correlation between T_gas and T_wall_wetted
ax8 = axes[3, 1]
ax8.plot(hd.T_vessel_wetted - 273.15, hd.T_gas - 273.15, 'bo', markersize=2, alpha=0.5)
ax8.set_xlabel('T_wall_wetted [°C]')
ax8.set_ylabel('T_gas [°C]')
ax8.set_title('Gas Temperature vs Wetted Wall Temperature')
ax8.grid(True, alpha=0.3)
# Add diagonal line for reference
min_T = min(np.min(hd.T_vessel_wetted - 273.15), np.min(hd.T_gas - 273.15))
max_T = max(np.max(hd.T_vessel_wetted - 273.15), np.max(hd.T_gas - 273.15))
ax8.plot([min_T, max_T], [min_T, max_T], 'r--', alpha=0.5, label='T_gas = T_wall')
ax8.legend()

plt.tight_layout()
plt.savefig('psv_opening_diagnostic.pdf', dpi=300, bbox_inches='tight')
plt.savefig('psv_opening_diagnostic.png', dpi=150, bbox_inches='tight')

print("\n✓ Plots saved: psv_opening_diagnostic.pdf/png")

# Detailed analysis at key times
print("\n" + "="*80)
print("DETAILED ANALYSIS")
print("="*80)

# Find times when q_wetted drops significantly
q_wetted_avg = np.mean(hd.q_inner_wetted[hd.time_array > 200])
q_wetted_drops = hd.q_inner_wetted < 0.1 * q_wetted_avg
drop_times = hd.time_array[q_wetted_drops]

if len(drop_times) > 0:
    print(f"\nTimes when q_wetted drops below 10% of average:")
    for t in drop_times[::50]:  # Sample every 50th point
        idx = np.argmin(np.abs(hd.time_array - t))
        print(f"\n  t = {t:.1f} s:")
        print(f"    P: {hd.P[idx]/1e5:.2f} bar")
        print(f"    Mass rate: {hd.mass_rate[idx]:.2f} kg/s")
        print(f"    T_gas: {hd.T_gas[idx]:.2f} K")
        print(f"    T_liquid: {hd.T_liquid[idx]:.2f} K")
        print(f"    T_wall_wetted: {hd.T_vessel_wetted[idx]:.2f} K")
        print(f"    ΔT (wall-liq): {hd.T_vessel_wetted[idx] - hd.T_liquid[idx]:.2f} K")
        print(f"    h_wetted: {hd.h_inside_wetted[idx]:.1f} W/m²K")
        print(f"    q_wetted: {hd.q_inner_wetted[idx]/1000:.2f} kW/m²")
        print(f"    m_liquid: {hd.m_liquid[idx]:.1f} kg")
        print(f"    liquid_level: {hd.liquid_level[idx]:.4f} m")

# Analysis of T_gas following T_wall_wetted
print(f"\n" + "="*80)
print("CORRELATION ANALYSIS: T_gas vs T_wall_wetted")
print("="*80)

# Calculate correlation coefficient
correlation = np.corrcoef(hd.T_gas, hd.T_vessel_wetted)[0, 1]
print(f"\nCorrelation coefficient: {correlation:.4f}")

# Check if T_gas is close to T_wall_wetted
diff_gas_wall = hd.T_gas - hd.T_vessel_wetted
avg_diff = np.mean(diff_gas_wall)
std_diff = np.std(diff_gas_wall)

print(f"\nT_gas - T_wall_wetted:")
print(f"  Mean: {avg_diff:.2f} K")
print(f"  Std dev: {std_diff:.2f} K")
print(f"  Min: {np.min(diff_gas_wall):.2f} K")
print(f"  Max: {np.max(diff_gas_wall):.2f} K")

# Look at gas-side heat transfer
print(f"\n" + "="*80)
print("GAS-SIDE HEAT TRANSFER ANALYSIS")
print("="*80)

idx_400 = np.argmin(np.abs(hd.time_array - 400))
print(f"\nAt t = 400s:")
print(f"  T_gas: {hd.T_gas[idx_400]:.2f} K ({hd.T_gas[idx_400]-273.15:.1f}°C)")
print(f"  T_wall (unwetted): {hd.T_vessel[idx_400]:.2f} K ({hd.T_vessel[idx_400]-273.15:.1f}°C)")
print(f"  T_wall (wetted): {hd.T_vessel_wetted[idx_400]:.2f} K ({hd.T_vessel_wetted[idx_400]-273.15:.1f}°C)")
print(f"  ΔT (wall_unwet - gas): {hd.T_vessel[idx_400] - hd.T_gas[idx_400]:.2f} K")
print(f"  h_inside (gas-side): {hd.h_inside[idx_400]:.1f} W/m²K")
print(f"  q_inner (gas-side): {hd.q_inner[idx_400]/1000:.1f} kW/m²")

print("\n" + "="*80)
