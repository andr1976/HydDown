#!/usr/bin/env python
"""Compare wetted internal heat transfer coefficient between equilibrium and NEM."""

import sys
sys.path.insert(0, 'src')
from hyddown.hdclass import HydDown
import yaml
import numpy as np
import matplotlib.pyplot as plt

print("\n" + "="*80)
print("WETTED HEAT TRANSFER COEFFICIENT: EQUILIBRIUM vs NEM (h_gl=500)")
print("="*80)

# ============================================================================
# Run Equilibrium Simulation
# ============================================================================
with open('src/hyddown/examples/nem_propane_psv_fire.yml') as f:
    input_eq = yaml.load(f, Loader=yaml.FullLoader)

input_eq['calculation']['non_equilibrium'] = False
input_eq['calculation']['end_time'] = 400.0
input_eq['calculation']['h_gas_liquid'] = 500  # Not used in equilibrium, but set for consistency

print("\nRunning EQUILIBRIUM simulation...")
hd_eq = HydDown(input_eq)
hd_eq.run()
print("✓ Equilibrium complete")

# ============================================================================
# Run NEM Simulation with h_gl = 500
# ============================================================================
with open('src/hyddown/examples/nem_propane_psv_fire.yml') as f:
    input_nem = yaml.load(f, Loader=yaml.FullLoader)

input_nem['calculation']['non_equilibrium'] = True
input_nem['calculation']['end_time'] = 400.0
input_nem['calculation']['h_gas_liquid'] = 500  # Moderate coupling

print("Running NEM simulation (h_gl=500)...")
hd_nem = HydDown(input_nem)
hd_nem.run()
print("✓ NEM complete")

# ============================================================================
# Create Comparison Plots
# ============================================================================
fig, axes = plt.subplots(3, 2, figsize=(14, 12))

# Temperature profiles
ax1 = axes[0, 0]
ax1.plot(hd_eq.time_array, hd_eq.T_fluid - 273.15, 'k-', linewidth=2, label='Equilibrium: T_fluid')
ax1.plot(hd_nem.time_array, hd_nem.T_gas - 273.15, 'r-', linewidth=2, label='NEM: T_gas')
ax1.plot(hd_nem.time_array, hd_nem.T_liquid - 273.15, 'b-', linewidth=2, label='NEM: T_liquid')
ax1.set_xlabel('Time [s]')
ax1.set_ylabel('Temperature [°C]')
ax1.set_title('Phase Temperatures')
ax1.legend()
ax1.grid(True, alpha=0.3)

# Wetted wall temperature
ax2 = axes[0, 1]
ax2.plot(hd_eq.time_array, hd_eq.T_vessel_wetted - 273.15, 'k-', linewidth=2, label='Equilibrium')
ax2.plot(hd_nem.time_array, hd_nem.T_vessel_wetted - 273.15, 'r-', linewidth=2, label='NEM')
ax2.set_xlabel('Time [s]')
ax2.set_ylabel('T_wall wetted [°C]')
ax2.set_title('Wetted Wall Temperature')
ax2.legend()
ax2.grid(True, alpha=0.3)

# Wetted heat transfer coefficient
ax3 = axes[1, 0]
ax3.plot(hd_eq.time_array, hd_eq.h_inside_wetted, 'k-', linewidth=2, label='Equilibrium')
ax3.plot(hd_nem.time_array, hd_nem.h_inside_wetted, 'r-', linewidth=2, label='NEM')
ax3.set_xlabel('Time [s]')
ax3.set_ylabel('h_inside_wetted [W/m²K]')
ax3.set_title('Wetted Internal Heat Transfer Coefficient')
ax3.legend()
ax3.grid(True, alpha=0.3)

# Wetted heat flux
ax4 = axes[1, 1]
ax4.plot(hd_eq.time_array, hd_eq.q_inner_wetted/1000, 'k-', linewidth=2, label='Equilibrium')
ax4.plot(hd_nem.time_array, hd_nem.q_inner_wetted/1000, 'r-', linewidth=2, label='NEM')
ax4.set_xlabel('Time [s]')
ax4.set_ylabel('q_inner_wetted [kW/m²]')
ax4.set_title('Wetted Internal Heat Flux')
ax4.legend()
ax4.grid(True, alpha=0.3)

# Wetted heat transfer rate
ax5 = axes[2, 0]
ax5.plot(hd_eq.time_array, hd_eq.Q_inner_wetted/1000, 'k-', linewidth=2, label='Equilibrium')
ax5.plot(hd_nem.time_array, hd_nem.Q_inner_wetted/1000, 'r-', linewidth=2, label='NEM')
ax5.set_xlabel('Time [s]')
ax5.set_ylabel('Q_inner_wetted [kW]')
ax5.set_title('Wetted Internal Heat Transfer Rate')
ax5.legend()
ax5.grid(True, alpha=0.3)

# Temperature driving force (wall - liquid)
ax6 = axes[2, 1]
dT_eq = hd_eq.T_vessel_wetted - hd_eq.T_fluid
dT_nem = hd_nem.T_vessel_wetted - hd_nem.T_liquid
ax6.plot(hd_eq.time_array, dT_eq, 'k-', linewidth=2, label='Equilibrium: ΔT_wall-fluid')
ax6.plot(hd_nem.time_array, dT_nem, 'r-', linewidth=2, label='NEM: ΔT_wall-liquid')
ax6.set_xlabel('Time [s]')
ax6.set_ylabel('ΔT [K]')
ax6.set_title('Temperature Driving Force (Wall - Liquid)')
ax6.legend()
ax6.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('wetted_htc_equilibrium_vs_nem_h500.pdf', dpi=300, bbox_inches='tight')
plt.savefig('wetted_htc_equilibrium_vs_nem_h500.png', dpi=150, bbox_inches='tight')

print("\n✓ Plots saved:")
print("  - wetted_htc_equilibrium_vs_nem_h500.pdf")
print("  - wetted_htc_equilibrium_vs_nem_h500.png")

# ============================================================================
# Statistical Comparison
# ============================================================================
print("\n" + "="*80)
print("STATISTICAL COMPARISON AT KEY TIMES")
print("="*80)

times_to_check = [150, 200, 250, 300, 350, 400]

print(f"\n{'Time':<8} {'Model':<12} {'T_liq [°C]':<12} {'T_wall [°C]':<12} {'ΔT [K]':<10} {'h_wet [W/m²K]':<15} {'q_wet [kW/m²]':<12}")
print("-"*95)

for t in times_to_check:
    # Equilibrium
    idx_eq = np.argmin(np.abs(hd_eq.time_array - t))
    T_liq_eq = hd_eq.T_fluid[idx_eq] - 273.15
    T_wall_eq = hd_eq.T_vessel_wetted[idx_eq] - 273.15
    dT_eq = hd_eq.T_vessel_wetted[idx_eq] - hd_eq.T_fluid[idx_eq]
    h_wet_eq = hd_eq.h_inside_wetted[idx_eq]
    q_wet_eq = hd_eq.q_inner_wetted[idx_eq] / 1000

    # NEM
    idx_nem = np.argmin(np.abs(hd_nem.time_array - t))
    T_liq_nem = hd_nem.T_liquid[idx_nem] - 273.15
    T_wall_nem = hd_nem.T_vessel_wetted[idx_nem] - 273.15
    dT_nem = hd_nem.T_vessel_wetted[idx_nem] - hd_nem.T_liquid[idx_nem]
    h_wet_nem = hd_nem.h_inside_wetted[idx_nem]
    q_wet_nem = hd_nem.q_inner_wetted[idx_nem] / 1000

    print(f"{t:<8.0f} {'Equil':<12} {T_liq_eq:<12.1f} {T_wall_eq:<12.1f} {dT_eq:<10.2f} {h_wet_eq:<15.1f} {q_wet_eq:<12.1f}")
    print(f"{'':<8} {'NEM':<12} {T_liq_nem:<12.1f} {T_wall_nem:<12.1f} {dT_nem:<10.2f} {h_wet_nem:<15.1f} {q_wet_nem:<12.1f}")
    print(f"{'':<8} {'Difference':<12} {T_liq_nem-T_liq_eq:<12.1f} {T_wall_nem-T_wall_eq:<12.1f} {dT_nem-dT_eq:<10.2f} {h_wet_nem-h_wet_eq:<15.1f} {q_wet_nem-q_wet_eq:<12.1f}")
    print()

# Average values over last 100s
idx_eq_avg = hd_eq.time_array > 300
idx_nem_avg = hd_nem.time_array > 300

print("="*80)
print("AVERAGE VALUES (t > 300s)")
print("="*80)

print(f"\n{'Parameter':<35} {'Equilibrium':<15} {'NEM':<15} {'Difference':<15}")
print("-"*80)
print(f"{'T_liquid [°C]':<35} {np.mean(hd_eq.T_fluid[idx_eq_avg])-273.15:<15.1f} {np.mean(hd_nem.T_liquid[idx_nem_avg])-273.15:<15.1f} {np.mean(hd_nem.T_liquid[idx_nem_avg])-np.mean(hd_eq.T_fluid[idx_eq_avg]):<15.2f}")
print(f"{'T_wall_wetted [°C]':<35} {np.mean(hd_eq.T_vessel_wetted[idx_eq_avg])-273.15:<15.1f} {np.mean(hd_nem.T_vessel_wetted[idx_nem_avg])-273.15:<15.1f} {np.mean(hd_nem.T_vessel_wetted[idx_nem_avg])-np.mean(hd_eq.T_vessel_wetted[idx_eq_avg]):<15.2f}")
print(f"{'ΔT (wall-liquid) [K]':<35} {np.mean(hd_eq.T_vessel_wetted[idx_eq_avg]-hd_eq.T_fluid[idx_eq_avg]):<15.2f} {np.mean(hd_nem.T_vessel_wetted[idx_nem_avg]-hd_nem.T_liquid[idx_nem_avg]):<15.2f} {np.mean(hd_nem.T_vessel_wetted[idx_nem_avg]-hd_nem.T_liquid[idx_nem_avg])-np.mean(hd_eq.T_vessel_wetted[idx_eq_avg]-hd_eq.T_fluid[idx_eq_avg]):<15.2f}")
print(f"{'h_inside_wetted [W/m²K]':<35} {np.mean(hd_eq.h_inside_wetted[idx_eq_avg]):<15.1f} {np.mean(hd_nem.h_inside_wetted[idx_nem_avg]):<15.1f} {np.mean(hd_nem.h_inside_wetted[idx_nem_avg])-np.mean(hd_eq.h_inside_wetted[idx_eq_avg]):<15.1f}")
print(f"{'q_inner_wetted [kW/m²]':<35} {np.mean(hd_eq.q_inner_wetted[idx_eq_avg])/1000:<15.1f} {np.mean(hd_nem.q_inner_wetted[idx_nem_avg])/1000:<15.1f} {(np.mean(hd_nem.q_inner_wetted[idx_nem_avg])-np.mean(hd_eq.q_inner_wetted[idx_eq_avg]))/1000:<15.1f}")

print("\n" + "="*80)
