#!/usr/bin/env python
"""Sensitivity study for h_gl (gas-liquid heat transfer coefficient)."""

import yaml
import sys
sys.path.insert(0, 'src')
from hyddown.hdclass import HydDown
import numpy as np
import matplotlib.pyplot as plt

print("\n" + "="*80)
print("H_GL SENSITIVITY STUDY")
print("="*80)
print("\nSettings:")
print("  relax_factor = 0.8")
print("  heat_multiplier = 50.0")
print("  Energy method: (h+u)/2")

# Load base input
with open('src/hyddown/examples/nem_propane_psv_fire.yml') as f:
    input_base = yaml.load(f, Loader=yaml.FullLoader)

# Run equilibrium first
print("\n" + "="*80)
print("Running EQUILIBRIUM simulation...")
print("="*80)
input_eq = input_base.copy()
input_eq['calculation'] = input_base['calculation'].copy()
input_eq['calculation']['non_equilibrium'] = False
hd_eq = HydDown(input_eq)
hd_eq.run()
print("✓ Equilibrium complete")

# Test different h_gl values
h_gl_values = [50, 100, 200, 500, 1000, 2000]
results = {}

for h_gl in h_gl_values:
    print(f"\n{'='*80}")
    print(f"Running NEM with h_gl = {h_gl} W/m²K")
    print(f"{'='*80}")

    input_test = input_base.copy()
    input_test['calculation'] = input_base['calculation'].copy()
    input_test['calculation']['non_equilibrium'] = True
    input_test['calculation']['h_gas_liquid'] = h_gl

    try:
        hd = HydDown(input_test)
        hd.run()

        results[h_gl] = {
            'time': hd.time_array.copy(),
            'P': hd.P.copy(),
            'T_gas': hd.T_gas.copy(),
            'T_liquid': hd.T_liquid.copy(),
            'T_fluid': hd.T_fluid.copy(),
            'm_total': hd.mass_fluid.copy(),
            'm_gas': hd.m_gas.copy(),
            'm_liquid': hd.m_liquid.copy(),
        }

        i_final = len(hd.time_array) - 1
        print(f"✓ Complete")
        print(f"  Final P: {hd.P[i_final]/1e5:.2f} bar")
        print(f"  Final T_gas: {hd.T_gas[i_final]:.2f} K")
        print(f"  Final T_liquid: {hd.T_liquid[i_final]:.2f} K")
        print(f"  Final ΔT: {hd.T_gas[i_final] - hd.T_liquid[i_final]:.2f} K")
        print(f"  Final mass: {hd.mass_fluid[i_final]:.2f} kg")
        print(f"  Max P: {np.max(hd.P)/1e5:.2f} bar")

    except Exception as e:
        print(f"✗ Failed: {e}")

# Create plots
print(f"\n{'='*80}")
print("Creating plots...")
print(f"{'='*80}")

fig = plt.figure(figsize=(16, 10))

# Subplot 1: Pressure
ax1 = plt.subplot(2, 3, 1)
ax1.plot(hd_eq.time_array, hd_eq.P/1e5, 'k-', linewidth=2, label='Equilibrium')
for h_gl in h_gl_values:
    if h_gl in results:
        ax1.plot(results[h_gl]['time'], results[h_gl]['P']/1e5,
                label=f'h_gl={h_gl}')
ax1.axhline(14.3, color='r', linestyle='--', alpha=0.3, label='PSV set')
ax1.set_xlabel('Time [s]')
ax1.set_ylabel('Pressure [bar]')
ax1.set_title('Pressure Evolution')
ax1.legend(fontsize=8)
ax1.grid(True, alpha=0.3)

# Subplot 2: Gas Temperature
ax2 = plt.subplot(2, 3, 2)
ax2.plot(hd_eq.time_array, hd_eq.T_fluid, 'k-', linewidth=2, label='Equilibrium')
for h_gl in h_gl_values:
    if h_gl in results:
        ax2.plot(results[h_gl]['time'], results[h_gl]['T_gas'],
                label=f'h_gl={h_gl}')
ax2.set_xlabel('Time [s]')
ax2.set_ylabel('Gas Temperature [K]')
ax2.set_title('Gas Temperature')
ax2.legend(fontsize=8)
ax2.grid(True, alpha=0.3)

# Subplot 3: Liquid Temperature
ax3 = plt.subplot(2, 3, 3)
ax3.plot(hd_eq.time_array, hd_eq.T_fluid, 'k-', linewidth=2, label='Equilibrium')
for h_gl in h_gl_values:
    if h_gl in results:
        ax3.plot(results[h_gl]['time'], results[h_gl]['T_liquid'],
                label=f'h_gl={h_gl}')
ax3.set_xlabel('Time [s]')
ax3.set_ylabel('Liquid Temperature [K]')
ax3.set_title('Liquid Temperature')
ax3.legend(fontsize=8)
ax3.grid(True, alpha=0.3)

# Subplot 4: Temperature Difference
ax4 = plt.subplot(2, 3, 4)
for h_gl in h_gl_values:
    if h_gl in results:
        dT = results[h_gl]['T_gas'] - results[h_gl]['T_liquid']
        ax4.plot(results[h_gl]['time'], dT, label=f'h_gl={h_gl}')
ax4.axhline(0, color='k', linestyle='--', alpha=0.3)
ax4.set_xlabel('Time [s]')
ax4.set_ylabel('ΔT = T_gas - T_liquid [K]')
ax4.set_title('Temperature Stratification')
ax4.legend(fontsize=8)
ax4.grid(True, alpha=0.3)

# Subplot 5: Total Mass
ax5 = plt.subplot(2, 3, 5)
ax5.plot(hd_eq.time_array, hd_eq.mass_fluid, 'k-', linewidth=2, label='Equilibrium')
for h_gl in h_gl_values:
    if h_gl in results:
        ax5.plot(results[h_gl]['time'], results[h_gl]['m_total'],
                label=f'h_gl={h_gl}')
ax5.set_xlabel('Time [s]')
ax5.set_ylabel('Total Mass [kg]')
ax5.set_title('Total Mass Inventory')
ax5.legend(fontsize=8)
ax5.grid(True, alpha=0.3)

# Subplot 6: Final Results Summary
ax6 = plt.subplot(2, 3, 6)
ax6.axis('off')

# Create summary table
summary_text = "FINAL RESULTS SUMMARY\n"
summary_text += "="*50 + "\n\n"
summary_text += f"{'h_gl':<8} {'P':<8} {'T_gas':<8} {'T_liq':<8} {'ΔT':<8} {'Mass':<8}\n"
summary_text += f"{'[W/m²K]':<8} {'[bar]':<8} {'[K]':<8} {'[K]':<8} {'[K]':<8} {'[kg]':<8}\n"
summary_text += "-"*50 + "\n"

# Equilibrium
i_eq = len(hd_eq.time_array) - 1
summary_text += f"{'Equil':<8} {hd_eq.P[i_eq]/1e5:<8.2f} {hd_eq.T_fluid[i_eq]:<8.1f} {hd_eq.T_fluid[i_eq]:<8.1f} {'0.0':<8} {hd_eq.mass_fluid[i_eq]:<8.0f}\n"

# NEM results
for h_gl in h_gl_values:
    if h_gl in results:
        i_final = len(results[h_gl]['time']) - 1
        P_final = results[h_gl]['P'][i_final] / 1e5
        T_gas = results[h_gl]['T_gas'][i_final]
        T_liq = results[h_gl]['T_liquid'][i_final]
        dT = T_gas - T_liq
        mass = results[h_gl]['m_total'][i_final]
        summary_text += f"{h_gl:<8} {P_final:<8.2f} {T_gas:<8.1f} {T_liq:<8.1f} {dT:<8.1f} {mass:<8.0f}\n"

ax6.text(0.1, 0.5, summary_text, fontsize=9, family='monospace',
         verticalalignment='center', transform=ax6.transAxes)

plt.tight_layout()
plt.savefig('h_gl_sensitivity_study.pdf')
print("✓ Plot saved to 'h_gl_sensitivity_study.pdf'")

print(f"\n{'='*80}")
print("STUDY COMPLETE")
print(f"{'='*80}\n")
