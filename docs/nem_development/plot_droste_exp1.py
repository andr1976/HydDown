#!/usr/bin/env python
"""
Plot results for Droste & Schoen Experiment 1 validation
"""
import sys
sys.path.insert(0, 'src')

from hyddown.hdclass import HydDown
import matplotlib.pyplot as plt
import numpy as np
import yaml

# Experimental data from Droste & Schoen (1988) - Test 1
# Extracted from Figure 2(a) - Liquid temperature vs time
exp_temp_time = np.array([0, 120, 240, 360, 480, 600, 720])  # seconds
exp_temp_degC = np.array([10, 30, 45, 57, 65, 69, 72])  # °C

# Extracted from Figure 2(b) - Pressure vs time
exp_pres_time = np.array([0, 60, 120, 180, 240, 300, 340, 400, 500, 600, 720])  # seconds
exp_pres_bar = np.array([5.5, 7, 9, 11, 13, 15, 16.4, 18, 21, 23, 24.5])  # bar

# Key experimental results
exp_psv_opening_time = 340  # seconds (5'40")
exp_rupture_time = 720  # seconds (12'00")
exp_rupture_pressure = 24.5  # bar
exp_rupture_temp = 72  # °C

# Load and run simulation
print("Loading simulation...")
with open("nem_droste_exp1.yml") as infile:
    input_dict = yaml.load(infile, Loader=yaml.FullLoader)

hd = HydDown(input_dict)
hd.run(disable_pbar=True)

print("\nSimulation Results:")
print(f"  End time: {hd.time_array[-1]:.1f} s")
print(f"  Final pressure: {hd.P[-1]/1e5:.2f} bar")
print(f"  Final T_liquid: {hd.T_liquid[-1]-273.15:.1f} °C")
print(f"  Final T_gas: {hd.T_gas[-1]-273.15:.1f} °C")

# Find PSV opening time in simulation
psv_opened = np.where(np.abs(hd.mass_rate) > 1e-6)[0]
if len(psv_opened) > 0:
    sim_psv_opening_time = hd.time_array[psv_opened[0]]
    print(f"  PSV opening time: {sim_psv_opening_time:.1f} s (exp: {exp_psv_opening_time} s)")
else:
    print("  PSV did not open during simulation")

# Create comparison plots
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))

# Plot 1: Pressure comparison
ax1.plot(hd.time_array, hd.P/1e5, 'b-', linewidth=2, label='NEM Simulation')
ax1.plot(exp_pres_time, exp_pres_bar, 'ro', markersize=8, label='Droste & Schoen Exp 1')
ax1.axhline(y=15.6, color='g', linestyle='--', label='PSV Set Pressure (15.6 bar)')
ax1.axvline(x=exp_psv_opening_time, color='gray', linestyle=':', alpha=0.5, label=f'Exp PSV Opening ({exp_psv_opening_time}s)')
if len(psv_opened) > 0:
    ax1.axvline(x=sim_psv_opening_time, color='purple', linestyle=':', alpha=0.5, label=f'Sim PSV Opening ({sim_psv_opening_time:.0f}s)')
ax1.set_ylabel('Pressure [bar]', fontsize=12)
ax1.set_xlabel('Time [s]', fontsize=12)
ax1.set_title('Droste & Schoen (1988) - Experiment 1: Pressure Evolution', fontsize=14, fontweight='bold')
ax1.legend(loc='best', fontsize=10)
ax1.grid(True, alpha=0.3)
ax1.set_xlim(0, max(hd.time_array[-1], 720))

# Plot 2: Liquid temperature comparison
ax2.plot(hd.time_array, hd.T_liquid - 273.15, 'b-', linewidth=2, label='NEM T_liquid')
ax2.plot(hd.time_array, hd.T_gas - 273.15, 'c--', linewidth=1.5, label='NEM T_gas', alpha=0.7)
ax2.plot(exp_temp_time, exp_temp_degC, 'ro', markersize=8, label='Droste & Schoen Exp 1 (liquid)')
ax2.axvline(x=exp_psv_opening_time, color='gray', linestyle=':', alpha=0.5, label=f'Exp PSV Opening ({exp_psv_opening_time}s)')
if len(psv_opened) > 0:
    ax2.axvline(x=sim_psv_opening_time, color='purple', linestyle=':', alpha=0.5, label=f'Sim PSV Opening ({sim_psv_opening_time:.0f}s)')
ax2.set_ylabel('Temperature [°C]', fontsize=12)
ax2.set_xlabel('Time [s]', fontsize=12)
ax2.set_title('Droste & Schoen (1988) - Experiment 1: Temperature Evolution', fontsize=14, fontweight='bold')
ax2.legend(loc='best', fontsize=10)
ax2.grid(True, alpha=0.3)
ax2.set_xlim(0, max(hd.time_array[-1], 720))

plt.tight_layout()
plt.savefig('droste_exp1_validation.png', dpi=150, bbox_inches='tight')
print("\nPlot saved as: droste_exp1_validation.png")

# Additional diagnostic plot: Mass flow and fill level
fig2, (ax3, ax4) = plt.subplots(2, 1, figsize=(12, 10))

# Mass flow through PSV
ax3.plot(hd.time_array, hd.mass_rate, 'g-', linewidth=2, label='PSV Mass Flow Rate')
ax3.set_ylabel('Mass Flow Rate [kg/s]', fontsize=12)
ax3.set_xlabel('Time [s]', fontsize=12)
ax3.set_title('PSV Mass Flow Rate', fontsize=14, fontweight='bold')
ax3.legend(loc='best', fontsize=10)
ax3.grid(True, alpha=0.3)

# Fill level
ax4.plot(hd.time_array, hd.liquid_level, 'b-', linewidth=2, label='Liquid Level')
ax4.axhline(y=0.5, color='r', linestyle='--', label='Initial Fill (50%)')
ax4.set_ylabel('Fill Level [fraction]', fontsize=12)
ax4.set_xlabel('Time [s]', fontsize=12)
ax4.set_title('Tank Fill Level Evolution', fontsize=14, fontweight='bold')
ax4.legend(loc='best', fontsize=10)
ax4.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('droste_exp1_diagnostics.png', dpi=150, bbox_inches='tight')
print("Diagnostic plot saved as: droste_exp1_diagnostics.png")

plt.show()
