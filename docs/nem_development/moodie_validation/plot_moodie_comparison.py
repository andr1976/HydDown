"""
Create comparison plots for Moodie (1998) 22% fill validation
"""

import os
import sys

_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
_REPO_ROOT = os.path.abspath(os.path.join(_SCRIPT_DIR, '..', '..', '..'))
sys.path.insert(0, os.path.join(_REPO_ROOT, 'src'))

import yaml
import numpy as np
import matplotlib.pyplot as plt
from hyddown import HydDown

print("=" * 90)
print("MOODIE (1998) 22% FILL - CREATING COMPARISON PLOTS")
print("=" * 90)
print()

# Run NEM simulation
print("Running NEM simulation...")
with open(os.path.join(_REPO_ROOT, 'src/hyddown/examples/nem_propane_psv_fire.yml')) as f:
    input_dict = yaml.load(f, Loader=yaml.FullLoader)

hd = HydDown(input_dict)
hd.run()

print(f"  Simulation completed: {len(hd.time_array)} time steps")
print()

# Load experimental data
print("Loading experimental data...")
exp_pressure = np.loadtxt(os.path.join(_SCRIPT_DIR, 'Pres_22_perc_fill.txt'))
exp_gas_temp = np.loadtxt(os.path.join(_SCRIPT_DIR, 'Gas_temp_22_perc_fill.txt'))
exp_liq_temp = np.loadtxt(os.path.join(_SCRIPT_DIR, 'Liq_temp_22_perc_fill.txt'))
exp_peak_temp = np.loadtxt(os.path.join(_SCRIPT_DIR, 'Peak_temp_22_perc_fill.txt'))
exp_mass_loss = np.loadtxt(os.path.join(_SCRIPT_DIR, 'Mass_loss_22_perc_fill.txt'))

# Extract time series (pressure now correctly scaled in file)
exp_pres_time = exp_pressure[:, 0]
exp_pres_val = exp_pressure[:, 1]  # bara (correctly scaled)

exp_gas_time = exp_gas_temp[:, 0]
exp_gas_val = exp_gas_temp[:, 1]  # °C

exp_liq_time = exp_liq_temp[:, 0]
exp_liq_val = exp_liq_temp[:, 1]  # °C

exp_peak_time = exp_peak_temp[:, 0]
exp_peak_val = exp_peak_temp[:, 1]  # °C

exp_mass_time = exp_mass_loss[:, 0]
exp_mass_val = exp_mass_loss[:, 1]  # kg (accumulated mass loss)

# Convert simulation results
sim_time = hd.time_array
sim_pres = hd.P / 1e5  # Pa to bara
sim_gas_temp = hd.T_gas - 273.15  # K to °C
sim_liq_temp = hd.T_liquid - 273.15  # K to °C

# Max wall temperature (outer wall, unwetted region is typically hottest)
sim_wall_temp = np.maximum(hd.T_outer_wall, hd.T_outer_wall_wetted) - 273.15  # K to °C

# Calculate accumulated mass loss (initial mass - current mass)
initial_mass = hd.mass_fluid[0]
sim_mass_loss = initial_mass - hd.mass_fluid  # kg (accumulated mass discharged)

print("Creating comparison plots...")
print()

# Create main comparison plot: Pressure + Temperature + Mass Loss
fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(14, 13))

# Pressure subplot
ax1.plot(sim_time, sim_pres, 'b-', linewidth=2, label='NEM Simulation')
ax1.plot(exp_pres_time, exp_pres_val, 'ro', markersize=5, alpha=0.7, label='Moodie (1998) Exp')
ax1.set_xlabel('Time [s]', fontsize=12)
ax1.set_ylabel('Pressure [bara]', fontsize=12)
ax1.set_title('Moodie (1998) 22% Fill - NEM Validation: Pressure', fontsize=13, fontweight='bold')
ax1.set_xlim(0, 800)
ax1.grid(True, alpha=0.3)
ax1.legend(loc='best', fontsize=11)

# Temperature subplot (all temperatures together)
ax2.plot(sim_time, sim_gas_temp, 'r-', linewidth=2.5, label='NEM Gas Temperature', alpha=0.8)
ax2.plot(sim_time, sim_liq_temp, 'b-', linewidth=2.5, label='NEM Liquid Temperature', alpha=0.8)
ax2.plot(sim_time, sim_wall_temp, 'orange', linewidth=2.5, label='NEM Max Wall Temperature', alpha=0.8)

ax2.plot(exp_gas_time, exp_gas_val, 'rs', markersize=6, alpha=0.6, label='Exp Gas Temp', markerfacecolor='none', markeredgewidth=1.5)
ax2.plot(exp_liq_time, exp_liq_val, 'bs', markersize=6, alpha=0.6, label='Exp Liquid Temp', markerfacecolor='none', markeredgewidth=1.5)
ax2.plot(exp_peak_time, exp_peak_val, 'o', color='orange', markersize=6, alpha=0.6, label='Exp Peak/Wall Temp', markerfacecolor='none', markeredgewidth=1.5)

ax2.set_ylabel('Temperature [°C]', fontsize=12)
ax2.set_title('Moodie (1998) 22% Fill - NEM Validation: Temperatures', fontsize=13, fontweight='bold')
ax2.set_xlim(0, 800)
ax2.grid(True, alpha=0.3)
ax2.legend(loc='best', fontsize=10, ncol=2)

# Mass loss subplot
ax3.plot(sim_time, sim_mass_loss, 'b-', linewidth=2.5, label='NEM Simulation', alpha=0.8)
ax3.plot(exp_mass_time, exp_mass_val, 'ro', markersize=5, alpha=0.7, label='Moodie (1998) Exp')
ax3.set_xlabel('Time [s]', fontsize=12)
ax3.set_ylabel('Accumulated Mass Loss [kg]', fontsize=12)
ax3.set_title('Moodie (1998) 22% Fill - NEM Validation: Mass Discharge', fontsize=13, fontweight='bold')
ax3.set_xlim(0, 800)
ax3.grid(True, alpha=0.3)
ax3.legend(loc='best', fontsize=11)

plt.tight_layout()
plt.savefig('moodie_comparison.pdf', bbox_inches='tight')
plt.savefig('moodie_comparison.png', dpi=300, bbox_inches='tight')
print("Comparison plot saved: moodie_comparison.pdf and .png")
print()

# Create separate detailed temperature comparison
fig, ax = plt.subplots(1, 1, figsize=(14, 8))

ax.plot(sim_time, sim_gas_temp, 'r-', linewidth=2.5, label='NEM Gas Temperature', alpha=0.8)
ax.plot(sim_time, sim_liq_temp, 'b-', linewidth=2.5, label='NEM Liquid Temperature', alpha=0.8)
ax.plot(sim_time, sim_wall_temp, 'orange', linewidth=2.5, label='NEM Max Wall Temperature', alpha=0.8)

ax.plot(exp_gas_time, exp_gas_val, 'rs', markersize=7, alpha=0.7, label='Exp Gas Temp', markerfacecolor='none', markeredgewidth=2)
ax.plot(exp_liq_time, exp_liq_val, 'bs', markersize=7, alpha=0.7, label='Exp Liquid Temp', markerfacecolor='none', markeredgewidth=2)
ax.plot(exp_peak_time, exp_peak_val, 'o', color='orange', markersize=7, alpha=0.7, label='Exp Peak/Wall Temp', markerfacecolor='none', markeredgewidth=2)

ax.set_xlabel('Time [s]', fontsize=13)
ax.set_ylabel('Temperature [°C]', fontsize=13)
ax.set_title('Moodie (1998) 22% Fill - Temperature Comparison (Gas, Liquid, Wall)', fontsize=14, fontweight='bold')
ax.set_xlim(0, 800)
ax.grid(True, alpha=0.3)
ax.legend(loc='best', fontsize=12)

plt.tight_layout()
plt.savefig('moodie_temperature_comparison.pdf', bbox_inches='tight')
plt.savefig('moodie_temperature_comparison.png', dpi=300, bbox_inches='tight')
print("Temperature comparison plot saved: moodie_temperature_comparison.pdf and .png")
print()

# Calculate and print errors
print("=" * 90)
print("VALIDATION METRICS")
print("=" * 90)

# Pressure errors
sim_pres_interp = np.interp(exp_pres_time, sim_time, sim_pres)
pres_errors = np.abs(sim_pres_interp - exp_pres_val)
print(f"Pressure:")
print(f"  Mean absolute error: {np.mean(pres_errors):.2f} bar ({np.mean(pres_errors/exp_pres_val*100):.1f}%)")
print()

# Gas temperature errors
sim_gas_interp = np.interp(exp_gas_time, sim_time, sim_gas_temp)
gas_errors = np.abs(sim_gas_interp - exp_gas_val)
print(f"Gas Temperature:")
print(f"  Mean absolute error: {np.mean(gas_errors):.1f}°C")
print(f"  Max error: {np.max(gas_errors):.1f}°C")
print()

# Liquid temperature errors
sim_liq_interp = np.interp(exp_liq_time, sim_time, sim_liq_temp)
liq_errors = np.abs(sim_liq_interp - exp_liq_val)
print(f"Liquid Temperature:")
print(f"  Mean absolute error: {np.mean(liq_errors):.1f}°C")
print(f"  Max error: {np.max(liq_errors):.1f}°C")
print()

# Wall/Peak temperature errors
sim_wall_interp = np.interp(exp_peak_time, sim_time, sim_wall_temp)
wall_errors = np.abs(sim_wall_interp - exp_peak_val)
print(f"Peak/Wall Temperature:")
print(f"  Mean absolute error: {np.mean(wall_errors):.1f}°C")
print(f"  Max error: {np.max(wall_errors):.1f}°C")
print()

# Mass loss errors (only up to 800s for comparison)
exp_mass_800 = exp_mass_loss[exp_mass_time <= 800]
exp_mass_time_800 = exp_mass_time[exp_mass_time <= 800]
exp_mass_val_800 = exp_mass_val[exp_mass_time <= 800]
sim_mass_interp = np.interp(exp_mass_time_800, sim_time, sim_mass_loss)
mass_errors = np.abs(sim_mass_interp - exp_mass_val_800)
print(f"Accumulated Mass Loss:")
print(f"  Mean absolute error: {np.mean(mass_errors):.1f} kg")
print(f"  Max error: {np.max(mass_errors):.1f} kg")
print()

print("=" * 90)
print("PLOTS CREATED SUCCESSFULLY")
print("=" * 90)
