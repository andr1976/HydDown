"""
Validate NEM against Moodie (1998) with time-shifted experimental data
Uses scandpower_pool fire and rupture analysis for wall temperature
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
print("MOODIE (1998) 22% FILL - NEM VALIDATION (TIME-SHIFTED, RUPTURE ANALYSIS)")
print("=" * 90)
print()

# Run NEM simulation with scandpower_pool fire
print("Running NEM simulation with scandpower_pool fire...")
with open(os.path.join(_SCRIPT_DIR, 'nem_moodie_22perc.yml')) as f:
    input_dict = yaml.load(f, Loader=yaml.FullLoader)

hd = HydDown(input_dict)
hd.run()

print(f"  Simulation completed: {len(hd.time_array)} time steps")
print()

# Run rupture analysis
print("Running rupture analysis with scandpower_pool_peak fire...")
hd.analyze_rupture()
rupture_time = hd.rupture_time
print(f"  Rupture time: {rupture_time if rupture_time else 'No rupture predicted'}")
print()

# Load experimental data
print("Loading experimental data...")
exp_pressure = np.loadtxt(os.path.join(_SCRIPT_DIR, 'Pres_22_perc_fill.txt'))
exp_gas_temp = np.loadtxt(os.path.join(_SCRIPT_DIR, 'Gas_temp_22_perc_fill.txt'))
exp_liq_temp = np.loadtxt(os.path.join(_SCRIPT_DIR, 'Liq_temp_22_perc_fill.txt'))
exp_peak_temp = np.loadtxt(os.path.join(_SCRIPT_DIR, 'Peak_temp_22_perc_fill.txt'))
exp_mass_loss = np.loadtxt(os.path.join(_SCRIPT_DIR, 'Mass_loss_22_perc_fill.txt'))

# Time-shift experimental data by -120 seconds
TIME_SHIFT = -120  # seconds
print(f"Time-shifting experimental data by {TIME_SHIFT} seconds...")
print()

exp_pres_time = exp_pressure[:, 0] + TIME_SHIFT
exp_pres_val = exp_pressure[:, 1]  # bara

exp_gas_time = exp_gas_temp[:, 0] + TIME_SHIFT
exp_gas_val = exp_gas_temp[:, 1]  # °C

exp_liq_time = exp_liq_temp[:, 0] + TIME_SHIFT
exp_liq_val = exp_liq_temp[:, 1]  # °C

exp_peak_time = exp_peak_temp[:, 0] + TIME_SHIFT
exp_peak_val = exp_peak_temp[:, 1]  # °C

exp_mass_time = exp_mass_loss[:, 0] + TIME_SHIFT
exp_mass_val = exp_mass_loss[:, 1]  # kg

# Convert simulation results
sim_time = hd.time_array
sim_pres = hd.P / 1e5  # Pa to bara
sim_gas_temp = hd.T_gas - 273.15  # K to °C
sim_liq_temp = hd.T_liquid - 273.15  # K to °C

# Max wall temperature from pressure simulation
sim_wall_temp_pressure = np.maximum(hd.T_outer_wall, hd.T_outer_wall_wetted) - 273.15  # K to °C

# Rupture analysis wall temperature (unwetted region is hottest)
if hasattr(hd, 'peak_T_unwetted'):
    # Interpolate rupture wall temperature to simulation time grid
    rupture_wall_time = hd.peak_times
    rupture_wall_temp_K = hd.peak_T_unwetted
    sim_wall_temp_rupture = np.interp(sim_time, rupture_wall_time, rupture_wall_temp_K) - 273.15  # K to °C
    print(f"Using rupture analysis wall temperature (unwetted region)")
else:
    print("Warning: No rupture wall temperature found, using pressure simulation wall temp")
    sim_wall_temp_rupture = sim_wall_temp_pressure

# Calculate accumulated mass loss
initial_mass = hd.mass_fluid[0]
sim_mass_loss = initial_mass - hd.mass_fluid  # kg

print("Creating comparison plots...")
print()

# Create main comparison plot
fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(14, 13))

# Pressure subplot
ax1.plot(sim_time, sim_pres, 'b-', linewidth=2, label='HydDown NEM')
ax1.plot(exp_pres_time, exp_pres_val, 'ro', markersize=6, alpha=0.7, label='Moodie (1998)')
ax1.set_xlabel('Time [s]', fontsize=16)
ax1.set_ylabel('Pressure [bara]', fontsize=16)
ax1.set_title('Moodie (1998) 22% Fill - Pressure', fontsize=17, fontweight='bold')
ax1.set_xlim(0, 800)
ax1.grid(True, alpha=0.3)
ax1.legend(loc='best', fontsize=14)
ax1.tick_params(axis='both', labelsize=14)

# Temperature subplot
ax2.plot(sim_time, sim_gas_temp, 'r-', linewidth=2.5, label='NEM Gas Temperature', alpha=0.8)
ax2.plot(sim_time, sim_liq_temp, 'b-', linewidth=2.5, label='NEM Liquid Temperature', alpha=0.8)
ax2.plot(sim_time, sim_wall_temp_rupture, 'orange', linewidth=2.5, label='NEM Rupture Wall Temp (scandpower_pool_peak)', alpha=0.8)

ax2.plot(exp_gas_time, exp_gas_val, 'rs', markersize=6, alpha=0.6, label='Exp Gas Temp', markerfacecolor='none', markeredgewidth=1.5)
ax2.plot(exp_liq_time, exp_liq_val, 'bs', markersize=6, alpha=0.6, label='Exp Liquid Temp', markerfacecolor='none', markeredgewidth=1.5)
ax2.plot(exp_peak_time, exp_peak_val, 'o', color='orange', markersize=6, alpha=0.6, label='Exp Peak/Wall Temp', markerfacecolor='none', markeredgewidth=1.5)

ax2.set_ylabel('Temperature [°C]', fontsize=16)
ax2.set_title('Moodie (1998) 22% Fill - Temperatures', fontsize=17, fontweight='bold')
ax2.set_xlim(0, 800)
ax2.grid(True, alpha=0.3)
ax2.legend(loc='best', fontsize=13, ncol=2)
ax2.tick_params(axis='both', labelsize=14)

# Mass loss subplot
ax3.plot(sim_time, sim_mass_loss, 'b-', linewidth=2.5, label='HydDown NEM', alpha=0.8)
ax3.plot(exp_mass_time, exp_mass_val, 'ro', markersize=6, alpha=0.7, label='Moodie (1998)')
ax3.set_xlabel('Time [s]', fontsize=16)
ax3.set_ylabel('Accumulated Mass Loss [kg]', fontsize=16)
ax3.set_title('Moodie (1998) 22% Fill - Mass Discharge', fontsize=17, fontweight='bold')
ax3.set_xlim(0, 800)
ax3.grid(True, alpha=0.3)
ax3.legend(loc='best', fontsize=14)
ax3.tick_params(axis='both', labelsize=14)

plt.tight_layout()
plt.savefig('moodie_comparison_timeshift.pdf', bbox_inches='tight')
plt.savefig('moodie_comparison_timeshift.png', dpi=300, bbox_inches='tight')
print("Comparison plot saved: moodie_comparison_timeshift.pdf and .png")
print()

# Create temperature-only comparison
fig, ax = plt.subplots(1, 1, figsize=(14, 8))

ax.plot(sim_time, sim_gas_temp, 'r-', linewidth=2.5, label='NEM Gas Temperature', alpha=0.8)
ax.plot(sim_time, sim_liq_temp, 'b-', linewidth=2.5, label='NEM Liquid Temperature', alpha=0.8)
ax.plot(sim_time, sim_wall_temp_rupture, 'orange', linewidth=2.5, label='NEM Rupture Wall Temp (scandpower_pool_peak)', alpha=0.8)

ax.plot(exp_gas_time, exp_gas_val, 'rs', markersize=7, alpha=0.7, label='Exp Gas Temp', markerfacecolor='none', markeredgewidth=2)
ax.plot(exp_liq_time, exp_liq_val, 'bs', markersize=7, alpha=0.7, label='Exp Liquid Temp', markerfacecolor='none', markeredgewidth=2)
ax.plot(exp_peak_time, exp_peak_val, 'o', color='orange', markersize=7, alpha=0.7, label='Exp Peak/Wall Temp', markerfacecolor='none', markeredgewidth=2)

ax.set_xlabel('Time [s]', fontsize=16)
ax.set_ylabel('Temperature [°C]', fontsize=16)
ax.set_title('Moodie (1998) 22% Fill - Temperature Comparison', fontsize=17, fontweight='bold')
ax.set_xlim(0, 800)
ax.grid(True, alpha=0.3)
ax.legend(loc='best', fontsize=14)
ax.tick_params(axis='both', labelsize=14)

plt.tight_layout()
plt.savefig('moodie_temperature_comparison_timeshift.pdf', bbox_inches='tight')
plt.savefig('moodie_temperature_comparison_timeshift.png', dpi=300, bbox_inches='tight')
print("Temperature comparison plot saved: moodie_temperature_comparison_timeshift.pdf and .png")
print()

# Calculate errors (only for positive time values after shift)
print("=" * 90)
print("VALIDATION METRICS (TIME-SHIFTED DATA)")
print("=" * 90)

# Filter out negative times
valid_pres = exp_pres_time >= 0
valid_gas = exp_gas_time >= 0
valid_liq = exp_liq_time >= 0
valid_peak = exp_peak_time >= 0
valid_mass = exp_mass_time >= 0

# Pressure errors
if np.any(valid_pres):
    sim_pres_interp = np.interp(exp_pres_time[valid_pres], sim_time, sim_pres)
    pres_errors = np.abs(sim_pres_interp - exp_pres_val[valid_pres])
    print(f"Pressure:")
    print(f"  Mean absolute error: {np.mean(pres_errors):.2f} bar ({np.mean(pres_errors/exp_pres_val[valid_pres]*100):.1f}%)")
    print()

# Gas temperature errors
if np.any(valid_gas):
    sim_gas_interp = np.interp(exp_gas_time[valid_gas], sim_time, sim_gas_temp)
    gas_errors = np.abs(sim_gas_interp - exp_gas_val[valid_gas])
    print(f"Gas Temperature:")
    print(f"  Mean absolute error: {np.mean(gas_errors):.1f}°C")
    print(f"  Max error: {np.max(gas_errors):.1f}°C")
    print()

# Liquid temperature errors
if np.any(valid_liq):
    sim_liq_interp = np.interp(exp_liq_time[valid_liq], sim_time, sim_liq_temp)
    liq_errors = np.abs(sim_liq_interp - exp_liq_val[valid_liq])
    print(f"Liquid Temperature:")
    print(f"  Mean absolute error: {np.mean(liq_errors):.1f}°C")
    print(f"  Max error: {np.max(liq_errors):.1f}°C")
    print()

# Rupture wall temperature errors (compared to experimental peak temp)
if np.any(valid_peak):
    sim_wall_rupt_interp = np.interp(exp_peak_time[valid_peak], sim_time, sim_wall_temp_rupture)
    wall_errors = np.abs(sim_wall_rupt_interp - exp_peak_val[valid_peak])
    print(f"Rupture Wall Temperature (vs Exp Peak Temp):")
    print(f"  Mean absolute error: {np.mean(wall_errors):.1f}°C")
    print(f"  Max error: {np.max(wall_errors):.1f}°C")
    print()

# Mass loss errors
if np.any(valid_mass):
    sim_mass_interp = np.interp(exp_mass_time[valid_mass], sim_time, sim_mass_loss)
    mass_errors = np.abs(sim_mass_interp - exp_mass_val[valid_mass])
    print(f"Accumulated Mass Loss:")
    print(f"  Mean absolute error: {np.mean(mass_errors):.1f} kg")
    print(f"  Max error: {np.max(mass_errors):.1f} kg")
    print()

print("=" * 90)
print("PLOTS CREATED SUCCESSFULLY")
print("=" * 90)
