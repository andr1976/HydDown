"""
Validate NEM implementation against Moodie (1998) 22% fill experiment
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
print("MOODIE (1998) 22% FILL VALIDATION - NEM")
print("=" * 90)
print()

# Run NEM simulation
print("Running NEM simulation...")
with open(os.path.join(_REPO_ROOT, 'src/hyddown/examples/nem_propane_psv_fire.yml')) as f:
    input_dict = yaml.load(f, Loader=yaml.FullLoader)

hd = HydDown(input_dict)
hd.run()

print(f"  Simulation completed: {len(hd.time_array)} time steps")
print(f"  End time: {hd.time_array[-1]:.1f} s")
print()

# Load experimental data
print("Loading experimental data...")
exp_pressure = np.loadtxt(os.path.join(_SCRIPT_DIR, 'Pres_22_perc_fill.txt'))
exp_gas_temp = np.loadtxt(os.path.join(_SCRIPT_DIR, 'Gas_temp_22_perc_fill.txt'))
exp_liq_temp = np.loadtxt(os.path.join(_SCRIPT_DIR, 'Liq_temp_22_perc_fill.txt'))

print(f"  Pressure data points: {len(exp_pressure)}")
print(f"  Gas temperature data points: {len(exp_gas_temp)}")
print(f"  Liquid temperature data points: {len(exp_liq_temp)}")
print()

# Extract time series
exp_pres_time = exp_pressure[:, 0]
exp_pres_val = exp_pressure[:, 1]  # bara

exp_gas_time = exp_gas_temp[:, 0]
exp_gas_val = exp_gas_temp[:, 1]  # °C

exp_liq_time = exp_liq_temp[:, 0]
exp_liq_val = exp_liq_temp[:, 1]  # °C

# Convert simulation results to appropriate units
sim_time = hd.time_array
sim_pres = hd.P / 1e5  # Pa to bara
sim_gas_temp = hd.T_gas - 273.15  # K to °C (NEM model)
sim_liq_temp = hd.T_liquid - 273.15  # K to °C (NEM model)

# Find PSV opening time
psv_opening_time = None
for i, mdot in enumerate(hd.mass_rate):
    if abs(mdot) > 1e-6:
        psv_opening_time = hd.time_array[i]
        psv_opening_pres = hd.P[i] / 1e5
        break

if psv_opening_time:
    print(f"PSV Opening:")
    print(f"  Simulated time: {psv_opening_time:.1f} s at {psv_opening_pres:.2f} bara")
    # Find experimental PSV opening (when pressure peaks or starts to plateau)
    exp_max_pres_idx = np.argmax(exp_pres_val)
    exp_psv_opening = exp_pres_time[exp_max_pres_idx]
    exp_psv_pres = exp_pres_val[exp_max_pres_idx]
    print(f"  Experimental time: ~{exp_psv_opening:.1f} s at {exp_psv_pres:.2f} bara")
    print(f"  Error: {psv_opening_time - exp_psv_opening:+.1f} s ({(psv_opening_time - exp_psv_opening)/exp_psv_opening*100:+.1f}%)")
    print()

# Calculate errors at key time points
print("Validation Metrics:")
print("-" * 90)

# Interpolate simulation to experimental time points for pressure
sim_pres_interp = np.interp(exp_pres_time, sim_time, sim_pres)
pres_errors = np.abs(sim_pres_interp - exp_pres_val)
pres_rel_errors = pres_errors / exp_pres_val * 100

print(f"Pressure:")
print(f"  Mean absolute error: {np.mean(pres_errors):.3f} bar")
print(f"  Mean relative error: {np.mean(pres_rel_errors):.1f}%")
print(f"  Max absolute error: {np.max(pres_errors):.3f} bar")
print()

# Interpolate for gas temperature
sim_gas_interp = np.interp(exp_gas_time, sim_time, sim_gas_temp)
gas_errors = np.abs(sim_gas_interp - exp_gas_val)
print(f"Gas Temperature:")
print(f"  Mean absolute error: {np.mean(gas_errors):.1f}°C")
print(f"  Max absolute error: {np.max(gas_errors):.1f}°C")
print()

# Interpolate for liquid temperature
sim_liq_interp = np.interp(exp_liq_time, sim_time, sim_liq_temp)
liq_errors = np.abs(sim_liq_interp - exp_liq_val)
print(f"Liquid Temperature:")
print(f"  Mean absolute error: {np.mean(liq_errors):.1f}°C")
print(f"  Max absolute error: {np.max(liq_errors):.1f}°C")
print()

# Create validation plot
fig, axes = plt.subplots(3, 1, figsize=(12, 10))

# Pressure plot
ax = axes[0]
ax.plot(sim_time, sim_pres, 'b-', linewidth=2, label='NEM Simulation')
ax.plot(exp_pres_time, exp_pres_val, 'ro', markersize=6, label='Moodie (1998) Exp')
if psv_opening_time:
    ax.axvline(psv_opening_time, color='b', linestyle='--', alpha=0.5, label=f'PSV opens ({psv_opening_time:.0f}s)')
    ax.axvline(exp_psv_opening, color='r', linestyle='--', alpha=0.5, label=f'Exp PSV ({exp_psv_opening:.0f}s)')
ax.set_xlabel('Time [s]', fontsize=11)
ax.set_ylabel('Pressure [bara]', fontsize=11)
ax.set_title('Moodie (1998) 22% Fill - NEM Validation: Pressure', fontsize=12, fontweight='bold')
ax.grid(True, alpha=0.3)
ax.legend(loc='best', fontsize=10)

# Gas temperature plot
ax = axes[1]
ax.plot(sim_time, sim_gas_temp, 'b-', linewidth=2, label='NEM Simulation')
ax.plot(exp_gas_time, exp_gas_val, 'ro', markersize=6, label='Moodie (1998) Exp')
ax.set_xlabel('Time [s]', fontsize=11)
ax.set_ylabel('Gas Temperature [°C]', fontsize=11)
ax.set_title('Moodie (1998) 22% Fill - NEM Validation: Gas Temperature', fontsize=12, fontweight='bold')
ax.grid(True, alpha=0.3)
ax.legend(loc='best', fontsize=10)

# Liquid temperature plot
ax = axes[2]
ax.plot(sim_time, sim_liq_temp, 'b-', linewidth=2, label='NEM Simulation')
ax.plot(exp_liq_time, exp_liq_val, 'ro', markersize=6, label='Moodie (1998) Exp')
ax.set_xlabel('Time [s]', fontsize=11)
ax.set_ylabel('Liquid Temperature [°C]', fontsize=11)
ax.set_title('Moodie (1998) 22% Fill - NEM Validation: Liquid Temperature', fontsize=12, fontweight='bold')
ax.grid(True, alpha=0.3)
ax.legend(loc='best', fontsize=10)

plt.tight_layout()
plt.savefig('moodie_nem_validation.png', dpi=300, bbox_inches='tight')
print("Validation plot saved: moodie_nem_validation.png")
print()

# Create error plot
fig, axes = plt.subplots(3, 1, figsize=(12, 10))

# Pressure error
ax = axes[0]
ax.plot(exp_pres_time, pres_errors, 'b-', linewidth=2, marker='o', markersize=4)
ax.axhline(0, color='k', linestyle='--', alpha=0.3)
ax.set_xlabel('Time [s]', fontsize=11)
ax.set_ylabel('Absolute Error [bar]', fontsize=11)
ax.set_title(f'Pressure Error (Mean: {np.mean(pres_errors):.3f} bar)', fontsize=12, fontweight='bold')
ax.grid(True, alpha=0.3)

# Gas temperature error
ax = axes[1]
ax.plot(exp_gas_time, gas_errors, 'b-', linewidth=2, marker='o', markersize=4)
ax.axhline(0, color='k', linestyle='--', alpha=0.3)
ax.set_xlabel('Time [s]', fontsize=11)
ax.set_ylabel('Absolute Error [°C]', fontsize=11)
ax.set_title(f'Gas Temperature Error (Mean: {np.mean(gas_errors):.1f}°C)', fontsize=12, fontweight='bold')
ax.grid(True, alpha=0.3)

# Liquid temperature error
ax = axes[2]
ax.plot(exp_liq_time, liq_errors, 'b-', linewidth=2, marker='o', markersize=4)
ax.axhline(0, color='k', linestyle='--', alpha=0.3)
ax.set_xlabel('Time [s]', fontsize=11)
ax.set_ylabel('Absolute Error [°C]', fontsize=11)
ax.set_title(f'Liquid Temperature Error (Mean: {np.mean(liq_errors):.1f}°C)', fontsize=12, fontweight='bold')
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('moodie_nem_errors.png', dpi=300, bbox_inches='tight')
print("Error plot saved: moodie_nem_errors.png")
print()

print("=" * 90)
print("VALIDATION SUMMARY")
print("=" * 90)
print(f"NEM simulation captures the two-temperature behavior of the Moodie experiment.")
print(f"Pressure: {np.mean(pres_rel_errors):.1f}% mean error")
print(f"Gas Temperature: {np.mean(gas_errors):.1f}°C mean error")
print(f"Liquid Temperature: {np.mean(liq_errors):.1f}°C mean error")
if psv_opening_time:
    psv_error = abs(psv_opening_time - exp_psv_opening)
    print(f"PSV Opening: {psv_error:.1f}s error ({psv_error/exp_psv_opening*100:.1f}%)")
print("=" * 90)
