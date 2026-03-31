"""
Compare single-sided (calc) and two-sided (calc_two_sided) h_gl models
against Moodie (1998) 22% fill experiment with time-shifted data.
"""

import os
import sys
import copy

_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
_REPO_ROOT = os.path.abspath(os.path.join(_SCRIPT_DIR, '..', '..', '..'))
sys.path.insert(0, os.path.join(_REPO_ROOT, 'src'))

import yaml
import numpy as np
import matplotlib.pyplot as plt
from hyddown import HydDown

TIME_SHIFT = -120  # seconds

print("=" * 90)
print("MOODIE (1998) 22% FILL - h_gl MODEL COMPARISON")
print("=" * 90)
print()

# Load base input
with open(os.path.join(_SCRIPT_DIR, 'nem_moodie_22perc.yml')) as f:
    base_input = yaml.load(f, Loader=yaml.FullLoader)

# --- Run single-sided model ---
print("Running single-sided model (h_gas_liquid: calc)...")
input_single = copy.deepcopy(base_input)
input_single["calculation"]["h_gas_liquid"] = "calc"
hd_single = HydDown(input_single)
hd_single.run()
print(f"  Completed: {len(hd_single.time_array)} steps")

# --- Run two-sided model ---
print("Running two-sided model (h_gas_liquid: calc_two_sided)...")
input_two = copy.deepcopy(base_input)
input_two["calculation"]["h_gas_liquid"] = "calc_two_sided"
hd_two = HydDown(input_two)
hd_two.run()
print(f"  Completed: {len(hd_two.time_array)} steps")
print()

# --- Load experimental data ---
print("Loading experimental data...")
exp_pressure = np.loadtxt(os.path.join(_SCRIPT_DIR, 'Pres_22_perc_fill.txt'))
exp_gas_temp = np.loadtxt(os.path.join(_SCRIPT_DIR, 'Gas_temp_22_perc_fill.txt'))
exp_liq_temp = np.loadtxt(os.path.join(_SCRIPT_DIR, 'Liq_temp_22_perc_fill.txt'))
exp_peak_temp = np.loadtxt(os.path.join(_SCRIPT_DIR, 'Peak_temp_22_perc_fill.txt'))
exp_mass_loss = np.loadtxt(os.path.join(_SCRIPT_DIR, 'Mass_loss_22_perc_fill.txt'))

# Time-shift experimental data
exp_pres_time = exp_pressure[:, 0] + TIME_SHIFT
exp_pres_val = exp_pressure[:, 1]  # bara
exp_gas_time = exp_gas_temp[:, 0] + TIME_SHIFT
exp_gas_val = exp_gas_temp[:, 1]  # C
exp_liq_time = exp_liq_temp[:, 0] + TIME_SHIFT
exp_liq_val = exp_liq_temp[:, 1]  # C
exp_peak_time = exp_peak_temp[:, 0] + TIME_SHIFT
exp_peak_val = exp_peak_temp[:, 1]  # C
exp_mass_time = exp_mass_loss[:, 0] + TIME_SHIFT
exp_mass_val = exp_mass_loss[:, 1]  # kg
print()

# --- Helper to extract results ---
def extract_results(hd):
    return {
        "time": hd.time_array,
        "pres": hd.P / 1e5,
        "T_gas": hd.T_gas - 273.15,
        "T_liq": hd.T_liquid - 273.15,
        "h_gl": hd.h_gas_liquid,
        "mass_loss": hd.mass_fluid[0] - hd.mass_fluid,
    }

res_s = extract_results(hd_single)
res_t = extract_results(hd_two)

# --- Print h_gl statistics ---
print("=" * 90)
print("h_gl STATISTICS")
print("=" * 90)
for label, res in [("Single-sided", res_s), ("Two-sided", res_t)]:
    h = res["h_gl"]
    mask = h > 0
    if np.any(mask):
        print(f"  {label:16s}: min={np.min(h[mask]):.1f}  mean={np.mean(h[mask]):.1f}  max={np.max(h[mask]):.1f} W/m2K")
    else:
        print(f"  {label:16s}: no active h_gl values")
print()

# --- Compute errors (time-shifted, positive times only) ---
def compute_errors(res, exp_time, exp_val, var_key):
    valid = exp_time >= 0
    if not np.any(valid):
        return None
    sim_interp = np.interp(exp_time[valid], res["time"], res[var_key])
    errors = np.abs(sim_interp - exp_val[valid])
    return {"mae": np.mean(errors), "max": np.max(errors)}

print("=" * 90)
print("VALIDATION METRICS (time-shifted)")
print("=" * 90)

for label, res in [("Single-sided", res_s), ("Two-sided", res_t)]:
    print(f"\n  --- {label} ---")
    e = compute_errors(res, exp_pres_time, exp_pres_val, "pres")
    if e:
        print(f"  Pressure       MAE: {e['mae']:.2f} bar   Max: {e['max']:.2f} bar")
    e = compute_errors(res, exp_gas_time, exp_gas_val, "T_gas")
    if e:
        print(f"  Gas temp       MAE: {e['mae']:.1f} C     Max: {e['max']:.1f} C")
    e = compute_errors(res, exp_liq_time, exp_liq_val, "T_liq")
    if e:
        print(f"  Liquid temp    MAE: {e['mae']:.1f} C     Max: {e['max']:.1f} C")
    e = compute_errors(res, exp_mass_time, exp_mass_val, "mass_loss")
    if e:
        print(f"  Mass loss      MAE: {e['mae']:.1f} kg    Max: {e['max']:.1f} kg")
print()

# --- Plots ---
print("Creating comparison plots...")

fig, axes = plt.subplots(4, 1, figsize=(14, 16))

# Pressure
ax = axes[0]
ax.plot(res_s["time"], res_s["pres"], 'b-', lw=2, label='Single-sided (calc)')
ax.plot(res_t["time"], res_t["pres"], 'g--', lw=2, label='Two-sided (calc_two_sided)')
ax.plot(exp_pres_time, exp_pres_val, 'ro', ms=5, alpha=0.7, label='Moodie (1998)')
ax.set_ylabel('Pressure [bara]', fontsize=14)
ax.set_title('Pressure', fontsize=15, fontweight='bold')
ax.set_xlim(0, 800)
ax.grid(True, alpha=0.3)
ax.legend(fontsize=12)
ax.tick_params(labelsize=12)

# Gas temperature
ax = axes[1]
ax.plot(res_s["time"], res_s["T_gas"], 'b-', lw=2, label='Single-sided')
ax.plot(res_t["time"], res_t["T_gas"], 'g--', lw=2, label='Two-sided')
ax.plot(exp_gas_time, exp_gas_val, 'rs', ms=5, alpha=0.7, mfc='none', mew=1.5, label='Moodie Gas')
ax.set_ylabel('Gas Temperature [C]', fontsize=14)
ax.set_title('Gas Temperature', fontsize=15, fontweight='bold')
ax.set_xlim(0, 800)
ax.grid(True, alpha=0.3)
ax.legend(fontsize=12)
ax.tick_params(labelsize=12)

# Liquid temperature
ax = axes[2]
ax.plot(res_s["time"], res_s["T_liq"], 'b-', lw=2, label='Single-sided')
ax.plot(res_t["time"], res_t["T_liq"], 'g--', lw=2, label='Two-sided')
ax.plot(exp_liq_time, exp_liq_val, 'bs', ms=5, alpha=0.7, mfc='none', mew=1.5, label='Moodie Liquid')
ax.set_ylabel('Liquid Temperature [C]', fontsize=14)
ax.set_title('Liquid Temperature', fontsize=15, fontweight='bold')
ax.set_xlim(0, 800)
ax.grid(True, alpha=0.3)
ax.legend(fontsize=12)
ax.tick_params(labelsize=12)

# h_gl comparison
ax = axes[3]
ax.plot(res_s["time"], res_s["h_gl"], 'b-', lw=1.5, label='Single-sided')
ax.plot(res_t["time"], res_t["h_gl"], 'g--', lw=1.5, label='Two-sided')
ax.set_xlabel('Time [s]', fontsize=14)
ax.set_ylabel('h_gl [W/m2K]', fontsize=14)
ax.set_title('Gas-Liquid Heat Transfer Coefficient', fontsize=15, fontweight='bold')
ax.set_xlim(0, 800)
ax.grid(True, alpha=0.3)
ax.legend(fontsize=12)
ax.tick_params(labelsize=12)

plt.tight_layout()
plt.savefig(os.path.join(_SCRIPT_DIR, 'moodie_hgl_comparison.pdf'), bbox_inches='tight')
plt.savefig(os.path.join(_SCRIPT_DIR, 'moodie_hgl_comparison.png'), dpi=300, bbox_inches='tight')
print("Saved: moodie_hgl_comparison.pdf / .png")

print()
print("=" * 90)
print("DONE")
print("=" * 90)
