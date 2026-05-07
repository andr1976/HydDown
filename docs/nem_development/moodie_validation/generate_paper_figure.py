"""
Generate publication figure for Moodie (1988) validation using two-sided h_gl model.
Produces 3-panel plot matching paper style: pressure, temperatures, mass discharge.
"""

import os
import sys

_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
_REPO_ROOT = os.path.abspath(os.path.join(_SCRIPT_DIR, '..', '..', '..'))
sys.path.insert(0, os.path.join(_REPO_ROOT, 'src'))

import yaml
import numpy as np
import matplotlib.pyplot as plt
import scienceplots  # noqa: F401
from hyddown import HydDown

plt.style.use(['science', 'ieee'])
plt.rcParams['font.size'] = 9

TIME_SHIFT = -120  # seconds

# Run two-sided model
print("Running two-sided h_gl model...")
with open(os.path.join(_SCRIPT_DIR, 'nem_moodie_22perc.yml')) as f:
    input_dict = yaml.load(f, Loader=yaml.FullLoader)

input_dict["calculation"]["h_gas_liquid"] = "calc_two_sided"
hd = HydDown(input_dict)
hd.run()

# Load experimental data
exp_pressure = np.loadtxt(os.path.join(_SCRIPT_DIR, 'Pres_22_perc_fill.txt'))
exp_gas_temp = np.loadtxt(os.path.join(_SCRIPT_DIR, 'Gas_temp_22_perc_fill.txt'))
exp_liq_temp = np.loadtxt(os.path.join(_SCRIPT_DIR, 'Liq_temp_22_perc_fill.txt'))
exp_mass_loss = np.loadtxt(os.path.join(_SCRIPT_DIR, 'Mass_loss_22_perc_fill.txt'))

# Time-shift
exp_pres_time = exp_pressure[:, 0] + TIME_SHIFT
exp_pres_val = exp_pressure[:, 1]
exp_gas_time = exp_gas_temp[:, 0] + TIME_SHIFT
exp_gas_val = exp_gas_temp[:, 1]
exp_liq_time = exp_liq_temp[:, 0] + TIME_SHIFT
exp_liq_val = exp_liq_temp[:, 1]
exp_mass_time = exp_mass_loss[:, 0] + TIME_SHIFT
exp_mass_val = exp_mass_loss[:, 1]

# Simulation results
sim_time = hd.time_array
sim_pres = hd.P / 1e5
sim_gas_temp = hd.T_gas - 273.15
sim_liq_temp = hd.T_liquid - 273.15
sim_mass_loss = hd.mass_fluid[0] - hd.mass_fluid

# Compute validation metrics (positive times only)
valid_p = exp_pres_time >= 0
valid_g = exp_gas_time >= 0
valid_l = exp_liq_time >= 0
valid_m = exp_mass_time >= 0

p_interp = np.interp(exp_pres_time[valid_p], sim_time, sim_pres)
p_mae = np.mean(np.abs(p_interp - exp_pres_val[valid_p]))
p_rel = np.mean(np.abs(p_interp - exp_pres_val[valid_p]) / exp_pres_val[valid_p] * 100)

g_interp = np.interp(exp_gas_time[valid_g], sim_time, sim_gas_temp)
g_mae = np.mean(np.abs(g_interp - exp_gas_val[valid_g]))
g_max = np.max(np.abs(g_interp - exp_gas_val[valid_g]))

l_interp = np.interp(exp_liq_time[valid_l], sim_time, sim_liq_temp)
l_mae = np.mean(np.abs(l_interp - exp_liq_val[valid_l]))
l_max = np.max(np.abs(l_interp - exp_liq_val[valid_l]))

m_interp = np.interp(exp_mass_time[valid_m], sim_time, sim_mass_loss)
m_mae = np.mean(np.abs(m_interp - exp_mass_val[valid_m]))

print()
print("Validation metrics (two-sided h_gl):")
print(f"  Pressure:    MAE = {p_mae:.2f} bar ({p_rel:.1f}%)")
print(f"  Gas temp:    MAE = {g_mae:.1f} C (max {g_max:.1f} C)")
print(f"  Liquid temp: MAE = {l_mae:.1f} C (max {l_max:.1f} C)")
print(f"  Mass loss:   MAE = {m_mae:.1f} kg")
print()

# Create publication figure (3 panels)
fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(7.5, 8))

# Pressure
ax1.plot(sim_time, sim_pres, 'b-', linewidth=1.5, label='HydDown NEM')
ax1.plot(exp_pres_time, exp_pres_val, 'ro', markersize=4, alpha=0.7, label='Moodie et al. (1988)')
ax1.set_xlabel('Time [s]')
ax1.set_ylabel('Pressure [bara]')
ax1.set_title('Moodie et al. (1988) 22% Fill - Pressure')
ax1.set_xlim(0, 800)
ax1.grid(True, alpha=0.3)
ax1.legend(loc='best')

# Temperatures
ax2.plot(sim_time, sim_gas_temp, 'r-', linewidth=1.5, label='NEM Gas Temperature', alpha=0.8)
ax2.plot(sim_time, sim_liq_temp, 'b-', linewidth=1.5, label='NEM Liquid Temperature', alpha=0.8)
ax2.plot(exp_gas_time, exp_gas_val, 'rs', markersize=4, alpha=0.6,
         label='Exp Gas Temp', markerfacecolor='none', markeredgewidth=1.5)
ax2.plot(exp_liq_time, exp_liq_val, 'bs', markersize=4, alpha=0.6,
         label='Exp Liquid Temp', markerfacecolor='none', markeredgewidth=1.5)
ax2.set_ylabel(r'Temperature [$^\circ$C]')
ax2.set_title('Moodie et al. (1988) 22% Fill - Temperatures')
ax2.set_xlim(0, 800)
ax2.grid(True, alpha=0.3)
ax2.legend(loc='best')

# Mass discharge
ax3.plot(sim_time, sim_mass_loss, 'b-', linewidth=1.5, label='HydDown NEM', alpha=0.8)
ax3.plot(exp_mass_time, exp_mass_val, 'ro', markersize=4, alpha=0.7, label='Moodie et al. (1988)')
ax3.set_xlabel('Time [s]')
ax3.set_ylabel('Accumulated Mass Loss [kg]')
ax3.set_title('Moodie et al. (1988) 22% Fill - Mass Discharge')
ax3.set_xlim(0, 800)
ax3.grid(True, alpha=0.3)
ax3.legend(loc='best')

plt.tight_layout()

# Save to both locations
out_local = os.path.join(_SCRIPT_DIR, 'moodie_propane_fire_validation.pdf')
out_paper = os.path.join(_REPO_ROOT, '..', 'blowdown-paper', 'figures', 'moodie_propane_fire_validation.pdf')

plt.savefig(out_local, bbox_inches='tight')
print(f"Saved: {out_local}")

plt.savefig(out_paper, bbox_inches='tight')
print(f"Saved: {out_paper}")

plt.savefig(out_local.replace('.pdf', '.png'), dpi=300, bbox_inches='tight')
print("Done.")
