#!/usr/bin/env python3
"""
Detailed comparison of all array attributes between Equilibrium and NEM models.
"""

import yaml
import sys
import numpy as np
sys.path.insert(0, 'src')

from hyddown.hdclass import HydDown

print("="*80)
print("DETAILED ARRAY COMPARISON: Equilibrium vs NEM")
print("="*80)

# Load input
with open('src/hyddown/examples/nem_propane_psv_fire.yml') as f:
    input_data = yaml.load(f, Loader=yaml.FullLoader)

# Run equilibrium
print("\n[1/2] Running Equilibrium Model...")
input_eq = input_data.copy()
input_eq['calculation'] = input_data['calculation'].copy()
input_eq['calculation']['non_equilibrium'] = False

hdown_eq = HydDown(input_eq)
hdown_eq.run(disable_pbar=True)
print("✓ Complete")

# Run NEM
print("[2/2] Running NEM Model...")
input_nem = input_data.copy()
input_nem['calculation'] = input_data['calculation'].copy()
input_nem['calculation']['non_equilibrium'] = True

hdown_nem = HydDown(input_nem)
hdown_nem.run(disable_pbar=True)
print("✓ Complete")

print("\n" + "="*80)
print("ARRAY-BY-ARRAY COMPARISON")
print("="*80)

# Compare common arrays
common_arrays = [
    ('time_array', 'Time'),
    ('P', 'Pressure [bar]', 1e5),
    ('T_fluid', 'Fluid Temperature [K]'),
    ('mass_fluid', 'Total Mass [kg]'),
    ('mass_rate', 'Mass Flow Rate [kg/s]'),
    ('liquid_level', 'Liquid Level [m]'),
    ('Q_inner', 'Heat Transfer (dry wall) [kW]', 1e3),
    ('Q_inner_wetted', 'Heat Transfer (wetted wall) [kW]', 1e3),
]

for item in common_arrays:
    if len(item) == 2:
        attr, label = item
        scale = 1.0
    else:
        attr, label, scale = item

    if hasattr(hdown_eq, attr) and hasattr(hdown_nem, attr):
        eq_val = getattr(hdown_eq, attr) / scale
        nem_val = getattr(hdown_nem, attr) / scale

        print(f"\n{label}:")
        print(f"  Equilibrium:")
        print(f"    Initial: {eq_val[0]:.3f}, Final: {eq_val[-1]:.3f}")
        print(f"    Min: {np.min(eq_val):.3f}, Max: {np.max(eq_val):.3f}, Mean: {np.mean(eq_val):.3f}")

        print(f"  NEM:")
        print(f"    Initial: {nem_val[0]:.3f}, Final: {nem_val[-1]:.3f}")
        print(f"    Min: {np.min(nem_val):.3f}, Max: {np.max(nem_val):.3f}, Mean: {np.mean(nem_val):.3f}")

        print(f"  Difference (NEM - Eq):")
        diff = nem_val - eq_val
        print(f"    Final: {diff[-1]:.3f}")
        print(f"    Max abs diff: {np.max(np.abs(diff)):.3f}")
        print(f"    RMS diff: {np.sqrt(np.mean(diff**2)):.3f}")

# NEM-specific arrays
print("\n" + "="*80)
print("NEM-SPECIFIC ARRAYS")
print("="*80)

nem_only_arrays = [
    ('T_gas', 'Gas Temperature [K]'),
    ('T_liquid', 'Liquid Temperature [K]'),
    ('m_gas', 'Gas Mass [kg]'),
    ('m_liquid', 'Liquid Mass [kg]'),
    ('U_gas', 'Gas Specific Internal Energy [kJ/kg]', 1e3),
    ('U_liquid', 'Liquid Specific Internal Energy [kJ/kg]', 1e3),
    ('rho_gas', 'Gas Density [kg/m³]'),
    ('rho_liquid', 'Liquid Density [kg/m³]'),
    ('mdot_phase_transfer', 'Phase Transfer Rate [kg/s]'),
]

for item in nem_only_arrays:
    if len(item) == 2:
        attr, label = item
        scale = 1.0
    else:
        attr, label, scale = item

    if hasattr(hdown_nem, attr):
        val = getattr(hdown_nem, attr) / scale

        print(f"\n{label}:")
        print(f"  Initial: {val[0]:.3f}, Final: {val[-1]:.3f}")
        print(f"  Min: {np.min(val):.3f}, Max: {np.max(val):.3f}, Mean: {np.mean(val):.3f}")

        # Special analysis for phase transfer
        if attr == 'mdot_phase_transfer':
            evap = val[val < 0]
            cond = val[val > 0]
            print(f"  Evaporation timesteps: {len(evap)} ({len(evap)/len(val)*100:.1f}%)")
            print(f"  Condensation timesteps: {len(cond)} ({len(cond)/len(val)*100:.1f}%)")
            if len(evap) > 0:
                print(f"  Mean evap rate: {np.mean(evap):.3e} kg/s")
            if len(cond) > 0:
                print(f"  Mean cond rate: {np.mean(cond):.3e} kg/s")

# Energy analysis
print("\n" + "="*80)
print("ENERGY BALANCE ANALYSIS")
print("="*80)

dt = hdown_eq.tstep

# Equilibrium
Q_total_eq = np.sum((hdown_eq.Q_inner + hdown_eq.Q_inner_wetted) * dt)
mass_discharged_eq = hdown_eq.mass_fluid[0] - hdown_eq.mass_fluid[-1]

# NEM
Q_total_nem = np.sum((hdown_nem.Q_inner + hdown_nem.Q_inner_wetted) * dt)
mass_discharged_nem = hdown_nem.mass_fluid[0] - hdown_nem.mass_fluid[-1]

print(f"\nTotal Heat Input:")
print(f"  Equilibrium:      {Q_total_eq/1e6:.2f} MJ")
print(f"  Non-Equilibrium:  {Q_total_nem/1e6:.2f} MJ")
print(f"  Difference:       {(Q_total_nem - Q_total_eq)/1e6:.2f} MJ ({abs(Q_total_eq - Q_total_nem)/Q_total_eq*100:.1f}%)")

print(f"\nMass Discharge:")
print(f"  Equilibrium:      {mass_discharged_eq:.2f} kg ({mass_discharged_eq/hdown_eq.mass_fluid[0]*100:.1f}%)")
print(f"  Non-Equilibrium:  {mass_discharged_nem:.2f} kg ({mass_discharged_nem/hdown_nem.mass_fluid[0]*100:.1f}%)")

print(f"\nHeat Distribution (NEM):")
Q_gas = np.sum(hdown_nem.Q_inner * dt)
Q_liquid = np.sum(hdown_nem.Q_inner_wetted * dt)
print(f"  To gas (dry wall):     {Q_gas/1e6:.2f} MJ ({Q_gas/Q_total_nem*100:.1f}%)")
print(f"  To liquid (wet wall):  {Q_liquid/1e6:.2f} MJ ({Q_liquid/Q_total_nem*100:.1f}%)")

# Phase transfer analysis
total_phase_transfer = np.sum(hdown_nem.mdot_phase_transfer * dt)
print(f"\nPhase Transfer (NEM):")
print(f"  Net transfer: {total_phase_transfer:.2f} kg")
print(f"  Initial liquid: {hdown_nem.m_liquid[0]:.2f} kg")
print(f"  Final liquid: {hdown_nem.m_liquid[-1]:.2f} kg")
print(f"  Liquid evaporated: {hdown_nem.m_liquid[0] - hdown_nem.m_liquid[-1]:.2f} kg")

# PSV operation check
print("\n" + "="*80)
print("PSV OPERATION")
print("="*80)

Pset = input_data['valve']['set_pressure'] / 1e5  # bar
print(f"PSV Set Pressure: {Pset:.2f} bar")

eq_above_set = hdown_eq.P > input_data['valve']['set_pressure']
nem_above_set = hdown_nem.P > input_data['valve']['set_pressure']

print(f"\nEquilibrium:")
print(f"  Timesteps above set: {np.sum(eq_above_set)} ({np.sum(eq_above_set)/len(hdown_eq.P)*100:.1f}%)")
print(f"  Max pressure: {np.max(hdown_eq.P)/1e5:.2f} bar")

print(f"\nNon-Equilibrium:")
print(f"  Timesteps above set: {np.sum(nem_above_set)} ({np.sum(nem_above_set)/len(hdown_nem.P)*100:.1f}%)")
print(f"  Max pressure: {np.max(hdown_nem.P)/1e5:.2f} bar")

print("\n" + "="*80)
