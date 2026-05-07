#!/usr/bin/env python
"""Trace NEM state updates to verify fluid objects are correct."""

import yaml
import sys
sys.path.insert(0, 'src')
from hyddown.hdclass import HydDown
import numpy as np
import CoolProp.CoolProp as CP

# Load input and enable NEM
with open('src/hyddown/examples/nem_propane_psv_fire.yml') as f:
    input_dict = yaml.load(f, Loader=yaml.FullLoader)
input_dict['calculation']['non_equilibrium'] = True
input_dict['calculation']['end_time'] = 10.0  # Short run

# Run NEM simulation
hd = HydDown(input_dict)
hd.run()

print("\n" + "="*80)
print("NEM STATE VERIFICATION")
print("="*80)

# Check a few timesteps
for i in [5, 10, 15, 20]:
    print(f"\nTimestep i={i}, t={hd.time_array[i]:.2f}s:")
    print(f"  Stored values:")
    print(f"    P[{i}] = {hd.P[i]/1e5:.4f} bar")
    print(f"    T_gas[{i}] = {hd.T_gas[i]:.2f} K")
    print(f"    T_liquid[{i}] = {hd.T_liquid[i]:.2f} K")
    print(f"    rho_gas[{i}] = {hd.rho_gas[i]:.2f} kg/m³")
    print(f"    rho_liquid[{i}] = {hd.rho_liquid[i]:.2f} kg/m³")
    print(f"    U_gas[{i}] = {hd.U_gas[i]:.1f} J/kg (specific)")
    print(f"    U_liquid[{i}] = {hd.U_liquid[i]:.1f} J/kg (specific)")
    print(f"    m_gas[{i}] = {hd.m_gas[i]:.2f} kg")
    print(f"    m_liquid[{i}] = {hd.m_liquid[i]:.2f} kg")

    # Verify by recalculating with CoolProp
    try:
        fluid_gas_check = CP.AbstractState("HEOS", input_dict['initial']['fluid'])
        fluid_gas_check.update(CP.PUmass_INPUTS, hd.P[i], hd.U_gas[i])

        fluid_liquid_check = CP.AbstractState("HEOS", input_dict['initial']['fluid'])
        fluid_liquid_check.update(CP.PUmass_INPUTS, hd.P[i], hd.U_liquid[i])

        print(f"  Verification (recalculate from P, U):")
        print(f"    Gas: T={fluid_gas_check.T():.2f} K, rho={fluid_gas_check.rhomass():.2f} kg/m³, Q={fluid_gas_check.Q():.4f}")
        print(f"    Liquid: T={fluid_liquid_check.T():.2f} K, rho={fluid_liquid_check.rhomass():.2f} kg/m³, Q={fluid_liquid_check.Q():.4f}")

        # Check consistency
        T_gas_error = abs(fluid_gas_check.T() - hd.T_gas[i])
        T_liquid_error = abs(fluid_liquid_check.T() - hd.T_liquid[i])
        rho_gas_error = abs(fluid_gas_check.rhomass() - hd.rho_gas[i])
        rho_liquid_error = abs(fluid_liquid_check.rhomass() - hd.rho_liquid[i])

        if T_gas_error > 0.1 or T_liquid_error > 0.1:
            print(f"  ⚠️  Temperature mismatch!")
        if rho_gas_error > 0.01 or rho_liquid_error > 0.01:
            print(f"  ⚠️  Density mismatch!")

    except Exception as e:
        print(f"  ❌ Verification failed: {e}")

    # Check volume conservation
    V_gas = hd.m_gas[i] / hd.rho_gas[i] if hd.rho_gas[i] > 0 else 0.0
    V_liquid = hd.m_liquid[i] / hd.rho_liquid[i] if hd.rho_liquid[i] > 0 else 0.0
    V_total = V_gas + V_liquid
    V_error = abs(V_total - hd.vol) / hd.vol * 100
    print(f"  Volume check:")
    print(f"    V_gas = {V_gas:.4f} m³")
    print(f"    V_liquid = {V_liquid:.4f} m³")
    print(f"    V_total = {V_total:.4f} m³ (vessel = {hd.vol:.4f} m³)")
    print(f"    Error: {V_error:.2f}%")
    if V_error > 0.1:
        print(f"  ⚠️  Volume conservation violated!")

print("\n" + "="*80)

# Check if gas and liquid fluid objects are in correct state after simulation
print("\nFinal fluid object states:")
print(f"  hd.fluid_gas:")
print(f"    P = {hd.fluid_gas.p()/1e5:.4f} bar")
print(f"    T = {hd.fluid_gas.T():.2f} K")
print(f"    rho = {hd.fluid_gas.rhomass():.2f} kg/m³")
print(f"    U = {hd.fluid_gas.umass():.1f} J/kg")
print(f"    Q = {hd.fluid_gas.Q():.4f}")

print(f"  hd.fluid_liquid:")
print(f"    P = {hd.fluid_liquid.p()/1e5:.4f} bar")
print(f"    T = {hd.fluid_liquid.T():.2f} K")
print(f"    rho = {hd.fluid_liquid.rhomass():.2f} kg/m³")
print(f"    U = {hd.fluid_liquid.umass():.1f} J/kg")
print(f"    Q = {hd.fluid_liquid.Q():.4f}")

i_final = len(hd.time_array) - 1
print(f"\n  Expected from stored arrays (i={i_final}):")
print(f"    P[{i_final}] = {hd.P[i_final]/1e5:.4f} bar")
print(f"    T_gas[{i_final}] = {hd.T_gas[i_final]:.2f} K, rho_gas = {hd.rho_gas[i_final]:.2f} kg/m³")
print(f"    T_liquid[{i_final}] = {hd.T_liquid[i_final]:.2f} K, rho_liquid = {hd.rho_liquid[i_final]:.2f} kg/m³")

print("\n" + "="*80)
