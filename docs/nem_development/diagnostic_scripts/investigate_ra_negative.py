#!/usr/bin/env python
"""Investigate negative Ra number in h_gl correlation."""

import sys
sys.path.insert(0, 'src')
import CoolProp.CoolProp as CP
import numpy as np

print("\n" + "="*80)
print("INVESTIGATING NEGATIVE Ra NUMBER")
print("="*80)

# Create fluid objects
fluid_gas = CP.AbstractState("HEOS", "Propane")
fluid_liquid = CP.AbstractState("HEOS", "Propane")

# Test range of conditions
pressures = [5e5, 10e5, 15e5]  # Pa
T_gas_values = [280, 300, 350, 400, 500]  # K
T_liquid_values = [278, 300, 320]  # K

print("\nTesting beta (thermal expansion coefficient) across conditions:")
print(f"{'P [bar]':>8} {'Phase':>8} {'T [K]':>8} {'T [°C]':>8} {'Beta [1/K]':>15} {'Quality':>10} {'State':>12}")
print("-" * 85)

for P in pressures:
    for T in T_gas_values:
        try:
            # Try gas phase
            fluid_gas.update(CP.PT_INPUTS, P, T)
            beta_gas = fluid_gas.isobaric_expansion_coefficient()
            Q = fluid_gas.Q()
            phase = "Gas"

            # Check state
            if Q >= 0 and Q <= 1:
                state = "Two-phase"
            elif Q < 0:
                state = "Liquid"
            else:
                state = "Supercrit/Gas"

            print(f"{P/1e5:8.1f} {phase:>8} {T:8.1f} {T-273.15:8.1f} {beta_gas:15.6e} {Q:10.4f} {state:>12}")

            if beta_gas < 0:
                print(f"  >>> NEGATIVE BETA DETECTED! <<<")
        except Exception as e:
            print(f"{P/1e5:8.1f} {'Gas':>8} {T:8.1f} {T-273.15:8.1f} {'FAILED':>15} {'---':>10} {str(e)[:20]:>12}")

    for T in T_liquid_values:
        try:
            # Try liquid phase
            fluid_liquid.update(CP.PT_INPUTS, P, T)
            beta_liquid = fluid_liquid.isobaric_expansion_coefficient()
            Q = fluid_liquid.Q()
            phase = "Liquid"

            # Check state
            if Q >= 0 and Q <= 1:
                state = "Two-phase"
            elif Q < 0:
                state = "Subcooled"
            else:
                state = "Superheated"

            print(f"{P/1e5:8.1f} {phase:>8} {T:8.1f} {T-273.15:8.1f} {beta_liquid:15.6e} {Q:10.4f} {state:>12}")

            if beta_liquid < 0:
                print(f"  >>> NEGATIVE BETA DETECTED! <<<")
        except Exception as e:
            print(f"{P/1e5:8.1f} {'Liquid':>8} {T:8.1f} {T-273.15:8.1f} {'FAILED':>15} {'---':>10} {str(e)[:20]:>12}")

# Specific test: Near saturation conditions
print("\n" + "="*80)
print("TESTING NEAR SATURATION (where beta might be problematic)")
print("="*80)

P_test = 10e5  # 10 bar

# Get saturation temperature
fluid_test = CP.AbstractState("HEOS", "Propane")
fluid_test.update(CP.PQ_INPUTS, P_test, 0.0)  # Saturated liquid
T_sat = fluid_test.T()

print(f"\nAt P = {P_test/1e5:.1f} bar:")
print(f"Saturation temperature: {T_sat:.2f} K ({T_sat-273.15:.2f}°C)")

# Test liquid phase near saturation
print(f"\n{'ΔT from Tsat [K]':>18} {'T [K]':>8} {'Beta [1/K]':>15} {'Density [kg/m³]':>18} {'State':>12}")
print("-" * 75)

for dT in [-10, -5, -1, -0.1, -0.01, 0, 0.01, 0.1, 1, 5, 10]:
    T_test = T_sat + dT
    try:
        fluid_test.update(CP.PT_INPUTS, P_test, T_test)
        beta = fluid_test.isobaric_expansion_coefficient()
        rho = fluid_test.rhomass()
        Q = fluid_test.Q()

        if Q >= 0 and Q <= 1:
            state = f"Q={Q:.4f}"
        elif Q < 0:
            state = "Subcooled"
        else:
            state = "Superheated"

        marker = " <<<" if beta < 0 else ""
        print(f"{dT:18.2f} {T_test:8.2f} {beta:15.6e} {rho:18.2f} {state:>12}{marker}")
    except Exception as e:
        print(f"{dT:18.2f} {T_test:8.2f} {'FAILED':>15} {'---':>18} {str(e)[:30]}")

print("\n" + "="*80)
print("ANALYSIS")
print("="*80)

print("""
Beta (thermal expansion coefficient) can be negative for:
1. Near saturation conditions for liquids
2. Water-like substances with density maximum (not applicable to propane)
3. Numerical issues in CoolProp property calculation

For liquids very close to saturation:
- Beta can change sign or become very small
- This represents physical behavior near phase transition
- Using abs(beta) is appropriate for buoyancy calculations
""")

print("="*80)
