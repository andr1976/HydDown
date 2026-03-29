#!/usr/bin/env python
"""
Diagnostic script to understand NEM solver failures
"""
import numpy as np
import matplotlib.pyplot as plt
from CoolProp import CoolProp as CP

# Parameters from the failure case
m_gas = 89.359  # kg
m_liquid = 824.454  # kg
U_gas_end = 69317129.6  # J (approximate from earlier run)
U_liquid_end = 292003641.8  # J (approximate from earlier run)
V_vessel = 4.85  # m3

# Propane fluid objects
fluid_gas = CP.AbstractState("HEOS", "propane")
fluid_liquid = CP.AbstractState("HEOS", "propane")

# Previous pressure
P_prev = 1058804.51 * 2  # ~21 bar (approximate from failure)

# Current bounds
P_min = max(P_prev * 0.5, 1e5)
P_max = P_prev * 2.0

print(f"Diagnostic for NEM Solver Failure")
print(f"=" * 60)
print(f"m_gas: {m_gas:.3f} kg")
print(f"m_liquid: {m_liquid:.3f} kg")
print(f"U_gas_end: {U_gas_end:.1f} J")
print(f"U_liquid_end: {U_liquid_end:.1f} J")
print(f"V_vessel: {V_vessel:.3f} m³")
print(f"P_prev: {P_prev/1e5:.2f} bar")
print(f"P_min: {P_min/1e5:.2f} bar")
print(f"P_max: {P_max/1e5:.2f} bar")
print()

# Define volume residual function
def volume_residual(P):
    """
    For given pressure P, update phases with (P, U) and check volume constraint.
    Returns (V_gas + V_liquid - V_vessel) / V_vessel
    """
    try:
        U_gas_specific = U_gas_end / m_gas
        U_liquid_specific = U_liquid_end / m_liquid

        # Update gas with (P, U)
        fluid_gas.update(CP.PUmass_INPUTS, P, U_gas_specific)
        rho_gas = fluid_gas.rhomass()
        T_gas = fluid_gas.T()

        # Update liquid with (P, U)
        fluid_liquid.update(CP.PUmass_INPUTS, P, U_liquid_specific)
        rho_liquid = fluid_liquid.rhomass()
        T_liquid = fluid_liquid.T()

        # Calculate volumes
        V_gas = m_gas / rho_gas
        V_liquid = m_liquid / rho_liquid
        V_total = V_gas + V_liquid

        # Return normalized volume residual
        residual = (V_total - V_vessel) / V_vessel

        return residual, V_total, V_gas, V_liquid, rho_gas, rho_liquid, T_gas, T_liquid

    except Exception as e:
        return None, None, None, None, None, None, None, str(e)

# Test at bounds
print("Testing at bounds:")
print("-" * 60)

res_min, V_tot_min, V_g_min, V_l_min, rho_g_min, rho_l_min, T_g_min, T_l_min = volume_residual(P_min)
if res_min is not None:
    print(f"At P_min ({P_min/1e5:.2f} bar):")
    print(f"  Residual: {res_min:.6f}")
    print(f"  V_total: {V_tot_min:.4f} m³ (target: {V_vessel:.4f} m³)")
    print(f"  V_gas: {V_g_min:.4f} m³, V_liquid: {V_l_min:.4f} m³")
    print(f"  T_gas: {T_g_min-273.15:.2f} °C, T_liquid: {T_l_min-273.15:.2f} °C")
    print(f"  rho_gas: {rho_g_min:.2f} kg/m³, rho_liquid: {rho_l_min:.2f} kg/m³")
else:
    print(f"At P_min ({P_min/1e5:.2f} bar): CoolProp FAILED - {T_l_min}")

print()

res_max, V_tot_max, V_g_max, V_l_max, rho_g_max, rho_l_max, T_g_max, T_l_max = volume_residual(P_max)
if res_max is not None:
    print(f"At P_max ({P_max/1e5:.2f} bar):")
    print(f"  Residual: {res_max:.6f}")
    print(f"  V_total: {V_tot_max:.4f} m³ (target: {V_vessel:.4f} m³)")
    print(f"  V_gas: {V_g_max:.4f} m³, V_liquid: {V_l_max:.4f} m³")
    print(f"  T_gas: {T_g_max-273.15:.2f} °C, T_liquid: {T_l_max-273.15:.2f} °C")
    print(f"  rho_gas: {rho_g_max:.2f} kg/m³, rho_liquid: {rho_l_max:.2f} kg/m³")
else:
    print(f"At P_max ({P_max/1e5:.2f} bar): CoolProp FAILED - {T_l_max}")

print()

# Check if root is bracketed
if res_min is not None and res_max is not None:
    if np.sign(res_min) != np.sign(res_max):
        print("✓ Root IS bracketed (residuals have opposite signs)")
    else:
        print("✗ Root NOT bracketed (residuals have SAME sign)")
        print(f"  → Need to adjust bounds or starting guess")

        # Try to find where the root might be
        print("\nSearching for root outside bounds...")
        P_test_low = np.linspace(1e5, P_min, 20)
        P_test_high = np.linspace(P_max, P_max*2, 20)

        for P_test in P_test_low:
            res, _, _, _, _, _, _, _ = volume_residual(P_test)
            if res is not None and np.sign(res) != np.sign(res_min):
                print(f"  Found sign change at P = {P_test/1e5:.2f} bar (below P_min)")
                break

        for P_test in P_test_high:
            res, _, _, _, _, _, _, _ = volume_residual(P_test)
            if res is not None and np.sign(res) != np.sign(res_max):
                print(f"  Found sign change at P = {P_test/1e5:.2f} bar (above P_max)")
                break

# Plot residual function
print("\nGenerating residual plot...")
P_range = np.linspace(max(1e5, P_min*0.5), min(P_max*1.5, 100e5), 200)
residuals = []
V_totals = []
successful = []

for P in P_range:
    res, V_tot, _, _, _, _, _, _ = volume_residual(P)
    if res is not None:
        residuals.append(res)
        V_totals.append(V_tot)
        successful.append(True)
    else:
        residuals.append(np.nan)
        V_totals.append(np.nan)
        successful.append(False)

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))

# Plot residual
ax1.plot(P_range/1e5, residuals, 'b-', linewidth=2, label='Volume Residual')
ax1.axhline(y=0, color='r', linestyle='--', label='Target (residual=0)')
ax1.axvline(x=P_min/1e5, color='g', linestyle=':', label=f'P_min ({P_min/1e5:.1f} bar)')
ax1.axvline(x=P_max/1e5, color='purple', linestyle=':', label=f'P_max ({P_max/1e5:.1f} bar)')
ax1.set_xlabel('Pressure [bar]', fontsize=12)
ax1.set_ylabel('Volume Residual [(V-V_vessel)/V_vessel]', fontsize=12)
ax1.set_title('NEM Solver: Volume Residual vs Pressure', fontsize=14, fontweight='bold')
ax1.legend()
ax1.grid(True, alpha=0.3)

# Plot total volume
ax2.plot(P_range/1e5, V_totals, 'b-', linewidth=2, label='Calculated Volume')
ax2.axhline(y=V_vessel, color='r', linestyle='--', label=f'Target Volume ({V_vessel} m³)')
ax2.axvline(x=P_min/1e5, color='g', linestyle=':', label=f'P_min ({P_min/1e5:.1f} bar)')
ax2.axvline(x=P_max/1e5, color='purple', linestyle=':', label=f'P_max ({P_max/1e5:.1f} bar)')
ax2.set_xlabel('Pressure [bar]', fontsize=12)
ax2.set_ylabel('Total Volume [m³]', fontsize=12)
ax2.set_title('NEM Solver: Total Volume vs Pressure', fontsize=14, fontweight='bold')
ax2.legend()
ax2.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('nem_solver_diagnostics.png', dpi=150, bbox_inches='tight')
print("Plot saved as: nem_solver_diagnostics.png")

# Recommendations
print("\n" + "=" * 60)
print("RECOMMENDATIONS:")
print("=" * 60)
if res_min is not None and res_max is not None:
    if np.sign(res_min) == np.sign(res_max):
        if res_min > 0 and res_max > 0:
            print("• Both bounds give V_total > V_vessel")
            print("  → Need HIGHER pressure (try expanding P_max)")
        else:
            print("• Both bounds give V_total < V_vessel")
            print("  → Need LOWER pressure (try expanding P_min)")

    # Check if bounds are reasonable
    if abs(res_min) > 0.1 or abs(res_max) > 0.1:
        print("• Residuals at bounds are large (>10%)")
        print("  → Consider using previous pressure as initial guess")
        print("  → Or use adaptive bounds based on residual signs")
else:
    print("• CoolProp update failing at bounds")
    print("  → Check thermodynamic state validity")
    print("  → May need to constrain energy/pressure ranges")
