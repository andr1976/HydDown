#!/usr/bin/env python
"""Diagnostic script to investigate Rohsenow nucleate boiling properties in NEM."""

import sys
sys.path.insert(0, 'src')
from hyddown.hdclass import HydDown
from hyddown import transport as tp
import yaml
import numpy as np
from CoolProp import CoolProp as CP

print("\n" + "="*80)
print("ROHSENOW NUCLEATE BOILING PROPERTIES DIAGNOSTIC FOR NEM")
print("="*80)

# Load NEM input
with open('src/hyddown/examples/nem_propane_psv_fire.yml') as f:
    input_data = yaml.load(f, Loader=yaml.FullLoader)

# Shorten for quick test
input_data['calculation']['end_time'] = 400.0

print("\nRunning NEM simulation to t=400s (where stratification develops)...")
hd = HydDown(input_data)
hd.run()

# Find a timestep where we have good stratification and liquid is present
idx = np.argmin(np.abs(hd.time_array - 300))

print(f"\n" + "="*80)
print(f"ANALYSIS AT t = {hd.time_array[idx]:.1f} s")
print("="*80)

print(f"\nPhase Temperatures:")
print(f"  T_gas:    {hd.T_gas[idx]:.2f} K ({hd.T_gas[idx]-273.15:.2f} °C)")
print(f"  T_liquid: {hd.T_liquid[idx]:.2f} K ({hd.T_liquid[idx]-273.15:.2f} °C)")
print(f"  T_wall (wetted): {hd.T_vessel_wetted[idx]:.2f} K ({hd.T_vessel_wetted[idx]-273.15:.2f} °C)")
print(f"  ΔT (gas-liquid): {hd.T_gas[idx] - hd.T_liquid[idx]:.2f} K")
print(f"  ΔT (wall-liquid): {hd.T_vessel_wetted[idx] - hd.T_liquid[idx]:.2f} K")

P = hd.P[idx]

print(f"\nPressure: {P/1e5:.2f} bar")

# ============================================================================
# CURRENT IMPLEMENTATION (using equilibrium fluid object)
# ============================================================================
print(f"\n" + "="*80)
print("CURRENT IMPLEMENTATION: Properties used in Rohsenow correlation")
print("="*80)

print("\nIn hdclass.py around line 1250-1256:")
print("  hiw = tp.h_inside_wetted(")
print("      L, T_wall, T_fluid,")
print("      transport_fluid_wet,  # <-- fluid parameter")
print("      self.fluid)           # <-- master_fluid parameter (EQUILIBRIUM!)")

# Recreate what's happening in current code
hd.fluid.update(CP.PT_INPUTS, P, hd.T_fluid[idx])

print(f"\nself.fluid (equilibrium mixed fluid) at P={P/1e5:.2f} bar, T_fluid={hd.T_fluid[idx]:.2f} K:")
print(f"  Quality Q: {hd.fluid.Q():.4f}")
print(f"  Phase: {hd.fluid.phase()}")

# For equilibrium, get saturated properties at the pressure
# (equilibrium fluid might be superheated, so update to saturation)
hd.fluid.update(CP.PQ_INPUTS, P, 0.5)  # Two-phase at 50% quality

print(f"\nSaturated properties from self.fluid (EQUILIBRIUM) at P={P/1e5:.2f} bar:")
rhol_eq = hd.fluid.saturated_liquid_keyed_output(CP.iDmass)
rhog_eq = hd.fluid.saturated_vapor_keyed_output(CP.iDmass)
h_liq_eq = hd.fluid.saturated_liquid_keyed_output(CP.iHmass)
h_vap_eq = hd.fluid.saturated_vapor_keyed_output(CP.iHmass)
Hvap_eq = h_vap_eq - h_liq_eq
sigma_eq = hd.fluid.surface_tension()

print(f"  rhol (saturated liquid density): {rhol_eq:.2f} kg/m³")
print(f"  rhog (saturated vapor density):  {rhog_eq:.2f} kg/m³")
print(f"  Hvap (latent heat):              {Hvap_eq/1000:.2f} kJ/kg")
print(f"  sigma (surface tension):         {sigma_eq*1000:.4f} mN/m")

# ============================================================================
# WHAT IT SHOULD BE (using liquid phase object)
# ============================================================================
print(f"\n" + "="*80)
print("WHAT IT SHOULD BE: Properties from liquid phase object")
print("="*80)

print(f"\nself.fluid_liquid at P={P/1e5:.2f} bar, rho_liquid={hd.rho_liquid[idx]:.2f} kg/m³:")

# Update liquid fluid object to its actual state
hd.fluid_liquid.update(CP.DmassUmass_INPUTS, hd.rho_liquid[idx], hd.U_liquid[idx])

print(f"  T_liquid: {hd.fluid_liquid.T():.2f} K")
print(f"  Quality Q: {hd.fluid_liquid.Q():.4f}")
print(f"  Phase: {hd.fluid_liquid.phase()}")

print(f"\nSaturated properties from self.fluid_liquid (LIQUID PHASE):")
rhol_nem = hd.fluid_liquid.saturated_liquid_keyed_output(CP.iDmass)
rhog_nem = hd.fluid_liquid.saturated_vapor_keyed_output(CP.iDmass)
h_liq_nem = hd.fluid_liquid.saturated_liquid_keyed_output(CP.iHmass)
h_vap_nem = hd.fluid_liquid.saturated_vapor_keyed_output(CP.iHmass)
Hvap_nem = h_vap_nem - h_liq_nem
sigma_nem = hd.fluid_liquid.surface_tension()

print(f"  rhol (saturated liquid density): {rhol_nem:.2f} kg/m³")
print(f"  rhog (saturated vapor density):  {rhog_nem:.2f} kg/m³")
print(f"  Hvap (latent heat):              {Hvap_nem/1000:.2f} kJ/kg")
print(f"  sigma (surface tension):         {sigma_nem*1000:.4f} mN/m")

# ============================================================================
# COMPARISON
# ============================================================================
print(f"\n" + "="*80)
print("COMPARISON: Equilibrium vs NEM Liquid Phase Properties")
print("="*80)

print(f"\n{'Property':<30} {'Equilibrium':<15} {'NEM Liquid':<15} {'Difference':<15}")
print("-"*75)
print(f"{'rhol [kg/m³]':<30} {rhol_eq:<15.2f} {rhol_nem:<15.2f} {rhol_nem-rhol_eq:<15.2f}")
print(f"{'rhog [kg/m³]':<30} {rhog_eq:<15.2f} {rhog_nem:<15.2f} {rhog_nem-rhog_eq:<15.2f}")
print(f"{'Hvap [kJ/kg]':<30} {Hvap_eq/1000:<15.2f} {Hvap_nem/1000:<15.2f} {(Hvap_nem-Hvap_eq)/1000:<15.2f}")
print(f"{'sigma [mN/m]':<30} {sigma_eq*1000:<15.4f} {sigma_nem*1000:<15.4f} {(sigma_nem-sigma_eq)*1000:<15.4f}")

# ============================================================================
# IMPACT ON ROHSENOW CORRELATION
# ============================================================================
print(f"\n" + "="*80)
print("IMPACT ON ROHSENOW HEAT TRANSFER COEFFICIENT")
print("="*80)

# Get transport properties at liquid temperature
T_film = (hd.T_liquid[idx] + hd.T_vessel_wetted[idx]) / 2
hd.transport_fluid_wet.update(CP.PT_INPUTS, P, T_film)

kl = hd.transport_fluid_wet.conductivity()
mul = hd.transport_fluid_wet.viscosity()
Cpl = hd.transport_fluid_wet.cpmass()

from ht import Rohsenow

# Current implementation (equilibrium properties)
h_boil_eq = Rohsenow(
    rhol=rhol_eq,
    rhog=rhog_eq,
    mul=mul,
    kl=kl,
    Cpl=Cpl,
    Hvap=Hvap_eq,
    sigma=sigma_eq,
    Te=max(hd.T_vessel_wetted[idx] - hd.T_liquid[idx], 0),
    Csf=0.013,
    n=1.7,
)

# Corrected (liquid phase properties)
h_boil_nem = Rohsenow(
    rhol=rhol_nem,
    rhog=rhog_nem,
    mul=mul,
    kl=kl,
    Cpl=Cpl,
    Hvap=Hvap_nem,
    sigma=sigma_nem,
    Te=max(hd.T_vessel_wetted[idx] - hd.T_liquid[idx], 0),
    Csf=0.013,
    n=1.7,
)

print(f"\nRohsenow heat transfer coefficient:")
print(f"  Current (equilibrium):  h_boil = {h_boil_eq:.1f} W/m²K")
print(f"  Corrected (NEM liquid): h_boil = {h_boil_nem:.1f} W/m²K")
print(f"  Difference:             {h_boil_nem - h_boil_eq:.1f} W/m²K ({100*(h_boil_nem-h_boil_eq)/h_boil_eq:.1f}%)")

print(f"\n" + "="*80)
print("CONCLUSION")
print("="*80)
print("""
The Rohsenow nucleate boiling correlation is currently using saturated properties
from the EQUILIBRIUM fluid object (self.fluid) instead of the LIQUID PHASE object
(self.fluid_liquid).

For NEM, the liquid phase can be at a different temperature than saturation,
so the saturated properties should come from the liquid phase object evaluated
at its actual pressure and density/internal energy state.

KEY ISSUE in hdclass.py lines 1250-1256 and 1282-1288:
  hiw = tp.h_inside_wetted(
      L, T_wall, T_fluid,
      self.transport_fluid_wet,  # OK - transport properties
      self.fluid)                # WRONG - should be self.fluid_liquid for NEM!

RECOMMENDATION:
Add NEM-specific logic to pass self.fluid_liquid when calculating wetted heat
transfer coefficient for the liquid phase.
""")

print(f"\n" + "="*80)
