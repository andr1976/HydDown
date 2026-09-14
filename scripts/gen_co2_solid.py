"""Generate solid-CO2 property tables (offline, from thermopack GERG), re-referenced into
CoolProp's h/s basis via the physical sublimation latent heat. Writes src/hyddown/co2_solid.py
with embedded numpy arrays + accessors so the RUNTIME needs only CoolProp + numpy (no thermopack)."""
import os
import numpy as np
from scipy import optimize
import CoolProp as CP
from CoolProp.CoolProp import AbstractState
from thermopack.multiparameter import multiparam

# output: src/hyddown/co2_solid.py, relative to this script (scripts/gen_co2_solid.py)
_OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src", "hyddown", "co2_solid.py")

g = multiparam("CO2", "GERG2008"); g.init_solid("CO2")
z = np.array([1.0]); M = g.compmoleweight(1) / 1000.0
LIQ, VAP = g.LIQPH, g.VAPPH

def gsol(T, P):
    return (float(g.solid_enthalpy(T, P, z)[0]) - T * float(g.solid_entropy(T, P, z)[0])) / M
def gflu(T, P, ph):
    return (float(g.enthalpy(T, P, z, ph)[0]) - T * float(g.entropy(T, P, z, ph)[0])) / M

# tcPR/GERG self-consistent triple point (g_solid == g_liquid on the bubble line)
Ttp = optimize.brentq(lambda T: gsol(T, float(g.bubble_pressure(T, z)[0]))
                      - gflu(T, float(g.bubble_pressure(T, z)[0]), LIQ), 208, 220)
Ptp = float(g.bubble_pressure(Ttp, z)[0])

def Psub(T):  # sublimation pressure: g_vapour == g_solid
    return optimize.brentq(lambda P: gflu(T, P, VAP) - gsol(T, P), 1.0, Ptp * 1.0001)

# CoolProp forced-gas vapour (extrapolates below triple) for the basis anchor
gas = AbstractState("HEOS", "CO2"); gas.specify_phase(CP.iphase_gas)
def cp_gas(P, T):
    gas.update(CP.PT_INPUTS, P, T); return gas.hmass(), gas.smass()

T_grid = np.linspace(150.0, Ttp, 160)
P_sub = np.empty_like(T_grid); v_s = np.empty_like(T_grid)
h_s = np.empty_like(T_grid); s_s = np.empty_like(T_grid)
for i, T in enumerate(T_grid):
    P = Psub(T); P_sub[i] = P
    hs_tp = float(g.solid_enthalpy(T, P, z)[0]) / M
    ss_tp = float(g.solid_entropy(T, P, z)[0]) / M
    v_s[i] = float(g.solid_volume(T, P, z)[0]) / M
    hg_tp = float(g.enthalpy(T, P, z, VAP)[0]) / M
    sg_tp = float(g.entropy(T, P, z, VAP)[0]) / M
    L_sub = hg_tp - hs_tp           # basis-independent latent heat of sublimation
    ds_sub = sg_tp - ss_tp
    hg_cp, sg_cp = cp_gas(P, T)     # CoolProp vapour at (T, P_sub)
    h_s[i] = hg_cp - L_sub          # solid enthalpy in CoolProp basis
    s_s[i] = sg_cp - ds_sub         # solid entropy   in CoolProp basis

# sanity: at the triple point solid must be colder (lower h) than CoolProp sat-liquid
from CoolProp.CoolProp import PropsSI
hl_tp = PropsSI("Hmass", "P", Ptp, "Q", 0, "CO2")
hg_tp = PropsSI("Hmass", "P", Ptp, "Q", 1, "CO2")
print(f"triple T={Ttp:.4f}K P={Ptp/1e5:.5f}bar")
print(f"@triple h_solid={h_s[-1]:.0f}  h_liq(CP)={hl_tp:.0f}  h_gas(CP)={hg_tp:.0f}  (solid<liq<gas: {h_s[-1]<hl_tp<hg_tp})")
print(f"L_fus@triple = h_liq-h_solid = {hl_tp-h_s[-1]:.0f} J/kg   L_sub = {hg_tp-h_s[-1]:.0f} J/kg")

# write the module
def arr(a): return "np.array([" + ", ".join(f"{x:.8g}" for x in a) + "])"
src = f'''"""Solid CO2 (dry-ice) property tables for HydDown's CoolProp-only CO2 release model.

Generated OFFLINE by scratchpad/gen_solid.py from thermopack GERG-2008 + the Hammer solid-CO2
model, then re-referenced into CoolProp's enthalpy/entropy basis via the physical latent heat of
sublimation at each temperature. The RUNTIME needs only CoolProp + numpy - NOT thermopack.

All quantities are mass-specific, tabulated along the sublimation line P_sub(T) (the solid is
nearly incompressible, so the weak pressure dependence off the line is neglected). Values are in
CoolProp's HEOS::CO2 reference basis, so they combine directly with CoolProp fluid properties.
"""
import numpy as np

T_TRIPLE = {Ttp:.6f}      # K   (GERG self-consistent triple point)
P_TRIPLE = {Ptp:.6f}      # Pa

_T   = {arr(T_grid)}
_PSUB = {arr(P_sub)}
_VS  = {arr(v_s)}
_HS  = {arr(h_s)}
_SS  = {arr(s_s)}


def P_sub(T):
    """Sublimation (solid+vapour equilibrium) pressure [Pa] at temperature T [K]."""
    return float(np.interp(T, _T, _PSUB))


def T_sub(P):
    """Sublimation temperature [K] at pressure P [Pa] (inverse of P_sub)."""
    return float(np.interp(P, _PSUB, _T))


def v_solid(T):
    """Solid CO2 specific volume [m3/kg] at temperature T [K]."""
    return float(np.interp(T, _T, _VS))


def h_solid(T):
    """Solid CO2 specific enthalpy [J/kg] (CoolProp basis) at temperature T [K]."""
    return float(np.interp(T, _T, _HS))


def s_solid(T):
    """Solid CO2 specific entropy [J/kg/K] (CoolProp basis) at temperature T [K]."""
    return float(np.interp(T, _T, _SS))


def u_solid(T, P):
    """Solid CO2 specific internal energy [J/kg]: u = h - P v."""
    return h_solid(T) - P * v_solid(T)


def cp_solid(T):
    """Solid CO2 specific heat [J/kg/K] = d h_solid / dT."""
    return float(np.interp(T, 0.5 * (_T[:-1] + _T[1:]), np.diff(_HS) / np.diff(_T)))
'''
with open(_OUT, "w") as f:
    f.write(src)
print("wrote src/hyddown/co2_solid.py")
print(f"cp_solid@200K ~ {float(np.interp(200,0.5*(T_grid[:-1]+T_grid[1:]),np.diff(h_s)/np.diff(T_grid))):.1f} J/kgK")
