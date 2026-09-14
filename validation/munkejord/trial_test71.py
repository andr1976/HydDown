"""
Trial: HydDown tcPR CO2 engine vs SINTEF/Munkejord Test 71 (8 mm no-riser, 122.6 bar, 25.2 C).
Dense-phase start -> two-phase (CoolProp) -> triple point / dry ice (thermopack CO2ReleaseModel).
Lumped SS316 wall heat model. Cd = 1 (smooth nozzle, frictionless HEM, per paper).
"""
import sys, math
import numpy as np
sys.path.insert(0, r"C:/Users/AndersAndreasen/Documents/GitHub/HydDown/src")
from CoolProp.CoolProp import PropsSI
from hyddown.co2_release import CO2ReleaseModel

# ---- geometry / material ----
D = 0.273; H = 1.0
V = math.pi/4*D**2*H                      # 0.05853 m3
d_th = 0.008; A = math.pi/4*d_th**2       # 8 mm nozzle throat
Cd = 1.0
P0 = 122.6e5; T0 = 25.2+273.15
p_back = 0.90e5                            # Trondheim ambient (~0.87 bar measured)
m_steel = 220.0; cp_steel = 500.0         # lumped effective wall
A_in = math.pi*D*H + math.pi/4*D**2       # inner surface ~0.916 m2
A_side = math.pi*D*H
h_boil = 3000.0; h_conv = 15.0            # wetted boiling / dry free-conv HTC
t_strat = 17.0                            # stratification time (paper)

import os
EOS = os.environ.get("EOS", "tcPR")
if EOS.upper() == "CP":
    from hyddown.co2_release_cp import CO2ReleaseModelCP
    rm = CO2ReleaseModelCP(back_pressure=p_back, atm_pressure=p_back)
elif EOS.lower().replace("-", "") == "tcpr":
    rm = CO2ReleaseModel(back_pressure=p_back, atm_pressure=p_back, eos="tcPR")
else:
    # build with a thermopack multiparameter EOS (GERG2008 / MEOS), bypassing the
    # tcPR-only __init__ but running HydDown's own init routines
    from thermopack.multiparameter import multiparam
    rm = CO2ReleaseModel.__new__(CO2ReleaseModel)
    rm.eos = multiparam("CO2", EOS); rm.eos.init_solid("CO2")
    rm.z = np.array([1.0]); rm.LIQ = rm.eos.LIQPH; rm.VAP = rm.eos.VAPPH
    rm.M = rm.eos.compmoleweight(1)/1000.0
    rm.p_back = p_back; rm.p_atm = p_back
    rm.liquid_nonequilibrium = 0.0; rm.liquid_ne_pressure_scaled = False; rm.liquid_ne_pref = None
    rm.P_TRIPLE = 5.18e5
    rm._init_atm_endpoints(); rm._init_triple_point(); rm._init_gas_table()
print(f"[EOS={EOS}] triple P={rm.P_TRIPLE_EOS/1e5:.4f} bar  frost={rm.T_frost:.2f} K")
Mw = rm.M
z = getattr(rm, "z", None); LIQ = getattr(rm, "LIQ", None); VAP = getattr(rm, "VAP", None)
_TP_EOS = getattr(rm, "eos", None) is not None   # True for thermopack, False for CoolProp-only
P_TP = rm.P_TRIPLE_EOS
def gas_enthalpy(T, P):
    """Backend-agnostic vapour enthalpy [J/kg]."""
    if _TP_EOS:
        return float(rm.eos.enthalpy(T, P, z, VAP)[0]) / Mw
    return rm._gasp(T, P)[0]
h_boil = float(os.environ.get("HBOIL", h_boil))
h_conv = float(os.environ.get("HCONV", h_conv))
m_steel = float(os.environ.get("MSTEEL", m_steel))

USE_ROHSENOW = os.environ.get("ROHSENOW", "0") == "1"
# series wall-conduction resistance [m2K/W]; "auto" = (t_wall/2)/k_steel (inner face -> mid-plane)
_rc = os.environ.get("RCOND", "0")
R_COND = (0.0254/2)/16.3 if _rc == "auto" else float(_rc)
from ht import Rohsenow
def rohsenow_h(P, Tw, Tf):
    """HydDown's h_inside_wetted boiling model (ht.Rohsenow, Csf=0.013, n=1.7, cap 3000)."""
    try:
        sigma = PropsSI('surface_tension', 'P', P, 'Q', 0, 'CO2')
    except Exception:
        return 0.0
    if sigma < 1e-6:
        return 0.0
    rhol = PropsSI('Dmass','P',P,'Q',0,'CO2'); rhog = PropsSI('Dmass','P',P,'Q',1,'CO2')
    mul  = PropsSI('V','P',P,'Q',0,'CO2');     kl   = PropsSI('L','P',P,'Q',0,'CO2')
    Cpl  = PropsSI('Cpmass','P',P,'Q',0,'CO2')
    Hvap = PropsSI('Hmass','P',P,'Q',1,'CO2') - PropsSI('Hmass','P',P,'Q',0,'CO2')
    Te   = max(Tw - Tf, 0.0)
    try:
        h = Rohsenow(rhol=rhol, rhog=rhog, mul=mul, kl=kl, Cpl=Cpl,
                     Hvap=Hvap, sigma=sigma, Te=Te, Csf=0.013, n=1.7)
    except Exception:
        return 0.0
    return 0.0 if math.isnan(h) else min(h, 3000.0)

def dense_stagnation(P, T, phase_liq=True):
    """Single dense-phase stagnation at (T,P), in the model's own h/s basis."""
    if _TP_EOS:  # thermopack: molar entropy, thermopack basis
        ph = LIQ if phase_liq else VAP
        h0 = float(rm.eos.enthalpy(T, P, z, ph)[0])/Mw
        s0 = float(rm.eos.entropy(T, P, z, ph)[0])
        v0 = float(rm.eos.specific_volume(T, P, z, ph)[0])
        return h0, s0, Mw/v0
    # CoolProp-only: mass entropy, CoolProp basis
    return (PropsSI("Hmass","T",T,"P",P,"CO2"), PropsSI("Smass","T",T,"P",P,"CO2"),
            PropsSI("Dmass","T",T,"P",P,"CO2"))

# ---- init (CoolProp dense-liquid basis above triple) ----
rho0 = PropsSI('Dmass','P',P0,'T',T0,'CO2'); M = rho0*V
u0 = PropsSI('Umass','P',P0,'T',T0,'CO2'); U = M*u0
Tw = T0
dt = 0.1; tmax = 320.0
n = int(tmax/dt)+1

t_hist=[]; P_hist=[]; T_hist=[]; m_hist=[]; mdot_hist=[]
ms_hist=[]; regime_hist=[]; hw_hist=[]
below = False; t_two=None

for k in range(n):
    t = k*dt
    if not below:
        rho = M/V; u = U/M
        try:
            P = PropsSI('P','Dmass',rho,'Umass',u,'CO2')
            T = PropsSI('T','Dmass',rho,'Umass',u,'CO2')
            Q = PropsSI('Q','Dmass',rho,'Umass',u,'CO2')
        except Exception:
            break
        # hand off to thermopack at/below the triple point
        if P <= P_TP*1.001:
            below = True
            m_l, m_g, U = rm.triple_LG_from_MV(M, V)
            # fall through to below-triple branch next iteration
            regime='handoff'
            # record and continue
        if not below:
            two_phase = (0.0 <= Q <= 1.0)
            if two_phase:
                if t_two is None: t_two = t
                # stratified vapour leaves (ramp liquid fraction 1->0 over t_strat)
                frac_liq = max(1.0 - (t - t_two)/t_strat, 0.0)
                hV = PropsSI('Hmass','P',P,'Q',1,'CO2')
                hL = PropsSI('Hmass','P',P,'Q',0,'CO2')
                h_out = frac_liq*hL + (1-frac_liq)*hV
                r = rm.hem_rate(P, 'gas', Cd, A)          # vapour HEM (thermopack sat-vap stag)
                mdot = r['mdot']
                # liquid volume fraction for wetted area
                rhoL=PropsSI('Dmass','P',P,'Q',0,'CO2'); rhoV=PropsSI('Dmass','P',P,'Q',1,'CO2')
                x=min(max(Q,0),1); vfl=(1-x)/rhoL/((1-x)/rhoL+x/rhoV); alpha_l=vfl
            else:
                # single dense phase leaves
                h_out = PropsSI('Hmass','P',P,'T',T,'CO2')
                h0,s0,rd = dense_stagnation(P, T, phase_liq=True)
                r = rm.hem_rate_from_stagnation(h0,s0,rd,T,P,Cd,A)
                mdot = r['mdot']; alpha_l = 1.0
            # lumped wall heat (Rohsenow boiling HTC if enabled, else fixed h_boil)
            h_wet = rohsenow_h(P, Tw, T) if (USE_ROHSENOW and two_phase) else h_boil
            if R_COND > 0 and h_wet > 0:      # series wall-conduction resistance (lumped proxy for thermesh)
                h_wet = 1.0/(1.0/h_wet + R_COND)
            hw_hist.append(h_wet)
            A_wet = alpha_l*A_side + math.pi/4*D**2
            Q_wall = h_wet*A_wet*(Tw-T) + h_conv*(A_in-A_wet)*(Tw-T)
            # update
            M -= mdot*dt; U += (Q_wall - mdot*h_out)*dt
            Tw -= Q_wall/(m_steel*cp_steel)*dt
            t_hist.append(t); P_hist.append(P/1e5); T_hist.append(T-273.15)
            m_hist.append(M); mdot_hist.append(mdot); ms_hist.append(0.0)
            regime_hist.append('two' if two_phase else 'dense')
            continue
    # ---- below / at triple point (thermopack) ----
    st = rm.vessel_state_below_triple(M, U, V)
    reg = st['regime']; P = st['P']; T = st['T']
    m_s = max(st['m_s'], 0.0)
    if reg == 'above':
        below = False   # warmed back up (shouldn't happen here)
        continue
    T_gas = T
    h_out = gas_enthalpy(T_gas, P)
    mdot = rm.gas_leak_rate(T_gas, P, Cd, A)['mdot']
    A_wet = math.pi/4*D**2
    Q_wall = h_conv*(A_in-A_wet)*(Tw-T)
    M -= mdot*dt; U += (Q_wall - mdot*h_out)*dt
    Tw -= Q_wall/(m_steel*cp_steel)*dt
    t_hist.append(t); P_hist.append(P/1e5); T_hist.append(T-273.15)
    m_hist.append(M); mdot_hist.append(mdot); ms_hist.append(m_s)
    regime_hist.append(reg)
    if mdot < 1e-6 and P <= p_back*1.02:
        break

t=np.array(t_hist); P=np.array(P_hist); Tf=np.array(T_hist)
m=np.array(m_hist); ms=np.array(ms_hist)
if hw_hist:
    hw=np.array(hw_hist); print(f"boiling HTC (two-phase): min {hw.min():.0f}  mean {hw.mean():.0f}  max {hw.max():.0f} W/m2K  [ROHSENOW={USE_ROHSENOW}]")
print(f"steps={len(t)}  end t={t[-1]:.1f}s  end P={P[-1]:.2f}bar  end m={m[-1]:.2f}kg  retained dry ice={ms[-1]:.2f}kg")
print(f"two-phase onset t={t_two}  min Tf={Tf.min():.1f}C")
np.savez(r"trial71_out.npz", t=t,P=P,Tf=Tf,m=m,ms=ms)
