"""Independent CoolProp + HEM hand-calculation of the choked-throat pressure for CO2.

HEM (homogeneous equilibrium) mass flux along the isentrope from the stagnation state:

    G(P) = rho(P) * sqrt( 2 * (h0 - h(P)) ),   s(P) = s0  (isentropic to the throat)

The choked throat is the P that MAXIMISES G between the back pressure and the stagnation
pressure (equivalently, the minimiser of -G -- this is exactly how HydDown's solver finds
it). This script proves two things with plain CoolProp, independent of HydDown:

  1. SUB-COOLED liquid (P0 > Psat(T0)): the flux keeps rising as the single-phase liquid
     accelerates (Bernoulli-like), then peaks right at the saturation (flashing-inception)
     pressure Psat(T0). So the throat pressure sits at ~Psat, BELOW the vessel pressure,
     and the subcooling raises the rate but not the choke pressure.

  2. SATURATED liquid (P0 = Psat): flashing starts immediately, the two-phase sound speed
     collapses, and the flux maximum sits just below P0 -- so throat ~= stagnation
     (this is the A2/P11 tank observation).

Finally it cross-checks the CoolProp-only result against HydDown's CO2ReleaseModelCP.
"""
import os, sys, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar
import CoolProp.CoolProp as CP

matplotlib.rcParams["font.family"] = "Arial"
matplotlib.rcParams["font.size"] = 10
NAVY, RED, AMBER, SLATE, GREY = "#002D40", "#D61F39", "#E6A740", "#82979F", "#4C4D4E"

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src"))

FLUID = "CO2"
P_TRIPLE = 5.18e5     # Pa
P_BACK = 1.01325e5    # Pa
CD = 1.0
D = 0.20272           # 8" Sch40 ID [m] (for a representative mdot cross-check)
AREA = math.pi * D ** 2 / 4.0

_st = CP.AbstractState("HEOS", FLUID)


def props_PT(P, T):
    _st.update(CP.PT_INPUTS, P, T)
    return _st.hmass(), _st.smass(), _st.rhomass()


def iso_props(P, s0):
    """(h, rho) on the isentrope s = s0 at pressure P (above the triple point)."""
    _st.update(CP.PSmass_INPUTS, P, s0)
    return _st.hmass(), _st.rhomass()


def G_of_P(P, h0, s0):
    h, rho = iso_props(P, s0)
    dh = h0 - h
    return rho * math.sqrt(2.0 * dh) if dh > 0 else 0.0


def choke(P0, h0, s0, p_lo=P_TRIPLE + 2e3):
    """Maximise G over [p_lo, P0] == minimise -G. Returns (P_throat, G_max)."""
    res = minimize_scalar(lambda P: -G_of_P(P, h0, s0),
                          bounds=(p_lo, P0 - 1.0), method="bounded",
                          options={"xatol": 1.0})
    return res.x, -res.fun


def sweep(P0, h0, s0, p_lo=P_TRIPLE + 2e3, n=400):
    Ps = np.linspace(p_lo, P0, n)
    Gs = np.array([G_of_P(P, h0, s0) for P in Ps])
    return Ps, Gs


def report(label, P0, T0, saturated=False):
    if saturated:
        # stagnation is the saturated liquid at P0 (on the dome: use a Q=0 flash, not PT)
        Psat = P0
        T0 = CP.PropsSI("T", "P", P0, "Q", 0, FLUID)
        h0 = CP.PropsSI("Hmass", "P", P0, "Q", 0, FLUID)
        s0 = CP.PropsSI("Smass", "P", P0, "Q", 0, FLUID)
        rho0 = CP.PropsSI("Dmass", "P", P0, "Q", 0, FLUID)
    else:
        Psat = CP.PropsSI("P", "T", T0, "Q", 0, FLUID)
        h0, s0, rho0 = props_PT(P0, T0)
    subcool = P0 - Psat
    Pth, Gmax = choke(P0, h0, s0)
    hth, rhoth = iso_props(Pth, s0)
    Tth = CP.PropsSI("T", "P", Pth, "Smass", s0, FLUID)
    mdot = CD * AREA * Gmax
    print("\n== %s ==" % label)
    print("  stagnation      P0 = %8.3f bar   T0 = %7.2f C" % (P0 / 1e5, T0 - 273.15))
    print("  saturation    Psat = %8.3f bar   (subcooling dP = %6.3f bar)" % (Psat / 1e5, subcool / 1e5))
    print("  rho0 = %7.1f kg/m3   h0 = %9.1f J/kg   s0 = %8.2f J/kg/K" % (rho0, h0, s0))
    print("  --> choke     Pth = %8.3f bar   Tth = %7.2f C   rho_th = %7.1f kg/m3" %
          (Pth / 1e5, Tth - 273.15, rhoth))
    print("      Gmax = %9.1f kg/m2/s   v_th = %6.1f m/s   mdot(8\" Cd1) = %8.1f kg/s" %
          (Gmax, Gmax / rhoth, mdot))
    print("      Pth / Psat = %.4f    (=> throat pressure pinned at the flashing-inception line)" %
          (Pth / Psat))
    return dict(label=label, P0=P0, T0=T0, Psat=Psat, h0=h0, s0=s0, rho0=rho0,
                Pth=Pth, Gmax=Gmax, Tth=Tth, rhoth=rhoth, mdot=mdot)


def crosscheck_hyddown(P0, T0):
    """Same subcooled state through HydDown's CoolProp HEM model."""
    from hyddown.co2_release_cp import CO2ReleaseModelCP
    m = CO2ReleaseModelCP(back_pressure=P_BACK, atm_pressure=P_BACK)
    r = m.dense_leak_rate(T0, P0, CD, AREA, phase="liquid")
    print("\n== HydDown CO2ReleaseModelCP.dense_leak_rate (subcooled) ==")
    print("  Pth = %8.3f bar   Tth = %7.2f C   rho_th = %7.1f kg/m3   v_th = %6.1f m/s   mdot = %8.1f kg/s"
          % (r["P_throat"] / 1e5, r["T_throat"] - 273.15, r.get("rho_throat", float("nan")),
             r.get("v_throat", float("nan")), r["mdot"]))


def main():
    # sub-cooled liquid: the P5 line condition (17.5 bara, -26 C)
    sub = report("Sub-cooled liquid  (P5: 17.5 bara, -26 C)", 17.5e5, -26.0 + 273.15)
    # sub-cooled, more strongly (P4-like, colder margin) to show Pth still pinned at Psat
    sub2 = report("Sub-cooled liquid  (deeper: 25 bara, -26 C)", 25.0e5, -26.0 + 273.15)
    # saturated liquid: the P11 tank stagnation (16.2 bara, saturated)
    sat = report("Saturated liquid   (P11 tank: 16.2 bara, Q=0)", 16.2e5, None, saturated=True)

    try:
        crosscheck_hyddown(17.5e5, -26.0 + 273.15)
    except Exception as e:
        print("\n(HydDown cross-check skipped: %s)" % repr(e)[:80])

    # figure: G(P) for the sub-cooled and saturated stagnation states
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.5))
    for a, c, col in ((ax[0], sub, NAVY), (ax[1], sat, RED)):
        Ps, Gs = sweep(c["P0"], c["h0"], c["s0"])
        a.plot(Ps / 1e5, Gs, color=col, lw=2.0)
        a.axvline(c["Psat"] / 1e5, color=SLATE, ls="--", lw=1.2, label="$P_{sat}(T_0)$")
        a.axvline(c["P0"] / 1e5, color=GREY, ls=":", lw=1.2, label="$P_0$ (stagnation)")
        a.plot(c["Pth"] / 1e5, c["Gmax"], "o", color=AMBER, ms=8, label="choke (max $G$)")
        a.set_xlabel("Throat pressure [bar]"); a.set_ylabel("HEM mass flux $G$ [kg/m$^2$/s]")
        a.grid(True, color=SLATE, alpha=0.25, lw=0.5); a.legend(frameon=False, fontsize=8)
    ax[0].set_title("(a) Sub-cooled liquid: choke at $P_{sat}$, below $P_0$", color=NAVY)
    ax[1].set_title("(b) Saturated liquid: choke just below $P_0$", color=NAVY)
    fig.suptitle("HEM choked-throat pressure for CO$_2$ (CoolProp isentrope)",
                 fontsize=12, fontweight="bold", color=NAVY)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    out_dir = ("C:\\Users\\AndersAndreasen\\OneDrive - ORS\\Projects 2026 - "
               "120.691_D_CarbonCuts_Ruby Development Project\\01 Working area\\CFD\\Source terms")
    os.makedirs(out_dir, exist_ok=True)
    pdf = os.path.join(out_dir, "HEM_choke_handcalc.pdf")
    fig.savefig(pdf); plt.close(fig)
    print("\nwrote %s" % pdf)


if __name__ == "__main__":
    main()
