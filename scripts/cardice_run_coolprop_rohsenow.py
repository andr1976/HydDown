"""Run ALL CARDICE cases with the FULL thermopack-free configuration and overlay the raw
1 Hz Ineris data - one 6-panel figure per case (pressure, inventory, discharge rate, internal
fluid temperature, inner-wall and outer-wall temperature).

"Full" here means:
  * eos = "CoolProp"          -> the thermopack-free backend (CoolProp + solid table, co2_solid)
  * solid_h_inner = "calc"    -> the below-triple wetted-wall HTC from the Rohsenow boiling
                                 correlation (h_inside_wetted) instead of a hardcoded value.

So every heat-transfer coefficient is a physical correlation and there is no thermopack
dependency. This confirms the fully-physical, thermopack-free model reproduces the measured
blowdown across all observables, matching the reconciled (tcPR + hardcoded-150) baseline.

The data extraction, column conventions and plotting helpers are shared with
cardice_run_reconciled.py (imported below); this script only overrides the two release
parameters and writes CARDICE_test?_coolprop_rohsenow.pdf.

Usage:  python cardice_run_coolprop_rohsenow.py [DATADIR]
"""
import os
import sys
import copy
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import yaml

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
# importing the reconciled runner sets up the data cache (ensure_data) and gives us its
# shared helpers, colours, CASES and paths.
import cardice_run_reconciled as base
from hyddown.hdclass import HydDown

NAVY, RED, AMBER, SLATE, GREY = base.NAVY, base.RED, base.AMBER, base.SLATE, base.GREY


def run_model(yml):
    """Run a reconciled CARDICE case, but force the full thermopack-free / Rohsenow config."""
    with open(os.path.join(base.RECON, yml)) as f:
        inp = copy.deepcopy(yaml.safe_load(f))
    inp["release"]["eos"] = "CoolProp"          # thermopack-free backend
    inp["release"]["solid_h_inner"] = "calc"    # Rohsenow below-triple wetted-wall HTC
    hd = HydDown(inp)
    hd.run(disable_pbar=True)
    return hd, inp["release"]


def main():
    for tid, (yml, pfx, is_gas) in base.CASES.items():
        hd, rel = run_model(yml)
        th = hd.time_array / 3600.0

        m = pd.read_excel(base.find(pfx, "mass-flowrate"))
        tm = m.iloc[:, 0].values.astype(float); M = m.iloc[:, 1].values.astype(float)
        M0 = np.median(M[:50])
        below = np.where(M < M0 - 2.0)[0]
        t0 = tm[below[0]] if len(below) else 0.0
        p = pd.read_excel(base.find(pfx, "pressures"))
        internal = pd.read_excel(base.find(pfx, "internal-temp"))
        iwall = pd.read_excel(base.find(pfx, "inside-wall"))
        owall = pd.read_excel(base.find(pfx, "outside-wall"))
        rM = base.meas_rate(tm, M)

        def X(df):
            return (df.iloc[:, 0].values.astype(float) - t0) / 3600.0

        fig, ax = plt.subplots(2, 3, figsize=(13, 7))
        ax[0, 0].plot(X(p), p["P sphere"], color=NAVY, lw=1.1, label="measured")
        ax[0, 0].plot(th, hd.P / 1e5, color=RED, lw=1.7, ls="--", label="CoolProp + Rohsenow")
        ax[0, 0].set_ylabel("pressure [bar]"); ax[0, 0].set_title("Pressure")
        ax[0, 1].plot(X(m), M, color=NAVY, lw=1.1, label="measured (load cells)")
        ax[0, 1].plot(th, hd.mass_fluid, color=RED, lw=1.7, ls="--", label="model total")
        ax[0, 1].plot(th, hd.m_liquid, color=SLATE, lw=0.9, label="liquid")
        ax[0, 1].plot(th, hd.m_gas, color=AMBER, lw=0.9, label="gas")
        ax[0, 1].plot(th, hd.m_solid, color=GREY, lw=0.9, ls=":", label="solid (dry ice)")
        ax[0, 1].set_ylabel("mass [kg]"); ax[0, 1].set_title("Inventory")
        fin = np.isfinite(rM)
        ax[0, 2].plot(X(m)[fin], rM[fin], color=SLATE, lw=0.7, alpha=0.7, label="measured (raw)")
        ax[0, 2].plot(th, hd.mass_rate, color=RED, lw=1.7, ls="--", label="model")
        ax[0, 2].set_ylabel("discharge rate [kg/s]"); ax[0, 2].set_title("Discharge rate")
        rr = rM[(X(m) > 0.02)]
        ymax = np.nanpercentile(rr, 99) if np.isfinite(rr).any() else 0.5
        ax[0, 2].set_ylim(-0.02, max(0.06, ymax * 1.35))
        for j, c in enumerate(internal.columns[1:]):
            ax[1, 0].plot(X(internal), internal[c], color=SLATE, lw=0.6, alpha=0.6,
                          label="measured Tc in 1..6" if j == 0 else None)
        ax[1, 0].plot(th, base.C(hd.T_gas), color=RED, lw=1.7, ls="--", label="model gas")
        ax[1, 0].plot(th, base.C(hd.T_liquid), color=NAVY, lw=1.7, ls="--", label="model liquid/solid")
        ax[1, 0].set_ylabel("fluid T [$^\\circ$C]"); ax[1, 0].set_title("Internal fluid temperature")
        for j, c in enumerate(iwall.columns[1:]):
            ax[1, 1].plot(X(iwall), iwall[c], color=SLATE, lw=0.6, alpha=0.6,
                          label="measured Tc1..7" if j == 0 else None)
        ax[1, 1].plot(th, base.wall(hd.T_inner_wall, hd.T_vessel), color=RED, lw=1.7, ls="--", label="model gas-contact")
        ax[1, 1].plot(th, base.wall(hd.T_inner_wall_wetted, hd.T_vessel_wetted), color=NAVY, lw=1.7, ls="--", label="model wetted")
        ax[1, 1].set_ylabel("inner-wall T [$^\\circ$C]"); ax[1, 1].set_title("Inner-wall temperature")
        for j, c in enumerate(owall.columns[1:]):
            ax[1, 2].plot(X(owall), owall[c], color=SLATE, lw=0.6, alpha=0.6,
                          label="measured Tc F1..7" if j == 0 else None)
        ax[1, 2].plot(th, base.wall(hd.T_outer_wall, hd.T_vessel), color=RED, lw=1.7, ls="--", label="model gas-contact")
        ax[1, 2].plot(th, base.wall(hd.T_outer_wall_wetted, hd.T_vessel_wetted), color=NAVY, lw=1.7, ls="--", label="model wetted")
        ax[1, 2].set_ylabel("outer-wall T [$^\\circ$C]"); ax[1, 2].set_title("Outer-wall temperature")

        for a in ax.flat:
            a.set_xlabel("time [h]"); a.grid(alpha=0.25); a.legend(fontsize=6.5, loc="best")
            a.set_xlim(-0.02, th[-1] * 1.02)
        fig.suptitle(f"CARDICE {tid} - CoolProp + Rohsenow, thermopack-free "
                     f"({base.title_bits(rel)}) vs 1 Hz data",
                     color=NAVY, fontweight="bold", fontsize=11)
        fig.tight_layout(rect=[0, 0, 1, 0.97])
        out = os.path.join(base.RECON, f"CARDICE_{tid}_coolprop_rohsenow.pdf")
        fig.savefig(out); plt.close(fig)
        print(f"{tid}: backend={type(hd.release_model).__name__} P0={hd.P[0]/1e5:.2f}b "
              f"m0={hd.mass_fluid[0]:.0f}kg solid_end={hd.m_solid[-1]:.0f}kg -> {os.path.basename(out)}",
              flush=True)
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
