"""Regenerate the CARDICE reconciled validation figures and overlay them on the INERIS E0x
1 Hz data (validation/reconciled/_data/). One 2x2 figure per test:

  * Pressure        - measured "P sphere" vs model P
  * Inventory       - measured CO2 mass ("Masse", final value = retained dry ice) vs model
                      total / dry ice / liquid / gas
  * Internal temp   - TWO measured bands: upper 3 TCs (gas) and lower 3 TCs (liquid/solid),
                      vs model gas and liquid/solid temperatures
  * Inside-wall     - TWO measured bands: upper 3 (gas wall) and lower 4 (wetted/dry-ice wall),
                      vs model gas-wall and wetted-wall temperatures

Time-shift: CARDICE blowdowns are slow (tiny orifice, huge sphere) and some liquid tests barely
decline for hours before the rapid drop, so a fixed-drop onset mis-fires badly. Each measured
trace is instead anchored to the MODEL at a mid-blowdown reference pressure (P fallen ~15% of its
span), which is robust for both the gas and liquid tests.

Model physics (branch co2-release-hem): Churchill-Chu gas-wall convection, Cooper plateau boiling
into the liquid, dry-ice wall tracking the sublimation temperature, gas<->solid coupling over the
mass-derived ice-bed interface area. Run from the repo root:  python scripts/validate_cardice_reconciled.py
"""
import glob, os
import pandas as pd, numpy as np, yaml
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from hyddown.hdclass import HydDown

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + "/"
D = REPO + "validation/reconciled/_data/"
OUT = REPO + "validation/reconciled/"
NAVY = "#002D40"; RED = "#D61F39"; AMBER = "#E6A740"; SLATE = "#82979F"; GREEN = "#118a6b"
TESTMAP = {5: "E05", 6: "E6", 7: "E07", 8: "E08", 9: "E09", 10: "E10"}
TITLE = {5: "T5 gas 20 bar", 6: "T6 liquid 20 bar", 7: "T7 gas 15 bar",
         8: "T8 liquid 15 bar", 9: "T9 gas 10 bar", 10: "T10 liquid 10 bar"}


def find(pfx, kind):
    return glob.glob(D + pfx + "-" + kind + "*")[0]


def rd(pfx, kind):
    return pd.read_excel(find(pfx, kind), header=None, engine="calamine").iloc[1:].apply(pd.to_numeric, errors="coerce")


def merge_across_gaps(t, vals, destep=False, thr=5.0):
    """Merge across acquisition gaps (e.g. T10's 385 s outage at ~2.2 h, stored as NaN rows).

    Drops NaN-time samples (carrying ``vals`` along so they stay aligned), then closes time gaps:
    where a step exceeds ``thr`` x the median step, all subsequent samples are shifted earlier so
    the gap collapses to one normal step. With ``destep=True`` the value step across each gap is
    also removed so the series is continuous - the process evolved during the outage (the fluid
    warmed ~1-2 C over T10's gap), which otherwise reads as a spurious jump. Use destep for the
    fluctuating temperature channels; time-only for monotonic P / mass.

    ``vals`` is 1-D or 2-D (rows = series) sharing ``t``'s length. Returns (t2, vals2)."""
    t = np.asarray(t, dtype=float)
    v = np.array(vals, dtype=float)
    v2 = np.atleast_2d(v)
    keep = np.isfinite(t)
    t = t[keep]; v2 = v2[:, keep]
    if len(t) >= 3:
        dt = np.diff(t)
        med = np.nanmedian(dt[dt > 0])
        isgap = np.isfinite(dt) & (dt > thr * med)
        t = t - np.concatenate([[0.0], np.cumsum(np.where(isgap, dt - med, 0.0))])
        if destep:
            for j in np.where(isgap)[0]:
                step = v2[:, j + 1] - v2[:, j]
                step[~np.isfinite(step)] = 0.0
                v2[:, j + 1:] -= step[:, None]
    return t, (v2 if v.ndim > 1 else v2[0])


def main():
    for t, pfx in TESTMAP.items():
        try:
            P = rd(pfx, "pressures"); M = rd(pfx, "mass-flowrate")
            IT = rd(pfx, "internal-temperatures"); WT = rd(pfx, "inside-wall-temperatures")
            okp = np.isfinite(P.iloc[:, 0].values) & np.isfinite(P.iloc[:, 1].values)
            tP, ps = merge_across_gaps(P.iloc[:, 0].values[okp], P.iloc[:, 1].values[okp])

            hd = HydDown(yaml.safe_load(open(OUT + "CARDICE_test%d.yml" % t)))
            hd.run(disable_pbar=True)
            tm_s = np.asarray(hd.time_array); Pm_bar = np.asarray(hd.P) / 1e5
            tm = tm_s / 3600.0
            ET = tm.max()
            # anchor the measured trace to the model at a mid-blowdown reference pressure
            p0 = np.nanmax(ps[:120]); Pref = p0 - max(1.0, 0.15 * (p0 - 1.0))

            def _cross(tt, pp):
                idx = np.where(pp < Pref)[0]; return tt[idx[0]] if len(idx) else tt[-1]
            shift = _cross(tP, ps) - _cross(tm_s, Pm_bar)   # seconds to subtract from measured time
            sh = lambda tt, _s=shift: (np.asarray(tt) - _s) / 3600.0
            alive = np.asarray(hd.mass_fluid) > 0.02
            cond = np.asarray(hd.m_liquid) + np.asarray(hd.m_solid)

            fig, ax = plt.subplots(2, 2, figsize=(12, 8))
            fig.suptitle("CARDICE %s  (orifice %.0f mm)  -  model (CoolProp / Churchill-Chu) vs INERIS 1 Hz data"
                         % (TITLE[t], hd.D_release * 1000), color=NAVY, fontsize=12)
            # Pressure
            a = ax[0, 0]; a.plot(sh(tP), ps, color=SLATE, lw=1.6, label="measured (P sphere)")
            a.plot(tm, hd.P / 1e5, color=RED, lw=1.6, ls="--", label="model")
            a.set_ylabel("pressure [bar]"); a.set_title("Pressure"); a.legend(fontsize=8); a.grid(alpha=.3)
            # Inventory
            a = ax[0, 1]
            tM, mm = merge_across_gaps(M.iloc[:, 0].values, M.iloc[:, 1].values)
            a.plot(sh(tM), mm, color=SLATE, lw=1.6, label="measured CO2 (Masse)")
            a.plot(tm, hd.mass_fluid, color=RED, lw=1.6, ls="--", label="model total")
            a.plot(tm, hd.m_solid, color="k", lw=1, ls=":", label="dry ice")
            a.plot(tm, hd.m_liquid, color=NAVY, lw=1, ls="-.", label="liquid")
            a.plot(tm, hd.m_gas, color=AMBER, lw=1, label="gas")
            a.set_ylabel("CO2 mass [kg]"); a.set_title("Inventory"); a.legend(fontsize=7.5); a.grid(alpha=.3)
            # Internal temperature - two bands: upper 3 TCs (gas) vs lower 3 (liquid/solid)
            a = ax[1, 0]
            tIT, Tc = merge_across_gaps(IT.iloc[:, 0].values, IT.iloc[:, 1:7].values.T, destep=True)  # 6 x N
            Tlo = Tc[0:3]   # Tc1,2,3 - lower / liquid+solid
            Tup = Tc[3:6]   # Tc4,5,6 - upper / gas
            a.fill_between(sh(tIT), np.nanmin(Tup, 0), np.nanmax(Tup, 0), color=AMBER, alpha=.25, label="meas gas band (upper 3)")
            a.fill_between(sh(tIT), np.nanmin(Tlo, 0), np.nanmax(Tlo, 0), color=SLATE, alpha=.40, label="meas liquid/solid band (lower 3)")
            a.plot(tm, np.where(alive, np.asarray(hd.T_gas) - 273.15, np.nan), color=RED, lw=1.5, ls="--", label="model gas")
            a.plot(tm, np.where(alive & (cond > 0.1), np.asarray(hd.T_liquid) - 273.15, np.nan), color=NAVY, lw=1.5, ls="-.", label="model liquid/solid")
            a.set_ylabel("internal T [C]"); a.set_xlabel("time [h]"); a.set_title("Internal temperature"); a.legend(fontsize=7); a.grid(alpha=.3)
            # Inside-wall - two bands: upper 3 (gas wall) vs lower 4 (wetted/dry-ice wall)
            a = ax[1, 1]
            tWT, Wc = merge_across_gaps(WT.iloc[:, 0].values, WT.iloc[:, 1:8].values.T, destep=True)  # 7 x N
            Wlo = Wc[0:4]   # W1-4 lower / wetted wall
            Wup = Wc[4:7]   # W5-7 upper / gas wall
            a.fill_between(sh(tWT), np.nanmin(Wup, 0), np.nanmax(Wup, 0), color=AMBER, alpha=.20, label="meas gas-wall band (upper 3)")
            a.fill_between(sh(tWT), np.nanmin(Wlo, 0), np.nanmax(Wlo, 0), color=GREEN, alpha=.28, label="meas wetted-wall band (lower 4)")
            a.plot(tm, hd.T_vessel - 273.15, color=RED, lw=1.5, ls="--", label="model gas-wall")
            a.plot(tm, hd.T_vessel_wetted - 273.15, color=NAVY, lw=1.5, ls="-.", label="model wetted wall")
            a.set_ylabel("inside-wall T [C]"); a.set_xlabel("time [h]"); a.set_title("Wall temperature"); a.legend(fontsize=7); a.grid(alpha=.3)

            for row in ax:
                for a in row:
                    a.set_xlim(0, ET)
            fig.tight_layout(rect=[0, 0, 1, 0.96])
            for e in ("pdf", "png"):
                fig.savefig(OUT + "CARDICE_test%d_reconciled.%s" % (t, e), dpi=150)
            plt.close(fig)
            meas_ret = M.iloc[:, 1].values; meas_ret = meas_ret[np.isfinite(meas_ret)][-1]
            print("T%-2d m0=%.0f retained model=%.1f vs meas %.1f | Pref-align shift %.0fs -> CARDICE_test%d_reconciled.pdf"
                  % (t, hd.mass_fluid[0], hd.m_solid[-1], meas_ret, shift, t), flush=True)
        except Exception as ex:
            print("T%d FAILED: %s" % (t, str(ex)[:150]), flush=True)
    print("ALL DONE", flush=True)


if __name__ == "__main__":
    main()
