"""Validate the opt-in 2-D conjugate wall (hyddown.wall2d) on the SINTEF dense-phase tests.

Runs each case with the 2-D wall (full stepped domain: bottom plate + shell + flange + lid) and,
for reference, with the lumped two-node wall, and overlays both against the measured inside- and
outer-wall bands. Writes one 2x2 figure per test to validation/munkejord/munke_<name>_2dwall.pdf
(+ .png for inspection); does NOT overwrite the production munke_<name>.pdf.

Run from the repo root:  python scripts/validate_2dwall.py
"""
import os, sys, zipfile, io, time
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "src"))
sys.path.insert(0, os.path.join(REPO, "scripts"))
from hyddown.hdclass import HydDown
import validate_munkejord as vm

OUT = os.path.join(REPO, "validation", "munkejord"); os.makedirs(OUT, exist_ok=True)
NAVY = "#002D40"; RED = "#D61F39"; AMBER = "#E6A740"; SLATE = "#82979F"; GREEN = "#118a6b"; DGREY = "#4C4D4E"
OWALL = ["TT101", "TT111", "TT121", "TT131", "TT141", "TT151"]   # outer wall, 6 heights

# name, P0_bar, T0_C, nozzle_mm, is_riser, end_time
CASES = [("Exp53", 119.5, 24.4, 8.0, True, 100.0),
         ("Exp72", 119.0, 24.9, 6.5, False, 400.0)]


def read_outer(z, name, ET):
    xlsx = [n for n in z.namelist() if n.split("/")[-1].startswith(name + "_")][0]
    D = pd.read_excel(io.BytesIO(z.read(xlsx)), sheet_name="Data", header=None, engine="calamine")
    names = [str(x) for x in D.iloc[0].tolist()]; col = {n: i for i, n in enumerate(names)}
    dat = D.iloc[3:].reset_index(drop=True)

    def raw(tn, vn):
        t = pd.to_numeric(dat[col[tn]], errors="coerce").values
        v = pd.to_numeric(dat[col[vn]], errors="coerce").values
        m = np.isfinite(t) & np.isfinite(v); return t[m], v[m]
    tp, p = raw("tPT162", "PT162"); p0 = np.nanmax(p[:200])
    below = np.where(p < p0 - 5.0)[0]; t0 = tp[below[0]] if len(below) else 0.0
    tmax = min(tp.max(), ET + t0); g = np.arange(0.0, tmax - t0, 0.5)
    owall = np.vstack([np.interp(g, raw("tTT", c)[0] - t0, raw("tTT", c)[1]) for c in OWALL])
    return g, owall


def run(P0, T0C, noz, riser, ET, wall2d):
    T0 = T0C + 273.15
    vessel = {"length": 1.0, "diameter": 0.273, "thickness": 0.0254, "heat_capacity": 500,
              "density": 7950, "orientation": "vertical", "type": "Flat-end"}
    if wall2d:
        vessel.update({"wall_model": "2d", "bottom_thickness": 0.050,
                       "flange_thickness": 0.083, "lid_thickness": 0.080, "lid_diameter": 0.580})
    d = {
        "vessel": vessel,
        "initial": {"temperature": T0, "pressure": P0 * 1e5, "fluid": "CO2"},
        "calculation": {"type": "energybalance", "time_step": 0.1, "end_time": ET,
                        "non_equilibrium": True, "h_gas_liquid": "calc_two_sided"},
        "valve": {"flow": "discharge", "type": "none", "back_pressure": 101325.},
        "release": {"type": "liquid" if riser else "gas", "diameter": noz / 1000., "discharge_coef": 1.0,
                    "liquid_nonequilibrium": 0.0, "back_pressure": 101325., "atm_pressure": 101325.,
                    "eos": "CoolProp", "solid_in_vessel": True, "solid_h_inner": "cooper",
                    "solid_h_gas_wall": "churchill", "solid_h_gas_liquid": 0.0, "solid_h_gas_solid": 3.0,
                    "discharge_location": 0.009 if riser else 0.0},
        "heat_transfer": {"type": "specified_h", "temp_ambient": T0, "h_outer": 0.0, "h_inner": "churchill"},
    }
    for dt in (0.1, 0.02):
        d["calculation"]["time_step"] = dt
        try:
            t0 = time.time(); hd = HydDown(d); hd.run(disable_pbar=True); return hd, time.time() - t0
        except Exception:
            if dt == 0.02:
                raise
    return hd, 0.0


def figure(name, P0, T0C, noz, riser, ET):
    z = zipfile.ZipFile(os.path.join(REPO, "background", "19589510.zip"))
    g, Pavg, W, fluid, iwall = vm.measured(z, name, ET)
    _, owall = read_outer(z, name, ET)
    hd2, t2 = run(P0, T0C, noz, riser, ET, wall2d=True)
    hd0, t0 = run(P0, T0C, noz, riser, ET, wall2d=False)
    print("%s: 2D dryice=%.2f (%.1fs)  lumped dryice=%.2f (%.1fs)"
          % (name, hd2.m_solid[-1], t2, hd0.m_solid[-1], t0))
    t = np.asarray(hd2.time_array); tl = np.asarray(hd0.time_array)

    fig, ax = plt.subplots(2, 2, figsize=(12, 8))
    kind = "riser / liquid" if riser else "no-riser / gas"
    fig.suptitle("Munkejord %s (%.0f bar, %.1f C, %.1f mm, %s) - 2D conjugate wall vs lumped"
                 % (name, P0, T0C, noz, kind), color=NAVY, fontsize=12)
    a = ax[0, 0]; a.plot(g, Pavg, color=SLATE, lw=2, label="measured")
    a.plot(t, hd2.P / 1e5, color=RED, lw=1.6, ls="--", label="model")
    a.set_ylabel("pressure [bar]"); a.set_title("Pressure"); a.legend(fontsize=8); a.grid(alpha=.3)

    a = ax[0, 1]; a.plot(g, W, color=SLATE, lw=2, label="measured")
    a.plot(t, hd2.mass_fluid, color=RED, lw=1.6, ls="--", label="model total")
    a.plot(t, hd2.m_liquid, color=NAVY, lw=1, ls="-.", label="liquid")
    a.plot(t, hd2.m_solid, color="k", lw=1, ls=":", label="dry ice")
    a.set_ylabel("mass [kg]"); a.set_title("Inventory"); a.legend(fontsize=7.5); a.grid(alpha=.3)

    a = ax[1, 0]
    a.fill_between(g, fluid.min(0), fluid.max(0), color=SLATE, alpha=0.3, label="measured band")
    a.plot(g, np.median(fluid, 0), color=SLATE, lw=1.3, label="measured midline")
    alive = np.asarray(hd2.mass_fluid) > 0.02
    cond = np.asarray(hd2.m_liquid) + np.asarray(hd2.m_solid)
    a.plot(t, np.where(alive, np.asarray(hd2.T_gas) - 273.15, np.nan), color=RED, lw=1.5, ls="--", label="model gas")
    a.plot(t, np.where(alive & (cond > 0.1), np.asarray(hd2.T_liquid) - 273.15, np.nan), color=NAVY, lw=1.5, ls="-.", label="model liq/sol")
    a.set_ylabel("fluid T [C]"); a.set_xlabel("time [s]"); a.set_title("Fluid temperature"); a.legend(fontsize=7.5); a.grid(alpha=.3)

    a = ax[1, 1]
    a.fill_between(g, iwall.min(0), iwall.max(0), color=GREEN, alpha=0.20, label="meas. inside")
    a.fill_between(g, owall.min(0), owall.max(0), color=AMBER, alpha=0.16, label="meas. outer")
    a.plot(t, hd2.T_inner_wall - 273.15, color=RED, lw=1.7, ls="-", label="2D gas in")
    a.plot(t, hd2.T_outer_wall - 273.15, color=RED, lw=1.2, ls=":", label="2D gas out")
    a.plot(t, hd2.T_inner_wall_wetted - 273.15, color=NAVY, lw=1.7, ls="-", label="2D wet in")
    a.plot(t, hd2.T_outer_wall_wetted - 273.15, color=NAVY, lw=1.2, ls=":", label="2D wet out")
    a.plot(tl, hd0.T_inner_wall_wetted - 273.15, color=DGREY, lw=1.2, ls="--", label="lumped wet")
    a.plot(tl, hd0.T_inner_wall - 273.15, color=DGREY, lw=1.0, ls="-.", label="lumped gas")
    a.set_ylabel("wall T [C]"); a.set_xlabel("time [s]"); a.set_title("Wall temperature (2D inner/outer vs lumped)")
    a.legend(fontsize=6.5, ncol=2); a.grid(alpha=.3)
    # clip to the physical range; the below-triple descent stepper (not the 2D wall) drives the
    # drained wetted node to a spurious deep-cold value once the vessel empties below the triple point.
    lo = min(np.nanmin(iwall), np.nanmin(owall)) - 12
    a.set_ylim(lo, max(np.nanmax(iwall), np.nanmax(owall)) + 6)

    for row in ax:
        for aa in row:
            aa.set_xlim(0, ET)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    for e in ("pdf", "png"):
        fig.savefig(os.path.join(OUT, "munke_%s_2dwall.%s" % (name, e)), dpi=150)
    plt.close(fig)
    print("wrote munke_%s_2dwall.pdf/png" % name)


def main():
    for name, P0, T0C, noz, riser, ET in CASES:
        try:
            figure(name, P0, T0C, noz, riser, ET)
        except Exception as ex:
            import traceback; traceback.print_exc()
            print("%s FAILED: %s" % (name, str(ex)[:150]))
    print("ALL DONE")


if __name__ == "__main__":
    main()
