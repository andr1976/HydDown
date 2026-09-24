"""Reproduce the Munkejord (2026) per-height temperature panels with the 2-D wall model.

Fig. 6(c),(e) are Test 53 at z = 0.77 m and 0.05 m; Fig. A.9(d),(f) are Test 72 at the same two
heights. Each panel shows the fluid-centre temperature (TT1x4, red), the inner cylinder wall
(TT1x2, orange) and the outer cylinder wall (TT1x1, blue); measured = solid, model = dashed. The
fluid-near-wall sensor (TT1x3, green) is dropped, as the model has no radial fluid gradient
(matching the paper, which also omits it).

Model equivalents at height z: inner/outer wall are sampled from the 2-D shell at z; the fluid is
the gas zone when z is above the condensate bed and the liquid/solid zone when below it.

Run from the repo root:  python scripts/plot_wall2d_sensors.py
"""
import os, sys, zipfile, io
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "src"))
from hyddown.hdclass import HydDown
OUT = os.path.join(REPO, "validation", "munkejord"); os.makedirs(OUT, exist_ok=True)
NAVY = "#002D40"; RED = "#D61F39"; AMBER = "#E6A740"; SLATE = "#82979F"

HEIGHTS = [0.77, 0.05]                      # panel heights [m]  (Fig c/d = 0.77, e/f = 0.05)
# measured channels at each height: (fluid centre TT1x4, inner wall TT1x2, outer wall TT1x1)
CH = {0.77: ("TT114", "TT112", "TT111"), 0.05: ("TT154", "TT152", "TT151")}
# name, P0_bar, T0_C, nozzle_mm, riser, end_time, panel-letters
CASES = [("Exp53", 119.5, 24.4, 8.0, True, 100.0, ("c", "e")),
         ("Exp72", 119.0, 24.9, 6.5, False, 400.0, ("d", "f"))]


def measured(z, name, ET):
    xlsx = [n for n in z.namelist() if n.split("/")[-1].startswith(name + "_")][0]
    D = pd.read_excel(io.BytesIO(z.read(xlsx)), sheet_name="Data", header=None, engine="calamine")
    names = [str(x) for x in D.iloc[0].tolist()]; col = {n: i for i, n in enumerate(names)}
    dat = D.iloc[3:].reset_index(drop=True)

    def raw(tn, vn):
        tt = pd.to_numeric(dat[col[tn]], errors="coerce").values
        vv = pd.to_numeric(dat[col[vn]], errors="coerce").values
        m = np.isfinite(tt) & np.isfinite(vv); return tt[m], vv[m]
    tp, p = raw("tPT162", "PT162"); p0 = np.nanmax(p[:200])
    below = np.where(p < p0 - 5.0)[0]; t0 = tp[below[0]] if len(below) else 0.0
    g = np.arange(0.0, min(tp.max(), ET + t0) - t0, 0.5)
    out = {}
    for h in HEIGHTS:
        for ch in CH[h]:
            tt, vv = raw("tTT", ch); out[ch] = np.interp(g, tt - t0, vv)
    return g, out


def run(P0, T0C, noz, riser, ET, dt=0.1):
    T0 = T0C + 273.15
    vessel = {"length": 1.0, "diameter": 0.273, "thickness": 0.0254, "heat_capacity": 500,
              "density": 7950, "orientation": "vertical", "type": "Flat-end", "wall_model": "2d",
              "bottom_thickness": 0.050, "flange_thickness": 0.083, "lid_thickness": 0.080,
              "lid_diameter": 0.580}
    d = {"vessel": vessel, "initial": {"temperature": T0, "pressure": P0 * 1e5, "fluid": "CO2"},
         "calculation": {"type": "energybalance", "time_step": dt, "end_time": ET,
                         "non_equilibrium": True, "h_gas_liquid": "calc_two_sided"},
         "valve": {"flow": "discharge", "type": "none", "back_pressure": 101325.},
         "release": {"type": "liquid" if riser else "gas", "diameter": noz / 1000., "discharge_coef": 1.0,
                     "liquid_nonequilibrium": 0.0, "back_pressure": 101325., "atm_pressure": 101325.,
                     "eos": "CoolProp", "solid_in_vessel": True, "solid_h_inner": "cooper",
                     "solid_h_gas_wall": "churchill", "solid_h_gas_liquid": 0.0, "solid_h_gas_solid": 3.0,
                     "discharge_location": 0.009 if riser else 0.0},
         "heat_transfer": {"type": "specified_h", "temp_ambient": T0, "h_outer": 0.0, "h_inner": "churchill"}}
    for step in (dt, 0.02):
        d["calculation"]["time_step"] = step
        try:
            hd = HydDown(d); hd.wall2d_report_heights = list(HEIGHTS); hd.run(disable_pbar=True); return hd
        except Exception:
            if step == 0.02:
                raise
    return hd


def model_fluid_at(hd, z):
    Tg = np.asarray(hd.T_gas) - 273.15
    Tl = np.asarray(hd.T_liquid) - 273.15
    vl, vs = hd.release_model.v_l, hd.release_model.v_s
    bedV = np.asarray(hd.m_liquid) * vl + np.asarray(hd.m_solid) * vs
    bedh = np.array([hd.inner_vol.h_from_V(v) if v > 1e-9 else 0.0 for v in bedV])
    bedh = np.nan_to_num(bedh, nan=0.0)
    fld = np.where(z <= bedh, Tl, Tg)
    alive = np.asarray(hd.mass_fluid) > 0.02
    return np.where(alive, fld, np.nan)


def figure(name, P0, T0C, noz, riser, ET, letters):
    zf = zipfile.ZipFile(os.path.join(REPO, "background", "19589510.zip"))
    g, meas = measured(zf, name, ET)
    hd = run(P0, T0C, noz, riser, ET)
    t = np.asarray(hd.time_array)
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.6), sharey=True)
    for k, (h, ax, letter) in enumerate(zip(HEIGHTS, axes, letters)):
        fch, iwch, owch = CH[h]
        # measured (solid)
        ax.plot(g, meas[fch], color=RED, lw=1.8, label="fluid centre (meas.)")
        ax.plot(g, meas[iwch], color=AMBER, lw=1.8, label="inner wall (meas.)")
        ax.plot(g, meas[owch], color=NAVY, lw=1.8, label="outer wall (meas.)")
        # model (dashed)
        ax.plot(t, model_fluid_at(hd, h), color=RED, lw=1.6, ls="--", label="fluid (model)")
        ax.plot(t, hd.T_wall_h_in[:, k] - 273.15, color=AMBER, lw=1.6, ls="--", label="inner wall (model)")
        ax.plot(t, hd.T_wall_h_out[:, k] - 273.15, color=NAVY, lw=1.6, ls="--", label="outer wall (model)")
        ax.set_title("(%s)  z = %.2f m" % (letter, h), fontsize=10, color=NAVY)
        ax.set_xlabel("time [s]"); ax.set_xlim(0, ET); ax.grid(alpha=.3)
    axes[0].set_ylabel("temperature [C]")
    axes[0].legend(fontsize=7.5, ncol=2, loc="lower right")
    kind = "riser / liquid" if riser else "no-riser / gas"
    fig.suptitle("%s (%.0f bar, %.1f C, %.1f mm, %s) - fluid & wall temperatures at sensor heights"
                 % (name, P0, T0C, noz, kind), color=NAVY, fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    for e in ("pdf", "png"):
        fig.savefig(os.path.join(OUT, "wall2d_sensors_%s.%s" % (name, e)), dpi=150)
    plt.close(fig)
    print("wrote wall2d_sensors_%s.pdf/png" % name)


def main():
    for name, P0, T0C, noz, riser, ET, letters in CASES:
        try:
            figure(name, P0, T0C, noz, riser, ET, letters)
        except Exception as ex:
            import traceback; traceback.print_exc(); print("%s FAILED: %s" % (name, str(ex)[:150]))
    print("ALL DONE")


if __name__ == "__main__":
    main()
