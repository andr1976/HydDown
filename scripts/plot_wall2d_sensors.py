"""Temperature-transmitter reconciliation for the SINTEF tests with the 2-D conjugate wall AND
the lumped two-node base model, in manuscript (SciencePlots) style.

Reproduces the Munkejord (2026) per-height panels (Fig. 6(c),(e) for Exp53; Fig. A.9(d),(f) for
Exp72) at z = 0.77 m and 0.05 m, the two windows stacked vertically. Each panel shows the
fluid-centre temperature (TT1x4), the inner cylinder wall (TT1x2) and the outer cylinder wall
(TT1x1); measured (solid), 2-D conjugate-wall model (dashed) and lumped two-node model (dotted).
The lumped model has no through-thickness gradient, so its inner and outer walls coincide (one
line). The fluid-near-wall sensor (TT1x3) is dropped (no radial fluid gradient in the model).

Model equivalents at height z: inner/outer wall sampled from the 2-D shell (2-D) or the gas/wetted
node (lumped); the fluid is the vapour zone above the condensate bed and the liquid/solid zone
below it. No figure title (manuscript style). Output to paper/figures/wall2d_sensors_exp{53,72}.pdf.

Requires scienceplots + a LaTeX toolchain. Run from the repo root.
"""
import os, sys, zipfile, io
import numpy as np, pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import scienceplots  # noqa: F401

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "src"))
from hyddown.hdclass import HydDown
FIGDIR = os.path.join(REPO, "paper", "figures"); os.makedirs(FIGDIR, exist_ok=True)
NAVY = "#002D40"; RED = "#D61F39"; AMBER = "#E6A740"; SLATE = "#82979F"

plt.style.use(["science", "nature"])
matplotlib.rcParams.update({
    "text.usetex": True, "font.size": 10, "axes.titlesize": 11, "axes.labelsize": 11,
    "legend.fontsize": 8, "xtick.labelsize": 9, "ytick.labelsize": 9, "figure.dpi": 150,
})

HEIGHTS = [0.77, 0.05]
CH = {0.77: ("TT114", "TT112", "TT111"), 0.05: ("TT154", "TT152", "TT151")}
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


def run(P0, T0C, noz, riser, ET, wall2d, dt=0.1):
    T0 = T0C + 273.15
    vessel = {"length": 1.0, "diameter": 0.273, "thickness": 0.0254, "heat_capacity": 500,
              "density": 7950, "orientation": "vertical", "type": "Flat-end"}
    if wall2d:
        vessel.update({"wall_model": "2d", "bottom_thickness": 0.050, "flange_thickness": 0.083,
                       "lid_thickness": 0.080, "lid_diameter": 0.580})
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
            hd = HydDown(d)
            if wall2d:
                hd.wall2d_report_heights = list(HEIGHTS)
            hd.run(disable_pbar=True); return hd
        except Exception:
            if step == 0.02:
                raise
    return hd


def _bed_height(hd):
    vl, vs = hd.release_model.v_l, hd.release_model.v_s
    bedV = np.asarray(hd.m_liquid) * vl + np.asarray(hd.m_solid) * vs
    bh = np.array([hd.inner_vol.h_from_V(v) if v > 1e-9 else 0.0 for v in bedV])
    return np.nan_to_num(bh, nan=0.0)


def fluid_at(hd, z):
    Tg = np.asarray(hd.T_gas) - 273.15
    Tl = np.asarray(hd.T_liquid) - 273.15
    fld = np.where(z <= _bed_height(hd), Tl, Tg)
    return np.where(np.asarray(hd.mass_fluid) > 0.02, fld, np.nan)


def lumped_wall_at(hd, z):
    Tg = np.asarray(hd.T_vessel) - 273.15
    Tw = np.asarray(hd.T_vessel_wetted) - 273.15
    return np.where(z <= _bed_height(hd), Tw, Tg)


def figure(name, P0, T0C, noz, riser, ET, letters):
    zf = zipfile.ZipFile(os.path.join(REPO, "background", "19589510.zip"))
    g, meas = measured(zf, name, ET)
    hd2 = run(P0, T0C, noz, riser, ET, wall2d=True)
    hd0 = run(P0, T0C, noz, riser, ET, wall2d=False)
    t2 = np.asarray(hd2.time_array); t0 = np.asarray(hd0.time_array)

    fig, axes = plt.subplots(2, 1, figsize=(3.5, 5.4), sharex=True)
    for k, (h, ax, letter) in enumerate(zip(HEIGHTS, axes, letters)):
        fch, iwch, owch = CH[h]
        # measured (solid)
        ax.plot(g, meas[fch], color=RED, lw=1.3)
        ax.plot(g, meas[iwch], color=AMBER, lw=1.3)
        ax.plot(g, meas[owch], color=NAVY, lw=1.3)
        # 2-D model (dashed)
        ax.plot(t2, fluid_at(hd2, h), color=RED, lw=1.2, ls="--")
        ax.plot(t2, hd2.T_wall_h_in[:, k] - 273.15, color=AMBER, lw=1.2, ls="--")
        ax.plot(t2, hd2.T_wall_h_out[:, k] - 273.15, color=NAVY, lw=1.2, ls="--")
        # lumped two-node model (dotted); inner=outer -> single wall line
        ax.plot(t0, fluid_at(hd0, h), color=RED, lw=1.1, ls=":")
        ax.plot(t0, lumped_wall_at(hd0, h), color=AMBER, lw=1.1, ls=":")
        # y-limits from measured + 2-D only (exclude lumped extremes)
        lo = np.nanmin([np.nanmin(meas[fch]), np.nanmin(fluid_at(hd2, h)),
                        np.nanmin(hd2.T_wall_h_in[:, k] - 273.15)]) - 6
        hi = np.nanmax([np.nanmax(meas[owch]), np.nanmax(fluid_at(hd2, h))]) + 4
        ax.set_ylim(lo, hi)
        ax.set_ylabel(r"Temperature [$^\circ$C]")
        ax.text(0.03, 0.05, r"(%s) $z = %.2f$~m" % (letter, h), transform=ax.transAxes,
                fontsize=9, va="bottom", ha="left")
        ax.set_xlim(0, ET)
    axes[-1].set_xlabel(r"Time [s]")

    # two-key legend on the top panel
    col = [Line2D([], [], color=RED, lw=1.4), Line2D([], [], color=AMBER, lw=1.4),
           Line2D([], [], color=NAVY, lw=1.4)]
    sty = [Line2D([], [], color="k", lw=1.4, ls="-"), Line2D([], [], color="k", lw=1.4, ls="--"),
           Line2D([], [], color="k", lw=1.4, ls=":")]
    # Single two-row legend below the panels (colours on top, line styles below), clear of the
    # curves. Handles interleaved so the column-major fill gives a colour row and a style row.
    handles = [col[0], sty[0], col[1], sty[1], col[2], sty[2]]
    labels = ["Fluid centre", "Measured", "Inner wall", "2-D wall", "Outer wall", "Lumped 2-node"]
    fig.tight_layout(rect=[0, 0.085, 1, 1])
    fig.subplots_adjust(hspace=0.08)
    fig.legend(handles, labels, loc="lower center", bbox_to_anchor=(0.5, 0.0), ncol=3,
               frameon=True, handlelength=1.8, columnspacing=1.3, fontsize=8)
    for e in ("pdf", "png"):
        try:
            fig.savefig(os.path.join(FIGDIR, "wall2d_sensors_%s.%s" % (name.lower(), e)),
                        dpi=200, bbox_inches="tight")
        except OSError as ex:
            print("  skip %s (%s)" % (e, str(ex)[:60]))
    plt.close(fig)
    print("wrote paper/figures/wall2d_sensors_%s.pdf/png" % name.lower())


def main():
    for name, P0, T0C, noz, riser, ET, letters in CASES:
        try:
            figure(name, P0, T0C, noz, riser, ET, letters)
        except Exception as ex:
            import traceback; traceback.print_exc(); print("%s FAILED: %s" % (name, str(ex)[:150]))
    print("ALL DONE")


if __name__ == "__main__":
    main()
