"""Wall cross-section temperature field (2-D conjugate wall) at several times.

Runs a SINTEF case with the 2-D wall, captures the r-z steel-temperature field at a set of
times, and draws a row of axisymmetric cross-sections (mirrored about the axis) as pcolormesh
panels with a shared colour bar. The cavity and exterior are left blank; the condensate level is
marked. Writes validation/munkejord/wall2d_field_<name>.pdf (+ .png).

Run from the repo root:  python scripts/plot_wall2d_field.py
"""
import os, sys
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "src"))
from hyddown.hdclass import HydDown
OUT = os.path.join(REPO, "validation", "munkejord"); os.makedirs(OUT, exist_ok=True)
NAVY = "#002D40"

# name, P0_bar, T0_C, nozzle_mm, riser, end_time, snapshot times [s]
CASES = [
    ("Exp72", 119.0, 24.9, 6.5, False, 400.0, [20, 80, 155, 250, 350]),
    ("Exp53", 119.5, 24.4, 8.0, True, 100.0, [5, 15, 30, 45, 80]),
]


def run(P0, T0C, noz, riser, ET, snaps, dt=0.1):
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
            hd = HydDown(d); hd.wall2d_snap_times = list(snaps); hd.run(disable_pbar=True); return hd
        except Exception:
            if step == 0.02:
                raise
    return hd


def field_grid(W, Tvec):
    """Map the unknown vector to a (nr, nz) array (NaN off-steel), masked for plotting."""
    fld = np.full((W.nr, W.nz), np.nan)
    for n, (ir, iz) in enumerate(W.steel_cells):
        fld[ir, iz] = Tvec[n]
    return np.ma.masked_invalid(fld)


def plot_case(name, P0, T0C, noz, riser, ET, snaps):
    hd = run(P0, T0C, noz, riser, ET, snaps)
    W = hd._wall2d
    snapshots = sorted(hd.wall2d_snapshots, key=lambda r: r["t"])
    if not snapshots:
        print("%s: no snapshots captured" % name); return
    # mirrored radial edges [-r_max .. r_max] and field
    re, ze = W.re, W.ze
    r_full = np.concatenate([-re[::-1], re[1:]])
    Tmin = min(np.nanmin(field_grid(W, s["T"])) for s in snapshots) - 273.15
    Tmax = max(np.nanmax(field_grid(W, s["T"])) for s in snapshots) - 273.15
    norm = Normalize(vmin=np.floor(Tmin / 5) * 5, vmax=np.ceil(Tmax / 5) * 5)
    cmap = plt.get_cmap("RdYlBu_r")

    n = len(snapshots)
    fig, axes = plt.subplots(1, n, figsize=(2.4 * n + 1.2, 6.6), sharey=True)
    if n == 1:
        axes = [axes]
    for ax, s in zip(axes, snapshots):
        fld = field_grid(W, s["T"]) - 273.15
        fld_full = np.ma.concatenate([fld[::-1, :], fld], axis=0)   # mirror about axis
        pc = ax.pcolormesh(r_full, ze, fld_full.T, cmap=cmap, norm=norm, shading="flat")
        # condensate level (gas/condensate boundary) across the bore
        lvl = s["level"]
        if lvl and lvl > 1e-3:
            ax.plot([-W.geom["r_in"], W.geom["r_in"]], [lvl, lvl], color="k", lw=1.1, ls="--")
        ax.set_aspect("equal")
        ax.set_xlim(-W.geom["r_lid"] * 1.03, W.geom["r_lid"] * 1.03)
        ax.set_title("t = %.0f s\nP = %.1f bar" % (s["t"], s["P"] / 1e5), fontsize=9, color=NAVY)
        ax.set_xticks([])
        ax.tick_params(labelsize=8)
    axes[0].set_ylabel("axial position z [m]", fontsize=9)
    kind = "riser / liquid" if riser else "no-riser / gas"
    fig.suptitle("%s (%.0f bar, %.1f C, %.1f mm, %s) - wall temperature field [C]"
                 % (name, P0, T0C, noz, kind), color=NAVY, fontsize=12, y=0.99)
    cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), ax=axes, fraction=0.03, pad=0.02)
    cb.set_label("wall temperature [C]", fontsize=9)
    for e in ("pdf", "png"):
        fig.savefig(os.path.join(OUT, "wall2d_field_%s.%s" % (name, e)), dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("%s: %d snapshots -> wall2d_field_%s.pdf/png  (T range %.0f..%.0f C)"
          % (name, len(snapshots), name, Tmin, Tmax))


def main():
    for name, P0, T0C, noz, riser, ET, snaps in CASES:
        try:
            plot_case(name, P0, T0C, noz, riser, ET, snaps)
        except Exception as ex:
            import traceback; traceback.print_exc(); print("%s FAILED: %s" % (name, str(ex)[:150]))
    print("ALL DONE")


if __name__ == "__main__":
    main()
