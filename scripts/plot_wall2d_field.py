"""Wall cross-section temperature field (2-D conjugate wall) at several times, for both SINTEF
experiments stacked in one manuscript figure (Exp53 on top, Exp72 below).

Each row runs a case with the 2-D wall, captures the r-z steel-temperature field at a set of times,
and draws mirrored axisymmetric cross-sections (pcolormesh) sharing a per-row colour bar. The
cavity and exterior are blank; the condensate level is marked. SciencePlots style, no figure title.

Writes paper/figures/wall2d_field.pdf (+ .png).  Run from the repo root.
"""
import os, sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import scienceplots  # noqa: F401

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "src"))
from hyddown.hdclass import HydDown
FIGDIR = os.path.join(REPO, "paper", "figures"); os.makedirs(FIGDIR, exist_ok=True)
NAVY = "#002D40"

plt.style.use(["science", "nature"])
matplotlib.rcParams.update({
    "text.usetex": True, "font.size": 9, "axes.labelsize": 10, "figure.dpi": 150,
})

# name, P0_bar, T0_C, nozzle_mm, riser, end_time, snapshot times [s], row label
CASES = [
    ("Exp53", 119.5, 24.4, 8.0, True, 100.0, [5, 15, 30, 45, 80], "Exp53 (riser / liquid)"),
    ("Exp72", 119.0, 24.9, 6.5, False, 400.0, [20, 80, 155, 250, 350], "Exp72 (no-riser / gas)"),
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
    fld = np.full((W.nr, W.nz), np.nan)
    for n, (ir, iz) in enumerate(W.steel_cells):
        fld[ir, iz] = Tvec[n]
    return np.ma.masked_invalid(fld)


def make_figure(mode, outname):
    """mode='true': physical coordinates, equal aspect. mode='grid': uniform cell spacing (each
    control volume drawn the same size, i.e. cell-index coordinates), which enlarges the thin
    shell and bottom-plate cells so their gradients are visible; the axis labels then map to the
    (non-uniform) physical positions of the cell edges."""
    ncol = max(len(c[6]) for c in CASES)
    figsize = (9.6, 6.2) if mode == "grid" else (7.4, 6.6)
    fig, axes = plt.subplots(len(CASES), ncol, figsize=figsize)
    cmap = plt.get_cmap("RdYlBu_r")
    for row, (name, P0, T0C, noz, riser, ET, snaps, rowlab) in enumerate(CASES):
        hd = run(P0, T0C, noz, riser, ET, snaps)
        W = hd._wall2d
        snapshots = sorted(hd.wall2d_snapshots, key=lambda r: r["t"])
        re, ze, g = W.re, W.ze, W.geom
        if mode == "grid":
            ir_edge, iz_edge = np.arange(len(re), dtype=float), np.arange(len(ze), dtype=float)
            re_d, ze_d = ir_edge, iz_edge
            r_in_d = np.interp(g["r_in"], re, ir_edge); r_lid_d = np.interp(g["r_lid"], re, ir_edge)
            z_flt = g["z_top_cavity"]; z_top = z_flt + g["t_lid"]   # flange top, lid top
            ztk_phys = [-g["t_bot"], 0.0, 0.25, 0.5, 1.0, z_flt, z_top]
            ztk_disp = np.interp(ztk_phys, ze, iz_edge)
        else:
            re_d, ze_d = re, ze
            r_in_d, r_lid_d = g["r_in"], g["r_lid"]
        r_full = np.concatenate([-re_d[::-1], re_d[1:]])
        Tmin = min(np.nanmin(field_grid(W, s["T"])) for s in snapshots) - 273.15
        Tmax = max(np.nanmax(field_grid(W, s["T"])) for s in snapshots) - 273.15
        norm = Normalize(vmin=np.floor(Tmin / 5) * 5, vmax=np.ceil(Tmax / 5) * 5)
        for j in range(ncol):
            ax = axes[row, j]
            if j < len(snapshots):
                s = snapshots[j]
                fld = field_grid(W, s["T"]) - 273.15
                fld_full = np.ma.concatenate([fld[::-1, :], fld], axis=0)
                ax.pcolormesh(r_full, ze_d, fld_full.T, cmap=cmap, norm=norm, shading="flat")
                lvl = s["level"]
                if lvl and lvl > 1e-3:
                    zl = np.interp(lvl, ze, ze_d) if mode == "grid" else lvl
                    ax.plot([-r_in_d, r_in_d], [zl, zl], color="k", lw=0.9, ls="--")
                ax.set_title(r"$t=%d$\,s, $P=%.1f$\,bar" % (s["t"], s["P"] / 1e5), fontsize=7.5)
                ax.set_aspect("auto" if mode == "grid" else "equal")
                ax.set_xlim(-r_lid_d * 1.03, r_lid_d * 1.03)
                if mode == "grid" and row == len(CASES) - 1:   # radial ticks on the Exp72 row only
                    ax.set_xticks([-r_lid_d, -r_in_d, 0.0, r_in_d, r_lid_d])
                    ax.set_xticklabels(["0.29", "0.14", "0", "0.14", "0.29"], fontsize=6)
                    ax.set_xlabel(r"$r$ [m]", fontsize=8)
                else:
                    ax.set_xticks([])
                if j == 0:
                    ax.set_ylabel(r"$z$ [m]")
                    if mode == "grid":
                        ax.set_yticks(ztk_disp); ax.set_yticklabels(["%.2f" % z for z in ztk_phys])
                    ax.tick_params(labelsize=7)
                else:
                    ax.set_yticks([])
            else:
                ax.axis("off")
        sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
        cb = fig.colorbar(sm, ax=axes[row, :].tolist(), fraction=0.018, pad=0.01)
        cb.set_label(r"wall $T$ [$^\circ$C]", fontsize=8)
        cb.ax.tick_params(labelsize=7)
        axes[row, 0].annotate(rowlab, xy=(-0.62, 0.5), xycoords="axes fraction",
                              rotation=90, ha="center", va="center", fontsize=9, color=NAVY)
    for e in ("pdf", "png"):
        fig.savefig(os.path.join(FIGDIR, "%s.%s" % (outname, e)), dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("wrote paper/figures/%s.pdf/png" % outname)


def main():
    make_figure(mode="true", outname="wall2d_field")
    make_figure(mode="grid", outname="wall2d_field_grid")


if __name__ == "__main__":
    main()
