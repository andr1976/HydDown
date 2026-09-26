"""Graphical abstract: CO2 phase diagram with the blowdown path (three-phase schematic),
the pressure-time 'fingerprint' of a gas-release blowdown (flash -> triple-point plateau ->
gas/solid tail), and the 2-D conjugate-wall temperature field. One coherent case (ECCSEL/SINTEF
Exp72, no-riser gas release). Writes paper/figures/graphical_abstract.pdf (+ .png).

Run from the repo root:  python scripts/plot_graphical_abstract.py
"""
import os, sys
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.colors import Normalize
import CoolProp.CoolProp as CP

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "src"))
from hyddown.hdclass import HydDown
FIGDIR = os.path.join(REPO, "paper", "figures"); os.makedirs(FIGDIR, exist_ok=True)
NAVY = "#002D40"; RED = "#D61F39"; AMBER = "#E6A740"; SLATE = "#82979F"
plt.rcParams.update({"font.family": "sans-serif", "font.size": 8.5})

TCRIT, PCRIT = 304.128, 73.77          # CO2 critical (K, bar)
TTRIP, PTRIP = 216.592, 5.185          # CO2 triple point (K, bar)


def run_exp72(ET=400.0):
    T0 = 24.9 + 273.15
    vessel = {"length": 1.0, "diameter": 0.273, "thickness": 0.0254, "heat_capacity": 500,
              "density": 7950, "orientation": "vertical", "type": "Flat-end", "wall_model": "2d",
              "bottom_thickness": 0.050, "flange_thickness": 0.083, "lid_thickness": 0.080,
              "lid_diameter": 0.580}
    d = {"vessel": vessel, "initial": {"temperature": T0, "pressure": 119.0e5, "fluid": "CO2"},
         "calculation": {"type": "energybalance", "time_step": 0.1, "end_time": ET,
                         "non_equilibrium": True, "h_gas_liquid": "calc_two_sided"},
         "valve": {"flow": "discharge", "type": "none", "back_pressure": 101325.},
         "release": {"type": "gas", "diameter": 0.0065, "discharge_coef": 1.0,
                     "liquid_nonequilibrium": 0.0, "back_pressure": 101325., "atm_pressure": 101325.,
                     "eos": "CoolProp", "solid_in_vessel": True, "solid_h_inner": "cooper",
                     "solid_h_gas_wall": "churchill", "solid_h_gas_liquid": 0.0, "solid_h_gas_solid": 3.0,
                     "discharge_location": 0.0},
         "heat_transfer": {"type": "specified_h", "temp_ambient": T0, "h_outer": 0.0, "h_inner": "churchill"}}
    hd = HydDown(d); hd.wall2d_snap_times = [20, 250, 350]; hd.run(disable_pbar=True); return hd


def field_grid(W, Tvec):
    fld = np.full((W.nr, W.nz), np.nan)
    for n, (ir, iz) in enumerate(W.steel_cells):
        fld[ir, iz] = Tvec[n]
    return np.ma.masked_invalid(fld)


def phase_diagram(ax, hd):
    # vaporisation line (liquid-vapour, Q=0) from triple to critical
    Tv = np.linspace(TTRIP, TCRIT, 120)
    Pv = np.array([CP.PropsSI("P", "T", t, "Q", 0, "CO2") for t in Tv]) / 1e5
    ax.plot(Tv - 273.15, Pv, color=NAVY, lw=1.6)
    # sublimation line (solid-vapour) from the model, below the triple point
    rm = hd.release_model
    Ps = np.linspace(1.0, PTRIP, 60)
    Ts = np.array([rm._T_sub_of_P(p * 1e5) for p in Ps])
    ax.plot(Ts - 273.15, Ps, color=NAVY, lw=1.6)
    # melting line (solid-liquid), steep positive slope (stylised)
    ax.plot([TTRIP - 273.15, TTRIP - 273.15 + 2.2], [PTRIP, 210], color=NAVY, lw=1.6)
    # points
    ax.plot(TTRIP - 273.15, PTRIP, "o", color=NAVY, ms=4)
    ax.plot(TCRIT - 273.15, PCRIT, "o", color=NAVY, ms=4)
    ax.annotate("triple point", (TTRIP - 273.15, PTRIP), (-95, 8.5), fontsize=6.8, color=NAVY)
    ax.annotate("critical\npoint", (TCRIT - 273.15, PCRIT), (34, 72), fontsize=6.8,
                color=NAVY, ha="left", va="center")
    # region labels
    ax.text(-85, 40, "solid", color=SLATE, fontsize=8, style="italic")
    ax.text(-5, 90, "liquid", color=SLATE, fontsize=8, style="italic")
    ax.text(20, 3, "vapour", color=SLATE, fontsize=8, style="italic")
    ax.text(38, 100, "s.c.", color=SLATE, fontsize=7, style="italic")
    # blowdown path: dense liquid -> flash onto sat line -> down sat line -> sublimation tail
    Tsat25 = CP.PropsSI("P", "T", 25 + 273.15, "Q", 0, "CO2") / 1e5
    path_T = [25]; path_P = [119]
    path_T += [25]; path_P += [Tsat25]              # vertical flash onto the saturation line
    msk = (Tv - 273.15) <= 25
    path_T += list((Tv[msk] - 273.15)[::-1]); path_P += list(Pv[msk][::-1])   # down the sat line
    path_T += list(Ts - 273.15)[::-1]; path_P += list(Ps)[::-1]               # sublimation tail
    ax.plot(path_T, path_P, color=RED, lw=2.4, zorder=5, solid_capstyle="round")
    ax.annotate("", xy=(path_T[-1], path_P[-1]), xytext=(path_T[-6], path_P[-6]),
                arrowprops=dict(arrowstyle="-|>", color=RED, lw=2.4))
    ax.text(2, 118, "blowdown", color=RED, fontsize=7.5, fontweight="bold")
    ax.set_yscale("log"); ax.set_ylim(0.9, 205); ax.set_xlim(-100, 55)
    ax.set_xlabel("Temperature [$^\\circ$C]", fontsize=8); ax.set_ylabel("Pressure [bar]", fontsize=8)
    ax.set_title("CO$_2$ phase diagram", color=NAVY, fontsize=9)
    ax.tick_params(labelsize=7)


def fingerprint(ax, hd):
    t = np.asarray(hd.time_array); P = np.asarray(hd.P) / 1e5
    ms = np.asarray(hd.m_solid)
    ax.plot(t, P, color=RED, lw=2.2, zorder=4)
    ax.set_ylabel("Pressure [bar]", color=RED, fontsize=8); ax.set_yscale("log")
    ax.set_ylim(0.8, 140); ax.set_xlim(0, 400); ax.tick_params(axis="y", labelcolor=RED, labelsize=7)
    ax.tick_params(axis="x", labelsize=7); ax.set_xlabel("Time [s]", fontsize=8)
    # dry-ice mass on a secondary axis
    ax2 = ax.twinx()
    ax2.fill_between(t, 0, ms, color=SLATE, alpha=0.25, zorder=1)
    ax2.plot(t, ms, color=SLATE, lw=1.3, zorder=2)
    ax2.set_ylabel("Retained dry ice [kg]", color=SLATE, fontsize=8)
    ax2.set_ylim(0, max(ms) * 2.0); ax2.tick_params(axis="y", labelcolor=SLATE, labelsize=7)
    # triple-point line + stage annotations
    ax.axhline(PTRIP, color=NAVY, lw=0.8, ls=":")
    ax.text(370, PTRIP * 1.15, "triple point", color=NAVY, fontsize=6.3, ha="right")
    ax.annotate("1. depressurisation", (85, 48), fontsize=6.8, color=NAVY, ha="center")
    ax.annotate("2. triple-point\nplateau", (200, 12), fontsize=6.8, color=NAVY, ha="center")
    ax.annotate("3. gas +\nsolid tail", xy=(338, 1.6), xytext=(356, 24), fontsize=6.8,
                color=NAVY, ha="center", va="center",
                arrowprops=dict(arrowstyle="-", color=NAVY, lw=0.6))
    ax.set_title("Gas-release blowdown fingerprint", color=NAVY, fontsize=9)


def wall_field(fig, gs_cell, hd):
    W = hd._wall2d
    snaps = sorted(hd.wall2d_snapshots, key=lambda r: r["t"])
    re, ze, g = W.re, W.ze, W.geom
    r_full = np.concatenate([-re[::-1], re[1:]])
    Tmin = min(np.nanmin(field_grid(W, s["T"])) for s in snaps) - 273.15
    Tmax = max(np.nanmax(field_grid(W, s["T"])) for s in snaps) - 273.15
    norm = Normalize(vmin=np.floor(Tmin / 5) * 5, vmax=np.ceil(Tmax / 5) * 5); cmap = plt.get_cmap("RdYlBu_r")
    inner = gridspec.GridSpecFromSubplotSpec(1, len(snaps) + 1, subplot_spec=gs_cell,
                                             width_ratios=[1] * len(snaps) + [0.13], wspace=0.06)
    pc = None
    for j, s in enumerate(snaps):
        ax = fig.add_subplot(inner[0, j])
        fld = field_grid(W, s["T"]) - 273.15
        fld_full = np.ma.concatenate([fld[::-1, :], fld], axis=0)
        pc = ax.pcolormesh(r_full, ze, fld_full.T, cmap=cmap, norm=norm, shading="flat")
        lvl = s["level"]
        if lvl and lvl > 1e-3:
            ax.plot([-g["r_in"], g["r_in"]], [lvl, lvl], color="k", lw=0.7, ls="--")
        ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
        ax.set_xlim(-g["r_lid"] * 1.03, g["r_lid"] * 1.03)
        ax.set_title("%ds" % s["t"], fontsize=7, color=NAVY, pad=2)
        if j == 0:
            ax.set_title("%ds" % s["t"], fontsize=7, color=NAVY, pad=2)
    cax = fig.add_subplot(inner[0, -1])
    cb = fig.colorbar(pc, cax=cax); cb.ax.tick_params(labelsize=6); cb.set_label("wall $T$ [$^\\circ$C]", fontsize=7)
    # a title over the strip
    fig.text(gs_cell.get_position(fig).x0 + 0.5 * gs_cell.get_position(fig).width,
             gs_cell.get_position(fig).y1 + 0.01, "2-D wall temperature field",
             ha="center", va="bottom", color=NAVY, fontsize=9)


def main():
    hd = run_exp72()
    fig = plt.figure(figsize=(10.4, 3.6))
    gs = gridspec.GridSpec(1, 3, width_ratios=[1.0, 1.2, 1.42], wspace=0.40,
                           left=0.055, right=0.975, top=0.85, bottom=0.155)
    phase_diagram(fig.add_subplot(gs[0, 0]), hd)
    fingerprint(fig.add_subplot(gs[0, 1]), hd)
    wall_field(fig, gs[0, 2], hd)
    for e in ("pdf", "png"):
        fig.savefig(os.path.join(FIGDIR, "graphical_abstract.%s" % e), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("wrote paper/figures/graphical_abstract.pdf/png")


if __name__ == "__main__":
    main()
