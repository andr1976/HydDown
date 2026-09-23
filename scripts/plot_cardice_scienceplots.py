"""Re-plot selected CARDICE tests (7 and 8, the two used in the main paper) with the
SciencePlots 'science' + 'nature' styles and LaTeX text rendering, so that CO2 is set
as CO$_2$. Axis labels are capitalised and enlarged and the legend text is enlarged.

The data loading, the HydDown run and the measured-trace time alignment are reused
unchanged from ``validate_cardice_reconciled.py`` (so the physics and the alignment are
identical to the standard validation figures); only the presentation differs. Output is
written straight to the paper figure directory: paper/figures/cardice_t{7,8}.pdf (+ .png).

Requires: scienceplots, and a LaTeX toolchain (latex + dvipng) for text.usetex.

Run from the repo root:  python scripts/plot_cardice_scienceplots.py
"""
import os, sys
import numpy as np, yaml
import matplotlib
import matplotlib.pyplot as plt
import scienceplots  # noqa: F401  registers the 'science' and 'nature' styles

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import validate_cardice_reconciled as vcr  # data / model / alignment helpers
from hyddown.hdclass import HydDown

NAVY, RED, AMBER, SLATE, GREEN = vcr.NAVY, vcr.RED, vcr.AMBER, vcr.SLATE, vcr.GREEN
OUT = vcr.OUT
FIGDIR = vcr.REPO + "paper/figures/"
TESTS = dict(vcr.TESTMAP)     # all six CARDICE pure-CO2 tests (5-10); the SI shows all of them
XMAX = {8: 1.5}               # per-test x-axis limit [h]; default is the model end time
WALL_LEGEND_LOC = {7: "lower left"}  # per-test wall-panel legend location; default "best"

# SciencePlots 'science' enables text.usetex; 'nature' overlays the Nature layout.
# The rc overrides below only enlarge the axis-label and legend fonts and set a
# paper-appropriate figure size; the LaTeX preamble from the styles is left intact so
# CO$_2$ and $^\circ$C render through LaTeX.
plt.style.use(["science", "nature"])
matplotlib.rcParams.update({
    "text.usetex": True,
    "font.size": 10,
    "axes.titlesize": 12,
    "axes.labelsize": 12,   # axis titles: moderately enlarged
    "legend.fontsize": 9,   # legend text: moderately enlarged
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "figure.dpi": 150,
})


def plot_test(t, pfx):
    P = vcr.rd(pfx, "pressures"); M = vcr.rd(pfx, "mass-flowrate")
    IT = vcr.rd(pfx, "internal-temperatures"); WT = vcr.rd(pfx, "inside-wall-temperatures")
    okp = np.isfinite(P.iloc[:, 0].values) & np.isfinite(P.iloc[:, 1].values)
    tP, ps = vcr.merge_across_gaps(P.iloc[:, 0].values[okp], P.iloc[:, 1].values[okp])

    d = yaml.safe_load(open(OUT + "CARDICE_test%d.yml" % t))
    hd = HydDown(d); hd.run(disable_pbar=True)
    tm_s = np.asarray(hd.time_array); Pm_bar = np.asarray(hd.P) / 1e5
    tm = tm_s / 3600.0; ET = tm.max()

    # Time-align the measured trace to the model (gas: mid-blowdown pressure reference;
    # liquid: inventory-decline onset), identical to validate_cardice_reconciled.py.
    tM, mm = vcr.merge_across_gaps(M.iloc[:, 0].values, M.iloc[:, 1].values)
    Mm = np.asarray(hd.mass_fluid)
    if d["release"]["type"] == "gas":
        p0 = np.nanmax(ps[:120]); Pref = p0 - max(1.0, 0.15 * (p0 - 1.0))

        def _cross(tt, pp):
            idx = np.where(pp < Pref)[0]; return tt[idx[0]] if len(idx) else tt[-1]
        shift = _cross(tP, ps) - _cross(tm_s, Pm_bar)
    else:
        def _onset(tt, yy):
            y0 = np.nanmax(yy[:120]); thr = max(3.0, 0.005 * y0)
            idx = np.where(yy < y0 - thr)[0]; return tt[idx[0]] if len(idx) else tt[0]
        shift = _onset(tM, mm) - _onset(tm_s, Mm)
    sh = lambda tt, _s=shift: (np.asarray(tt) - _s) / 3600.0
    alive = np.asarray(hd.mass_fluid) > 0.02
    cond = np.asarray(hd.m_liquid) + np.asarray(hd.m_solid)

    fig, ax = plt.subplots(2, 2, figsize=(7.4, 5.6))

    # Pressure
    a = ax[0, 0]
    a.plot(sh(tP), ps, color=SLATE, lw=1.6, label="measured")
    a.plot(tm, hd.P / 1e5, color=RED, lw=1.6, ls="--", label="model")
    a.set_ylabel("Pressure [bar]"); a.set_title("Pressure"); a.legend(); a.grid(alpha=.3)
    # Inventory
    a = ax[0, 1]
    a.plot(sh(tM), mm, color=SLATE, lw=1.6, label="measured")
    a.plot(tm, hd.mass_fluid, color=RED, lw=1.6, ls="--", label="model total")
    a.plot(tm, hd.m_solid, color="k", lw=1, ls=":", label="dry ice")
    a.plot(tm, hd.m_liquid, color=NAVY, lw=1, ls="-.", label="liquid")
    a.plot(tm, hd.m_gas, color=AMBER, lw=1, label="gas")
    a.set_ylabel(r"CO$_2$ mass [kg]"); a.set_title("Inventory"); a.legend(); a.grid(alpha=.3)
    # Internal temperature
    a = ax[1, 0]
    tIT, Tc = vcr.merge_across_gaps(IT.iloc[:, 0].values, IT.iloc[:, 1:7].values.T, destep=True)
    Tlo = Tc[0:3]; Tup = Tc[3:6]
    a.fill_between(sh(tIT), np.nanmin(Tup, 0), np.nanmax(Tup, 0), color=AMBER, alpha=.25, label="meas gas band")
    a.fill_between(sh(tIT), np.nanmin(Tlo, 0), np.nanmax(Tlo, 0), color=SLATE, alpha=.40, label="meas liquid/solid band")
    a.plot(tm, np.where(alive, np.asarray(hd.T_gas) - 273.15, np.nan), color=RED, lw=1.5, ls="--", label="model gas")
    a.plot(tm, np.where(alive & (cond > 0.1), np.asarray(hd.T_liquid) - 273.15, np.nan), color=NAVY, lw=1.5, ls="-.", label="model liquid/solid")
    a.set_ylabel(r"Internal T [$^\circ$C]"); a.set_xlabel("Time [h]"); a.set_title("Internal temperature"); a.legend(); a.grid(alpha=.3)
    # Inside-wall temperature
    a = ax[1, 1]
    tWT, Wc = vcr.merge_across_gaps(WT.iloc[:, 0].values, WT.iloc[:, 1:8].values.T, destep=True)
    Wlo = Wc[0:4]; Wup = Wc[4:7]
    a.fill_between(sh(tWT), np.nanmin(Wup, 0), np.nanmax(Wup, 0), color=AMBER, alpha=.20, label="meas gas-wall band")
    a.fill_between(sh(tWT), np.nanmin(Wlo, 0), np.nanmax(Wlo, 0), color=GREEN, alpha=.28, label="meas wetted-wall band")
    a.plot(tm, hd.T_vessel - 273.15, color=RED, lw=1.5, ls="--", label="model gas-wall")
    a.plot(tm, hd.T_vessel_wetted - 273.15, color=NAVY, lw=1.5, ls="-.", label="model wetted wall")
    a.set_ylabel(r"Inside-wall T [$^\circ$C]"); a.set_xlabel("Time [h]"); a.set_title("Wall temperature"); a.legend(loc=WALL_LEGEND_LOC.get(t, "best")); a.grid(alpha=.3)

    xmax = XMAX.get(t, ET)
    for row in ax:
        for a in row:
            a.set_xlim(0, xmax)
    fig.tight_layout()
    os.makedirs(FIGDIR, exist_ok=True)
    fig.savefig(FIGDIR + "cardice_t%d.pdf" % t)
    fig.savefig(FIGDIR + "cardice_t%d.png" % t, dpi=200)  # PNG for quick inspection
    plt.close(fig)
    print("wrote paper/figures/cardice_t%d.{pdf,png}  (align shift %.0f s)" % (t, shift), flush=True)


if __name__ == "__main__":
    for t, pfx in TESTS.items():
        plot_test(t, pfx)
    print("DONE", flush=True)
