"""Re-plot the two Munkejord/SINTEF dense-phase tests used in the main paper (Exp72 and Exp53)
with the SciencePlots 'science' + 'nature' styles and LaTeX text rendering, so CO2 is set as
CO$_2$. Axis labels are capitalised and moderately enlarged and the legend text is enlarged,
matching the CARDICE paper figures (see plot_cardice_scienceplots.py).

Data loading (SINTEF 1 Hz zip), the HydDown run and the per-test time window are reused unchanged
from ``validate_munkejord.py``; only the presentation differs. The measured temperatures are shown
as a single minimum-to-maximum band over the six vertical thermocouples plus the median. Output is
written to paper/figures/munke_exp{72,53}.pdf (+ .png for quick inspection).

Requires: scienceplots, and a LaTeX toolchain (latex + dvipng) for text.usetex.

Run from the repo root:  python scripts/plot_munkejord_scienceplots.py
"""
import os, sys, zipfile
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import scienceplots  # noqa: F401  registers the 'science' and 'nature' styles

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import validate_munkejord as vm  # data / model / time-window helpers
from hyddown.hdclass import HydDown  # noqa: F401  (imported by vm; kept explicit)

NAVY, RED, AMBER, SLATE, GREEN = vm.NAVY, vm.RED, vm.AMBER, vm.SLATE, vm.GREEN
FIGDIR = os.path.join(vm.REPO, "paper", "figures") + os.sep

TESTS = {name: "munke_" + name.lower() for (name, *_rest) in vm.TESTS}  # all nine SINTEF tests (SI)
WALL_LEGEND_LOC = {}   # per-test wall-panel legend location, e.g. {"Exp72": "lower left"}; default "best"

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


def specs(name):
    for n, P0, T0, noz, riser in vm.TESTS:
        if n == name:
            return P0, T0, noz, riser
    raise KeyError(name)


def plot_test(z, name, base):
    P0, T0, noz, riser = specs(name)
    ET = vm.endtime(name, noz)
    g, Pavg, W, fluid, iwall = vm.measured(z, name, ET)
    hd = vm.run_model(P0, T0, noz, riser, ET)
    t = np.asarray(hd.time_array)
    alive = np.asarray(hd.mass_fluid) > 0.02
    cond = np.asarray(hd.m_liquid) + np.asarray(hd.m_solid)

    fig, ax = plt.subplots(2, 2, figsize=(7.4, 5.6))
    # Pressure
    a = ax[0, 0]
    a.plot(g, Pavg, color=SLATE, lw=1.8, label="measured")
    a.plot(t, hd.P / 1e5, color=RED, lw=1.6, ls="--", label="model")
    a.set_ylabel("Pressure [bar]"); a.set_title("Pressure"); a.legend(); a.grid(alpha=.3)
    # Inventory
    a = ax[0, 1]
    a.plot(g, W, color=SLATE, lw=1.8, label="measured")
    a.plot(t, hd.mass_fluid, color=RED, lw=1.6, ls="--", label="model total")
    a.plot(t, hd.m_liquid, color=NAVY, lw=1, ls="-.", label="liquid")
    a.plot(t, hd.m_gas, color=AMBER, lw=1, label="gas")
    a.plot(t, hd.m_solid, color="k", lw=1, ls=":", label="dry ice")
    a.set_ylabel(r"CO$_2$ mass [kg]"); a.set_title("Inventory"); a.legend(); a.grid(alpha=.3)
    # Fluid temperature
    a = ax[1, 0]
    a.fill_between(g, fluid.min(0), fluid.max(0), color=SLATE, alpha=.30, label="measured band")
    a.plot(g, np.median(fluid, 0), color=SLATE, lw=1.4, label="measured midline")
    a.plot(t, np.where(alive, np.asarray(hd.T_gas) - 273.15, np.nan), color=RED, lw=1.6, ls="--", label="model gas")
    a.plot(t, np.where(alive & (cond > 0.1), np.asarray(hd.T_liquid) - 273.15, np.nan), color=NAVY, lw=1.6, ls="-.", label="model liquid/solid")
    a.set_ylabel(r"Fluid T [$^\circ$C]"); a.set_xlabel("Time [s]"); a.set_title("Fluid temperature"); a.legend(); a.grid(alpha=.3)
    # Wall temperature
    a = ax[1, 1]
    a.fill_between(g, iwall.min(0), iwall.max(0), color=GREEN, alpha=.25, label="measured band")
    a.plot(g, np.median(iwall, 0), color=GREEN, lw=1.4, label="measured midline")
    a.plot(t, hd.T_vessel - 273.15, color=RED, lw=1.6, ls="--", label="model gas-wall")
    a.plot(t, hd.T_vessel_wetted - 273.15, color=NAVY, lw=1.6, ls="-.", label="model wetted wall")
    a.set_ylabel(r"Inside-wall T [$^\circ$C]"); a.set_xlabel("Time [s]"); a.set_title("Wall temperature")
    a.legend(loc=WALL_LEGEND_LOC.get(name, "best")); a.grid(alpha=.3)

    for row in ax:
        for a in row:
            a.set_xlim(0, ET)
    fig.tight_layout()
    os.makedirs(FIGDIR, exist_ok=True)
    fig.savefig(FIGDIR + base + ".pdf")
    fig.savefig(FIGDIR + base + ".png", dpi=200)  # PNG for quick inspection
    plt.close(fig)
    print("wrote paper/figures/%s.{pdf,png}  (end %.0f s)" % (base, ET), flush=True)


if __name__ == "__main__":
    z = zipfile.ZipFile(os.path.join(vm.REPO, "background", "19589510.zip"))
    for name, base in TESTS.items():
        plot_test(z, name, base)
    print("DONE", flush=True)
