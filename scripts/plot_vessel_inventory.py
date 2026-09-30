"""Vessel inventory per phase (gas / liquid / dry-ice) over time for a QRA case.

Shows the bottom-outlet drainage, the residual liquid at the triple point and any dry ice
formed in the vessel. Usage: python scripts/plot_vessel_inventory.py [case_basename]
(default A2_P11_10in_liquid). Writes <case>_inventory.pdf to the OneDrive Source terms folder.
"""
import os, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

matplotlib.rcParams["font.family"] = "Arial"
matplotlib.rcParams["font.size"] = 10
NAVY, RED, AMBER, SLATE, GREY = "#002D40", "#D61F39", "#E6A740", "#82979F", "#4C4D4E"

HERE = os.path.dirname(os.path.abspath(__file__))
CSV_DIR = os.path.join(HERE, "..", "qra_output")
OUT_DIR = os.environ.get("QRA_PLOT_DIR",
          ("C:\\Users\\AndersAndreasen\\OneDrive - ORS\\Projects 2026 - "
           "120.691_D_CarbonCuts_Ruby Development Project\\01 Working area\\CFD\\Source terms"))
P_TRIPLE = 5.18  # bar


def plot(case):
    df = pd.read_csv(os.path.join(CSV_DIR, case + ".csv"))
    t = df["Time (s)"].values
    mg = df["Vessel gas mass (kg)"].values
    ml = df["Vessel liquid mass (kg)"].values
    ms = df["Vessel dry-ice mass (kg)"].values if "Vessel dry-ice mass (kg)" in df else np.zeros_like(t)
    P = df["Pressure (bar)"].values
    itr = int(np.argmax(P <= P_TRIPLE)) if np.any(P <= P_TRIPLE) else -1

    fig, ax = plt.subplots(1, 2, figsize=(12, 4.8))
    fig.suptitle("Vessel inventory per phase  --  %s" % case, fontsize=12, fontweight="bold", color=NAVY)

    for a in ax:
        a.plot(t, ml, color=NAVY, lw=2.0, label="Liquid")
        a.plot(t, mg, color=RED, lw=2.0, label="Gas")
        a.plot(t, ms, color=AMBER, lw=2.0, label="Dry ice (solid)")
        if itr > 0:
            a.axvline(t[itr], color=SLATE, ls="--", lw=1.2, label="triple point (%.1f bar)" % P_TRIPLE)
        a.set_xlabel("Time [s]"); a.grid(True, color=SLATE, alpha=0.25, lw=0.5)
        a.set_xlim(t.min(), t.max())
    ax[0].set_ylabel("Vessel mass [kg]"); ax[0].set_title("(a) Full inventory", color=NAVY)
    ax[0].legend(frameon=False, fontsize=8)
    # zoom on the tail so the small liquid heel / dry ice near the triple point is visible
    ax[1].set_yscale("log"); ax[1].set_ylabel("Vessel mass [kg]  (log)")
    ax[1].set_title("(b) Log scale (residual heel + dry ice)", color=NAVY)
    ax[1].legend(frameon=False, fontsize=8)

    fig.tight_layout(rect=(0, 0, 1, 0.95))
    os.makedirs(OUT_DIR, exist_ok=True)
    pdf = os.path.join(OUT_DIR, case + "_inventory.pdf")
    fig.savefig(pdf); plt.close(fig)
    # quick console summary
    if itr > 0:
        print("at triple point (t=%.0fs): liquid=%.1f kg  gas=%.1f kg  dry-ice=%.1f kg"
              % (t[itr], ml[itr], mg[itr], ms[itr]))
    print("final: liquid=%.1f  gas=%.1f  dry-ice=%.1f" % (ml[-1], mg[-1], ms[-1]))
    print("wrote %s" % pdf)


if __name__ == "__main__":
    plot(sys.argv[1] if len(sys.argv) > 1 else "A2_P11_10in_liquid")
