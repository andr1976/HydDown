"""Per-case source-term plots for the QRA / KFX CFD deliverable.

Reads the HydDown CSVs from qra_output/ and writes one multi-panel PDF per case to the
OneDrive 'Source terms' folder, alongside the CSVs and xlsx. Panels: pressures, release
mass flow, temperatures, dry-ice mass fractions, throat density and throat velocity.
"""
import os, glob
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

matplotlib.rcParams["font.family"] = "Arial"
matplotlib.rcParams["font.size"] = 9
matplotlib.rcParams["axes.titlesize"] = 10

NAVY, RED, AMBER, SLATE, GREY = "#002D40", "#D61F39", "#E6A740", "#82979F", "#4C4D4E"

HERE = os.path.dirname(os.path.abspath(__file__))
CSV_DIR = os.path.join(HERE, "..", "qra_output")
OUT_DIR = ("C:\\Users\\AndersAndreasen\\OneDrive - ORS\\Projects 2026 - "
           "120.691_D_CarbonCuts_Ruby Development Project\\01 Working area\\CFD\\Source terms")
os.makedirs(OUT_DIR, exist_ok=True)


def plot_one(csv_path):
    name = os.path.splitext(os.path.basename(csv_path))[0]
    df = pd.read_csv(csv_path)
    t = df["Time (s)"].values

    fig, ax = plt.subplots(3, 2, figsize=(11, 11))
    fig.suptitle("QRA source term  --  %s" % name, fontsize=12, fontweight="bold", color=NAVY)

    # (a) pressures
    a = ax[0, 0]
    a.plot(t, df["Pressure (bar)"], color=NAVY, lw=1.8, label="Vessel")
    if "Throat pressure (bar)" in df:
        a.plot(t, df["Throat pressure (bar)"], color=SLATE, lw=1.4, ls="--", label="Throat")
    a.set_ylabel("Pressure [bar]"); a.set_title("(a) Pressure"); a.legend(frameon=False)

    # (b) release mass flow
    a = ax[0, 1]
    a.plot(t, df["Release mass rate (kg/s)"], color=RED, lw=1.8)
    a.set_ylabel("Mass flow [kg/s]"); a.set_title("(b) Release mass flow")

    # (c) temperatures (degC)
    a = ax[1, 0]
    a.plot(t, df["Fluid temperature (oC)"], color=NAVY, lw=1.6, label="Vessel fluid")
    if "Vessel gas temperature (oC)" in df:
        a.plot(t, df["Vessel gas temperature (oC)"], color=AMBER, lw=1.2, ls="-.", label="Vessel gas")
        a.plot(t, df["Vessel liquid/solid temperature (oC)"], color=SLATE, lw=1.2, ls=":", label="Vessel liq/solid")
    if "Throat temperature (oC)" in df:
        a.plot(t, df["Throat temperature (oC)"], color=RED, lw=1.4, ls="--", label="Throat")
    if "Atmospheric temperature (oC)" in df:
        a.plot(t, df["Atmospheric temperature (oC)"], color=GREY, lw=1.2, label="Atmospheric")
    a.set_ylabel("Temperature [$^\\circ$C]"); a.set_title("(c) Temperatures"); a.legend(frameon=False, fontsize=7)

    # (d) dry-ice mass fractions
    a = ax[1, 1]
    if "Vessel dry-ice mass fraction (-)" in df:
        a.plot(t, df["Vessel dry-ice mass fraction (-)"], color=NAVY, lw=1.6, label="Vessel")
    if "Throat dry-ice mass fraction (-)" in df:
        a.plot(t, df["Throat dry-ice mass fraction (-)"], color=RED, lw=1.4, ls="--", label="Throat")
    if "Atmospheric dry-ice mass fraction (-)" in df:
        a.plot(t, df["Atmospheric dry-ice mass fraction (-)"], color=GREY, lw=1.2, label="Atmospheric")
    a.set_ylabel("Dry-ice mass fraction [-]"); a.set_title("(d) Solid (dry-ice) fraction")
    a.legend(frameon=False, fontsize=7)

    # (e) throat density
    a = ax[2, 0]
    if "Throat density (kg/m3)" in df:
        a.plot(t, df["Throat density (kg/m3)"], color=NAVY, lw=1.6)
    a.set_ylabel("Throat density [kg/m$^3$]"); a.set_xlabel("Time [s]"); a.set_title("(e) Throat density")

    # (f) throat velocity
    a = ax[2, 1]
    if "Throat velocity (m/s)" in df:
        a.plot(t, df["Throat velocity (m/s)"], color=RED, lw=1.6)
    a.set_ylabel("Throat velocity [m/s]"); a.set_xlabel("Time [s]"); a.set_title("(f) Throat velocity")

    for row in ax:
        for a in row:
            a.grid(True, color=SLATE, alpha=0.25, lw=0.5)
            a.set_xlim(t.min(), t.max())

    fig.tight_layout(rect=(0, 0, 1, 0.97))
    pdf = os.path.join(OUT_DIR, name + ".pdf")
    fig.savefig(pdf)
    plt.close(fig)
    return name


def main():
    n = 0
    for csv in sorted(glob.glob(os.path.join(CSV_DIR, "*.csv"))):
        if os.path.basename(csv).startswith("_"):
            continue
        nm = plot_one(csv)
        print("  %-26s -> %s.pdf" % (nm, nm))
        n += 1
    print("wrote %d PDF plots to %s" % (n, OUT_DIR))


if __name__ == "__main__":
    main()
