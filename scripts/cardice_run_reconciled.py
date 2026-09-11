"""Run the reconciled CARDICE discharge sensitivity set and overlay the raw 1 Hz Ineris
data: pressure, vessel inventory, discharge rate, internal fluid temperatures, and
inner/outer wall temperatures. Confirms the true-orifice + split-Cd (Cd_gas 0.89 /
Cd_liquid 0.62) + calibrated-N parameterisation reproduces the measured blowdown across
ALL observables.

Data source: background/CO2_blowdown-DJa.zip (folders cardice-05..10), extracted to a
DATADIR of E??-{mass-flowrate,pressures,internal-temperatures,inside-wall,outside-wall}
xlsx. Column conventions: internal Tc in 1..6 = bottom(liquid)->top(gas); wall Tc1..7 /
Tc F1..7 = bottom(wetted)->top(gas-contact).

Usage:  python run_reconciled.py [DATADIR]
Outputs CARDICE_test?_reconciled.pdf next to this script.
"""
import os
import re
import sys
import glob
import zipfile
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import yaml

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, os.path.join(REPO, "src"))
from hyddown.hdclass import HydDown

# Reconciled inputs, the output figures and the local data cache live together here:
RECON = os.path.join(REPO, "validation", "reconciled")
# Raw 1 Hz data lives in background/CO2_blowdown-DJa.zip (not committed - ~160 MB). If no
# DATADIR is given, extract the needed E??-* xlsx from that zip into a local cache once.
ZIP = os.path.join(REPO, "background", "CO2_blowdown-DJa.zip")
DEFAULT_CACHE = os.path.join(RECON, "_data")


def ensure_data(datadir):
    if datadir:
        return datadir
    os.makedirs(DEFAULT_CACHE, exist_ok=True)
    if not glob.glob(os.path.join(DEFAULT_CACHE, "E*mass-flowrate*.xlsx")):
        if not os.path.exists(ZIP):
            sys.exit(f"Raw data not found: place CO2_blowdown-DJa.zip in {os.path.dirname(ZIP)} "
                     f"or pass a DATADIR of extracted E??-* xlsx files.")
        want = re.compile(r"cardice-(0[5-9]|10)/.*"
                          r"(mass-flowrate|pressures|internal-temp|inside-wall|outside-wall).*\.xlsx$", re.I)
        with zipfile.ZipFile(ZIP) as z:
            for m in z.namelist():
                if want.search(m):
                    with open(os.path.join(DEFAULT_CACHE, os.path.basename(m)), "wb") as f:
                        f.write(z.read(m))
    return DEFAULT_CACHE


DATA = ensure_data(sys.argv[1] if len(sys.argv) > 1 else None)

NAVY, RED, AMBER, SLATE, GREY = "#002D40", "#D61F39", "#E6A740", "#82979F", "#4C4D4E"
plt.rcParams.update({"font.family": "Arial", "font.size": 8.5})

# test -> (yaml, data-file prefix E??, is_gas)
CASES = {
    "test5":  ("CARDICE_test5.yml",  "05", True),
    "test6":  ("CARDICE_test6.yml",  "6",  False),
    "test7":  ("CARDICE_test7.yml",  "07", True),
    "test8":  ("CARDICE_test8.yml",  "08", False),
    "test9":  ("CARDICE_test9.yml",  "09", True),
    "test10": ("CARDICE_test10.yml", "10", False),
}
K = 273.15


def find(prefix, kind):
    hits = glob.glob(os.path.join(DATA, f"E{prefix}-*{kind}*.xlsx"))
    if not hits:
        raise FileNotFoundError(f"E{prefix} {kind}")
    return hits[0]


def meas_rate(t, M):
    r = np.full_like(t, np.nan, dtype=float)
    for i in range(0, len(t), 3):
        k = (t >= t[i] - 30) & (t <= t[i] + 30)
        if k.sum() >= 10:
            r[i] = -np.polyfit(t[k], M[k], 1)[0]
    return r


def C(model_arr):
    a = np.array(model_arr, dtype=float) - K
    a[a < -120] = np.nan       # mask unfilled/parked entries
    return a


def wall(detailed, lumped):
    """Continuous wall temperature across the CoolProp->thermopack handover: the detailed
    inner/outer-face node while above the triple point, then the two-zone lumped wall
    node (T_vessel / T_vessel_wetted) below it, where the detailed node is no longer
    updated (it parks at 0 K)."""
    d = np.array(detailed, dtype=float) - K
    l = np.array(lumped, dtype=float) - K
    out = np.where(d < -120, l, d)
    out[out < -120] = np.nan
    return out


def run_model(yml):
    with open(os.path.join(RECON, yml)) as f:
        inp = yaml.safe_load(f)
    hd = HydDown(inp)
    hd.run(disable_pbar=True)
    return hd, inp["release"]


def title_bits(rel):
    """Gas/liquid Cd + N read from the release block, for an honest figure title."""
    cd_gas = rel.get("discharge_coef_gas", rel["discharge_coef"])
    if rel["type"] == "liquid":
        return f"Cd_gas {cd_gas:g} / Cd_liq {rel['discharge_coef']:g} + N {rel.get('liquid_nonequilibrium', 0):g}"
    return f"Cd_gas {cd_gas:g} (gas release)"


def main():
    for tid, (yml, pfx, is_gas) in CASES.items():
        hd, rel = run_model(yml)
        th = hd.time_array / 3600.0

        m = pd.read_excel(find(pfx, "mass-flowrate"))
        tm = m.iloc[:, 0].values.astype(float); M = m.iloc[:, 1].values.astype(float)
        # Time-shift ALL channels to the actual start of discharge: the onset of the
        # load-cell mass decline (first 2 kg drop below the pre-discharge plateau). This
        # removes the valve-opening dead period present in every record - gas tests carry
        # a long pre-blowdown hold, liquid tests a short ramp (e.g. T6) - so the measured
        # discharge start aligns with the model's t = 0.
        M0 = np.median(M[:50])
        below = np.where(M < M0 - 2.0)[0]
        t0 = tm[below[0]] if len(below) else 0.0
        p = pd.read_excel(find(pfx, "pressures"))
        internal = pd.read_excel(find(pfx, "internal-temp"))
        iwall = pd.read_excel(find(pfx, "inside-wall"))
        owall = pd.read_excel(find(pfx, "outside-wall"))
        rM = meas_rate(tm, M)

        def X(df):  # time in hours, shifted to blowdown start
            return (df.iloc[:, 0].values.astype(float) - t0) / 3600.0

        fig, ax = plt.subplots(2, 3, figsize=(13, 7))
        # 1 pressure
        ax[0, 0].plot(X(p), p["P sphere"], color=NAVY, lw=1.1, label="measured")
        ax[0, 0].plot(th, hd.P / 1e5, color=RED, lw=1.7, ls="--", label="reconciled")
        ax[0, 0].set_ylabel("pressure [bar]"); ax[0, 0].set_title("Pressure")
        # 2 inventory
        ax[0, 1].plot(X(m), M, color=NAVY, lw=1.1, label="measured (load cells)")
        ax[0, 1].plot(th, hd.mass_fluid, color=RED, lw=1.7, ls="--", label="reconciled total")
        ax[0, 1].plot(th, hd.m_liquid, color=SLATE, lw=0.9, label="liquid")
        ax[0, 1].plot(th, hd.m_gas, color=AMBER, lw=0.9, label="gas")
        ax[0, 1].plot(th, hd.m_solid, color=GREY, lw=0.9, ls=":", label="solid (dry ice)")
        ax[0, 1].set_ylabel("mass [kg]"); ax[0, 1].set_title("Inventory")
        # 3 discharge rate
        fin = np.isfinite(rM)
        ax[0, 2].plot(X(m)[fin], rM[fin], color=SLATE, lw=0.7, alpha=0.7, label="measured (raw)")
        ax[0, 2].plot(th, hd.mass_rate, color=RED, lw=1.7, ls="--", label="reconciled")
        ax[0, 2].set_ylabel("discharge rate [kg/s]"); ax[0, 2].set_title("Discharge rate")
        rr = rM[(X(m) > 0.02)]
        ymax = np.nanpercentile(rr, 99) if np.isfinite(rr).any() else 0.5
        ax[0, 2].set_ylim(-0.02, max(0.06, ymax * 1.35))
        # 4 fluid temps: 6 internal Tc (bottom->top) vs model gas / liquid-solid
        for j, c in enumerate(internal.columns[1:]):
            ax[1, 0].plot(X(internal), internal[c], color=SLATE, lw=0.6, alpha=0.6,
                          label="measured Tc in 1..6" if j == 0 else None)
        ax[1, 0].plot(th, C(hd.T_gas), color=RED, lw=1.7, ls="--", label="model gas")
        ax[1, 0].plot(th, C(hd.T_liquid), color=NAVY, lw=1.7, ls="--", label="model liquid/solid")
        ax[1, 0].set_ylabel("fluid T [$^\\circ$C]"); ax[1, 0].set_title("Internal fluid temperature")
        # 5 inner wall: 7 Tc (bottom wetted -> top gas-contact) vs model gas-contact / wetted
        for j, c in enumerate(iwall.columns[1:]):
            ax[1, 1].plot(X(iwall), iwall[c], color=SLATE, lw=0.6, alpha=0.6,
                          label="measured Tc1..7" if j == 0 else None)
        ax[1, 1].plot(th, wall(hd.T_inner_wall, hd.T_vessel), color=RED, lw=1.7, ls="--", label="model gas-contact")
        ax[1, 1].plot(th, wall(hd.T_inner_wall_wetted, hd.T_vessel_wetted), color=NAVY, lw=1.7, ls="--", label="model wetted")
        ax[1, 1].set_ylabel("inner-wall T [$^\\circ$C]"); ax[1, 1].set_title("Inner-wall temperature")
        # 6 outer wall
        for j, c in enumerate(owall.columns[1:]):
            ax[1, 2].plot(X(owall), owall[c], color=SLATE, lw=0.6, alpha=0.6,
                          label="measured Tc F1..7" if j == 0 else None)
        ax[1, 2].plot(th, wall(hd.T_outer_wall, hd.T_vessel), color=RED, lw=1.7, ls="--", label="model gas-contact")
        ax[1, 2].plot(th, wall(hd.T_outer_wall_wetted, hd.T_vessel_wetted), color=NAVY, lw=1.7, ls="--", label="model wetted")
        ax[1, 2].set_ylabel("outer-wall T [$^\\circ$C]"); ax[1, 2].set_title("Outer-wall temperature")

        for a in ax.flat:
            a.set_xlabel("time [h]"); a.grid(alpha=0.25); a.legend(fontsize=6.5, loc="best")
            a.set_xlim(-0.02, th[-1] * 1.02)
        fig.suptitle(f"CARDICE {tid} - reconciled (true orifice, {title_bits(rel)}) vs 1 Hz data",
                     color=NAVY, fontweight="bold", fontsize=11)
        fig.tight_layout(rect=[0, 0, 1, 0.97])
        out = os.path.join(RECON, f"CARDICE_{tid}_reconciled.pdf")
        fig.savefig(out); plt.close(fig)
        print(f"{tid}: P0={hd.P[0]/1e5:.2f}b m0={hd.mass_fluid[0]:.0f}kg "
              f"solid_end={hd.m_solid[-1]:.0f}kg dur={th[np.argmax(hd.P/1e5 < 2.0) if (hd.P/1e5<2).any() else -1]:.2f}h -> {os.path.basename(out)}")


if __name__ == "__main__":
    main()
