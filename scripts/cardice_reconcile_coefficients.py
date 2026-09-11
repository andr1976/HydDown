"""Derive the reconciled CARDICE discharge coefficients from the raw 1 Hz Ineris data
and the TRUE orifice sizes (SSRN 4292065 / GHGT-16 Drescher et al., Table 1: 3/4/6 mm).

Scheme
------
  Cd_gas    : single-phase (HEM) gas discharge coefficient, from the gas tests
              (T5/T7/T9). Paired with the measured sphere pressure; verified constant
              along each pressure decline. -> the data gives 0.89 (0.87-0.92); the
              reconciled set adopts the canonical sharp-edged value Cd_gas = 0.84.
  Cd_liquid : two-phase / flashing sharp-orifice coefficient (Darby / API 520), FIXED
              at 0.62 for every liquid discharge (a saturated liquid flashing through a
              short sharp orifice; the paper itself cites Darby for exactly this).
  N (HNE)   : delayed/metastable flashing boost on top of Cd_liquid, calibrated per test
              to the measured steady liquid-drain rate:
                  m_liq = Cd_liquid * A * sqrt((1-N)*G_HEM^2 + N*G_frozen^2)

Rates are the -dM/dt plateau (rolling 61 s slope), skipping the initial valve-opening
dead period so the estimate is not overcompensated by the start-up spike.

Run:  python reconcile_coefficients.py        # extracts data from background/ zip
"""
import os
import re
import sys
import glob
import math
import zipfile
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, os.path.join(REPO, "src"))
from hyddown.co2_release import CO2ReleaseModel

ZIP = os.path.join(REPO, "background", "CO2_blowdown-DJa.zip")
CACHE = os.path.join(REPO, "validation", "reconciled", "_data")
PB = 101325.0
CD_GAS = 0.84       # canonical value adopted by the reconciled set (data gives ~0.89)
CD_LIQUID = 0.62    # fixed Darby / API 520 two-phase sharp-orifice coefficient
A = lambda d_mm: math.pi * (d_mm / 1000.0) ** 2 / 4.0


def ensure_data():
    os.makedirs(CACHE, exist_ok=True)
    if glob.glob(os.path.join(CACHE, "E*mass-flowrate*.xlsx")):
        return
    want = re.compile(r"cardice-(0[5-9]|10)/.*(mass-flowrate|pressures).*\.xlsx$", re.I)
    with zipfile.ZipFile(ZIP) as z:
        for m in z.namelist():
            if want.search(m):
                with open(os.path.join(CACHE, os.path.basename(m)), "wb") as f:
                    f.write(z.read(m))


def find(pfx, kind):
    return glob.glob(os.path.join(CACHE, f"E{pfx}-*{kind}*.xlsx"))[0]


def load(pfx):
    m = pd.read_excel(find(pfx, "mass-flowrate"))
    p = pd.read_excel(find(pfx, "pressures"))
    return (m.iloc[:, 0].values.astype(float), m.iloc[:, 1].values.astype(float),
            p.iloc[:, 0].values.astype(float), p["P sphere"].values.astype(float) * 1e5)


def rate(t, M, tc, h=30):
    k = (t >= tc - h) & (t <= tc + h)
    return -np.polyfit(t[k], M[k], 1)[0] if k.sum() >= 10 else np.nan


def main():
    ensure_data()
    model = CO2ReleaseModel(back_pressure=PB, atm_pressure=PB, liquid_nonequilibrium=0.0)

    print("GAS  ->  Cd_gas = m_gas / (A_true * G_HEM(P))  [equilibrium HEM]")
    print(f"{'test':>4} {'P[bar]':>7} {'m_gas':>8} {'d[mm]':>6} {'Cd_gas':>7}")
    gas = [("T5", "05", 3.0), ("T7", "07", 4.0), ("T9", "09", 4.0)]
    cds = []
    for name, pfx, d in gas:
        t, M, pt, P = load(pfx)
        M0 = np.median(M[:50]); t0 = t[np.argmax(M < M0 - 3)] if (M < M0 - 3).any() else 0.0
        a, b = t0 + 100, t0 + 400
        mg = np.nanmedian([rate(t, M, tc) for tc in np.arange(a, b, 20)])
        Pw = np.median(P[(pt >= a) & (pt <= b)])
        G = model.hem_rate(Pw, "gas", 1.0, 1.0)["G"]
        cd = mg / (A(d) * G); cds.append(cd)
        print(f"{name:>4} {Pw/1e5:7.2f} {mg:8.4f} {d:6.1f} {cd:7.3f}")
    print(f"  -> Cd_gas mean {np.mean(cds):.3f}  (range {min(cds):.2f}-{max(cds):.2f});  used: {CD_GAS}\n")

    print(f"LIQUID ->  N calibrated to steady drain, Cd_liquid fixed at {CD_LIQUID}")
    print(f"{'test':>4} {'P0[bar]':>8} {'m_liq':>8} {'d[mm]':>6} {'G_HEM':>8} {'G_froz':>8} {'N':>6}")
    liq = [("T6", "6", 3.0, (400, 1800)), ("T8", "08", 4.0, (400, 2200)),
           ("T10", "10", 4.0, (400, 2000))]
    for name, pfx, d, win in liq:
        t, M, pt, P = load(pfx)
        ml = float(np.nanmedian([rate(t, M, tc) for tc in np.arange(*win, 20)]))
        P0 = np.median(P[(pt >= 0) & (pt <= 30)])
        h0, s0, rho0, T0 = model.stagnation(P0, "liquid")
        G_HEM = model.hem_rate_from_stagnation(h0, s0, rho0, T0, P0, 1.0, 1.0)["G"]
        G_frozen = math.sqrt(2.0 * rho0 * max(P0 - PB, 0.0))
        Gb = ml / (CD_LIQUID * A(d))
        N = max((Gb ** 2 - G_HEM ** 2) / (G_frozen ** 2 - G_HEM ** 2), 0.0)
        print(f"{name:>4} {P0/1e5:8.2f} {ml:8.4f} {d:6.1f} {G_HEM:8.0f} {G_frozen:8.0f} {N:6.3f}")


if __name__ == "__main__":
    main()
