"""Build one xlsx per QRA source-term case from the HydDown CSVs.

Each xlsx has a 'Source term' sheet (time series, temperatures in K) with the in-vessel /
isolatable-volume conditions, the choked-throat conditions (incl. pressure) and the
fully-expanded atmospheric conditions, plus an 'Info' sheet with the case parameters.

Reads CSVs from qra_output/ and writes xlsx to the OneDrive 'Source terms' folder.
"""
import os, glob
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
CSV_DIR = os.path.join(HERE, "..", "qra_output")
XLSX_DIR = ("C:\\Users\\AndersAndreasen\\OneDrive - ORS\\Projects 2026 - "
            "120.691_D_CarbonCuts_Ruby Development Project\\01 Working area\\CFD\\Source terms")
os.makedirs(XLSX_DIR, exist_ok=True)

# base-case metadata (P as given: barg; T in C). Keyed by CSV name without the _liquid/_mix suffix.
META = {
    "A1_P11_6in":       ("A-1",      "P11 storage tank (liquid outlet)",  1100.0, 15.2, -26.1, "6\"",  "-",      3.67e-4, "saturated (two-phase, 90% liq) + static head"),
    "A1_P5_8in":        ("A-1",      "P5 LCO2 combined offloading",         100.0, 16.5, -26.0, "8\"",  "Sch 40", 1.20e-5, "sub-cooled liquid"),
    "A_P4_8in":         ("A-1/A-2",  "P4 LCO2 header truck unloading",       15.6, 18.8, -26.0, "8\"",  "Sch 40", 1.83e-5, "sub-cooled liquid"),
    "A2_P11_10in":      ("A-2/B-2",  "P11 storage tank (liquid outlet)",   1100.0, 15.2, -26.1, "10\"", "-",      1.50e-3, "saturated (two-phase, 90% liq) + static head"),
    "A2_P6_10in":       ("A-2/B-2",  "P6 booster pump P-501A/B",            100.0, 22.2, -25.6, "10\"", "Sch 40", 1.52e-4, "sub-cooled liquid"),
    "B2_P8_8in_sch40":  ("B-2",      "P8 injection heater E-501 outlet",    100.0, 101.0, 0.0,  "8\"",  "Sch 40", 2.64e-4, "dense/sub-cooled (above critical P)"),
    "B2_P8_8in_sch120": ("B-2",      "P8 injection heater E-501 outlet",    100.0, 101.0, 0.0,  "8\"",  "Sch 120",2.64e-4, "dense/sub-cooled (above critical P)"),
}
K = 273.15


def build(csv_path):
    name = os.path.splitext(os.path.basename(csv_path))[0]
    rtype = "mix" if name.endswith("_mix") else "liquid"
    base = name[:-4] if name.endswith("_mix") else name[:-7]
    df = pd.read_csv(csv_path)
    out = pd.DataFrame()
    out["Time (s)"] = df["Time (s)"]
    # In-vessel / isolatable volume
    out["Vessel pressure (bar)"] = df["Pressure (bar)"]
    out["Vessel fluid temperature (K)"] = df["Fluid temperature (oC)"] + K
    if "Vessel gas temperature (oC)" in df:
        out["Vessel gas temperature (K)"] = df["Vessel gas temperature (oC)"] + K
        out["Vessel liquid/solid temperature (K)"] = df["Vessel liquid/solid temperature (oC)"] + K
    out["Vessel solid mass fraction (-)"] = df["Vessel dry-ice mass fraction (-)"]
    # Discharge mass flow (same at vessel / throat / atmosphere)
    out["Mass flow (kg/s)"] = df["Release mass rate (kg/s)"]
    # Choked-throat conditions
    out["Throat pressure (bar)"] = df["Throat pressure (bar)"]
    out["Throat temperature (K)"] = df["Throat temperature (oC)"] + K
    if "Throat density (kg/m3)" in df:
        out["Throat density (kg/m3)"] = df["Throat density (kg/m3)"]
        out["Throat velocity (m/s)"] = df["Throat velocity (m/s)"]
    out["Throat solid mass fraction (-)"] = df["Throat dry-ice mass fraction (-)"]
    # Fully-expanded atmospheric conditions
    out["Atmospheric temperature (K)"] = df["Atmospheric temperature (oC)"] + K
    out["Atmospheric solid mass fraction (-)"] = df["Atmospheric dry-ice mass fraction (-)"]

    m = META.get(base, ("?", base, np.nan, np.nan, np.nan, "?", "?", np.nan, "?"))
    info = pd.DataFrame({
        "Field": ["Scenario", "Segment", "Isolatable volume (m3)", "Pressure (barg)",
                  "Temperature (C)", "Line size", "Schedule", "Rupture frequency (/yr)",
                  "Initial state", "Release model", "Discharge coefficient Cd",
                  "Liquid non-equilibrium N", "Total released (kg)", "Peak mass flow (kg/s)",
                  "Modelled duration (s)"],
        "Value": [m[0], m[1], m[2], m[3], m[4], m[5], m[6], m[7], m[8],
                  "%s (HEM, adiabatic)" % rtype, 1.0, 0.0,
                  round(float(np.trapz(out["Mass flow (kg/s)"], out["Time (s)"])), 1),
                  round(float(out["Mass flow (kg/s)"].max()), 1),
                  round(float(out["Time (s)"].iloc[-1]), 1)],
    })

    xlsx = os.path.join(XLSX_DIR, name + ".xlsx")
    with pd.ExcelWriter(xlsx, engine="openpyxl") as xw:
        out.to_excel(xw, sheet_name="Source term", index=False)
        info.to_excel(xw, sheet_name="Info", index=False)
    return name, len(out)


def main():
    n = 0
    for csv in sorted(glob.glob(os.path.join(CSV_DIR, "*.csv"))):
        name, rows = build(csv)
        print("  %-26s -> %s.xlsx  (%d rows)" % (name, name, rows))
        n += 1
    print("wrote %d xlsx to %s" % (n, XLSX_DIR))


if __name__ == "__main__":
    main()
