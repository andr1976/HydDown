"""QRA source-term runs for the DNV KFX CFD analysis.

Builds and runs the CO2 depressurisation cases (see the scenario overview) and writes
one CSV per case (HydDown get_dataframe, in oC / bar) for the xlsx post-processing.

Vessel model: adiabatic external (h_outer = 0), NEM solver, HEM release. The storage
tank (P11) is a vertical hemispherical-headed cylinder with a bottom liquid outlet and
the time-dependent static head; the pipe/equipment segments are horizontal Sch-40 (or
Sch-120) lines with the static head off. Pipes are run with both a liquid and a mix
(bulk two-phase) release; the tank with a liquid release + head.

Usage:  python scripts/qra_source_terms.py [name_substring]   (default: all)
"""
import os, sys, copy, math
import numpy as np
import fluids

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src"))
from hyddown.hdclass import HydDown

OUT = os.path.join(HERE, "..", "qra_output")
os.makedirs(OUT, exist_ok=True)

# Pipe internal diameter [m] and wall [m] by (nominal, schedule)
PIPE_ID = {("6", 40): 0.15405, ("8", 40): 0.20272, ("10", 40): 0.25451, ("8", 120): 0.18256}
PIPE_WALL = {("6", 40): 0.00711, ("8", 40): 0.00818, ("10", 40): 0.00927, ("8", 120): 0.01826}
RHO_W, CP_W = 7850.0, 490.0          # LTCS wall
TANK_D, TANK_L, TANK_T = 6.785, 25.9, 0.039   # vertical, hemispherical heads, ~39 mm shell

# name, segment, V[m3], P[bara], T[C], nominal, sched, kind, releases, dt, end_time[s]
CASES = [
    ("A1_P11_6in",  "P11", 1100, 16.2, -26.1, "6",  40, "tank", ["liquid"],        0.5, 5000),
    ("A1_P5_8in",   "P5",  100,  17.5, -26.0, "8",  40, "pipe", ["liquid", "mix"], 0.1, 500),
    ("A_P4_8in",    "P4",  15.6, 19.8, -26.0, "8",  40, "pipe", ["liquid", "mix"], 0.05, 120),
    ("A2_P11_10in", "P11", 1100, 16.2, -26.1, "10", 40, "tank", ["liquid"],        0.5, 2500),
    ("A2_P6_10in",  "P6",  100,  23.2, -25.6, "10", 40, "pipe", ["liquid", "mix"], 0.1, 500),
    ("B2_P8_8in_sch40",  "P8", 100, 102.0, 0.0, "8", 40,  "pipe", ["liquid", "mix"], 0.05, 500),
    ("B2_P8_8in_sch120", "P8", 100, 102.0, 0.0, "8", 120, "pipe", ["liquid", "mix"], 0.05, 500),
]


def level_for_fill(D, L, vtype, horiz, fill):
    if vtype == "Hemispherical":
        t = fluids.TANK(D=D, L=L, sideA="spherical", sideB="spherical",
                        sideA_a=0.5 * D, sideB_a=0.5 * D, horizontal=horiz)
    else:
        t = fluids.TANK(D=D, L=L, horizontal=horiz)
    return t.h_from_V(fill * t.V_total), t.V_total


def build(name, seg, V, P_bara, T_C, nominal, sched, kind, rtype, dt, end_time):
    T = T_C + 273.15
    ID = PIPE_ID[(nominal, sched)]
    if kind == "tank":
        # Saturated storage tank: two-phase, 90% liquid fill, bottom outlet + static head.
        D, L, vtype, horiz, thick = TANK_D, TANK_L, "Hemispherical", False, TANK_T
        level, _ = level_for_fill(D, L, vtype, horiz, 0.90)
        static, orient = True, "vertical"
    else:
        # Sub-cooled liquid pipe/equipment (P > Psat): single-phase (no liquid_level), 100%
        # liquid; runs dense until it flashes into the two-phase region. Static head off.
        D, vtype, horiz, thick = ID, "Flat-end", True, PIPE_WALL[(nominal, sched)]
        L = V / (math.pi * D ** 2 / 4)
        level, static, orient = None, False, "horizontal"
    vessel = {"length": round(L, 4), "diameter": round(D, 5), "thickness": thick,
              "heat_capacity": CP_W, "density": RHO_W, "orientation": orient, "type": vtype}
    if level is not None:
        vessel["liquid_level"] = round(level, 5)
    return {
        "vessel": vessel,
        "initial": {"temperature": round(T, 2), "pressure": round(P_bara * 1e5, 1), "fluid": "CO2"},
        "calculation": {"type": "energybalance", "time_step": dt, "end_time": end_time,
                        "non_equilibrium": True, "h_gas_liquid": "calc_two_sided",
                        "flow_work_fraction": 0.0},  # pure-u phase transfer: no spurious gas superheat
        "valve": {"flow": "discharge", "type": "none", "back_pressure": 101325.0},
        "release": {"type": rtype, "diameter": round(ID, 5), "discharge_coef": 1.0,
                    "liquid_nonequilibrium": 0.0, "static_head": static,
                    "back_pressure": 101325.0, "atm_pressure": 101325.0, "eos": "CoolProp",
                    "solid_in_vessel": True, "solid_h_gas_liquid": 0.0, "solid_h_gas_solid": 3.0,
                    "discharge_location": 0.0},
        "heat_transfer": {"type": "specified_h", "temp_ambient": round(T, 2), "h_outer": 0.0,
                          "h_inner": "churchill"},
    }


def run_one(d, dts):
    """Return (hd, dt, ok). On failure return the partially-populated hd (the near-depletion
    NEM fragility trips only after the useful source term is essentially complete)."""
    last_hd = None
    for dt in dts:
        dd = copy.deepcopy(d); dd["calculation"]["time_step"] = dt
        hd = HydDown(dd)
        try:
            hd.run(disable_pbar=True); return hd, dt, True
        except Exception:
            last_hd = hd
    return last_hd, dts[-1], False


def export(hd, path):
    """Write the valid source-term rows.

    Two cuts: (i) trailing zeros left by a mid-run failure, and (ii) the triple-point validity
    floor. Once the vessel reaches the triple point (~5.18 bar) dry ice forms inside it
    (solid-in-vessel regime, out of scope) and the release model stops computing the throat
    state, leaving the throat columns at 0 (which would read as 0 bar / -273 oC). We end the
    source term at the last row with a valid (non-zero) throat pressure."""
    hd.isrun = True
    df = hd.get_dataframe()
    p = np.asarray(hd.P)
    if np.any(p > 1e3):
        n = int(np.max(np.where(p > 1e3)[0])) + 1
        df = df.iloc[:n]
    if "Throat pressure (bar)" in df:
        good = df["Throat pressure (bar)"].values > 1e-9
        if good.any():
            df = df.iloc[:int(np.max(np.where(good)[0])) + 1]
    df.to_csv(path, index=False)
    return df


def main():
    filt = sys.argv[1] if len(sys.argv) > 1 else ""
    for (name, seg, V, P, T, nom, sch, kind, releases, dt, et) in CASES:
        if filt and filt not in name:
            continue
        for rtype in releases:
            d = build(name, seg, V, P, T, nom, sch, kind, rtype, dt, et)
            tag = "%s_%s" % (name, rtype)
            # tank: prefer moderate dt (large-inventory NEM is less stable at tiny dt)
            dts = (dt, dt / 2) if kind == "tank" else (dt, dt / 2, dt / 5)
            try:
                hd, used, ok = run_one(d, dts)
                if hd is None:
                    print("FAIL %-26s (no valid steps)" % tag); continue
                path = os.path.join(OUT, tag + ".csv")
                df = export(hd, path)
                t_end = df["Time (s)"].iloc[-1]
                rate = df["Release mass rate (kg/s)"].values
                print("%-4s %-26s dt=%.3f  t_end=%7.1fs  rows=%5d  released=%9.1f kg  peak=%8.1f kg/s"
                      % ("OK" if ok else "PART", tag, used, t_end, len(df),
                         float(np.trapz(rate, df["Time (s)"].values)), float(np.max(rate))))
            except Exception as e:
                print("FAIL %-26s %s" % (tag, repr(e)[:90]))


if __name__ == "__main__":
    main()
