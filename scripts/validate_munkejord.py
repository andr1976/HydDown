"""Run all 9 Munkejord dense-phase tests (CoolProp) and compare to the SINTEF 1 Hz data
(background/19589510.zip), writing one 2x2 figure per test to validation/munkejord/.

Temperatures are shown as a measured band (min-max over the 6 vertical positions) + midline,
for both fluid (TT1x4 centre) and inside wall (TT1x2). Each channel is read on ITS OWN time
column (fluid/wall temps -> tTT; pressures -> tPT16x; weight -> tWeight) to avoid the clock
mismatch. Vessel pressure = (PT162+PT163)/2 per the paper. Blowdown onset is taken from the
PT163 dense flash and each trace is shifted so onset = model t=0.

NOTE (branch co2-release-hem): the dry-ice-region wall currently tracks the sublimation line,
which is right for CARDICE's thick ice bed but over-cools Munkejord's thin-ice / thicker-relative
wall (see the parked phenomenological 2D lateral-conduction lever). So the model wetted-wall trace
here reads too cold in the gas/solid tail until that lever is added.

Run from the repo root:  python scripts/validate_munkejord.py
"""
import zipfile, io, os
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from hyddown.hdclass import HydDown

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(REPO, "validation", "munkejord"); os.makedirs(OUT, exist_ok=True)
NAVY = "#002D40"; RED = "#D61F39"; AMBER = "#E6A740"; SLATE = "#82979F"; GREEN = "#118a6b"

# name, P0_bar, T0_C, nozzle_mm, is_riser (liquid draw)
TESTS = [
    ("Exp71", 122.6, 25.2, 8.0, False), ("Exp72", 119.0, 24.9, 6.5, False), ("Exp75", 119.0, 25.0, 4.5, False),
    ("Exp52", 119.9, 15.4, 8.0, True), ("Exp53", 119.5, 24.4, 8.0, True), ("Exp56", 119.1, 15.2, 6.5, True),
    ("Exp57", 116.8, 24.5, 6.5, True), ("Exp45", 119.9, 14.5, 4.5, True), ("Exp46", 116.7, 24.4, 4.5, True),
]
FLUID = ["TT104", "TT114", "TT124", "TT134", "TT144", "TT154"]   # y=4 centre, 6 heights
IWALL = ["TT102", "TT112", "TT122", "TT132", "TT142", "TT152"]   # y=2 inside wall, 6 heights


# per-test end-time overrides [s] (zoom into the active blowdown for fast-emptying riser tests)
END_OVERRIDE = {"Exp53": 100.0, "Exp57": 150.0, "Exp46": 300.0}


def endtime(name, noz):
    if name in END_OVERRIDE:
        return END_OVERRIDE[name]
    return {8.0: 250.0, 6.5: 400.0, 4.5: 600.0}.get(noz, 300.0)


def measured(z, name, ET):
    xlsx = [n for n in z.namelist() if n.split("/")[-1].startswith(name + "_")][0]
    D = pd.read_excel(io.BytesIO(z.read(xlsx)), sheet_name="Data", header=None, engine="calamine")
    names = [str(x) for x in D.iloc[0].tolist()]; col = {n: i for i, n in enumerate(names)}
    dat = D.iloc[3:].reset_index(drop=True)

    def raw(tn, vn):
        t = pd.to_numeric(dat[col[tn]], errors="coerce").values
        v = pd.to_numeric(dat[col[vn]], errors="coerce").values
        m = np.isfinite(t) & np.isfinite(v); return t[m], v[m]
    # Vessel pressure: normally (PT162 + PT163)/2, but the bottom sensor PT163 is DEAD (reads 0)
    # in some riser tests (Exp52, Exp57). Average only the valid sensors - matching the paper,
    # which plots PT162 alone for those (Fig. A.9(a)). Blowdown onset (the dense flash within the
    # first sample) is taken from the valid vessel sensor: the early peak, then the first >5 bar drop.
    def valid(vn):
        _t, _v = raw("t" + vn, vn); return len(_v) > 0 and np.nanmax(_v) > 5.0
    v162, v163 = valid("PT162"), valid("PT163")
    onset_sensor = "PT163" if v163 else "PT162"
    tp, p = raw("t" + onset_sensor, onset_sensor); p0 = np.nanmax(p[:200])
    below = np.where(p < p0 - 5.0)[0]; t0 = tp[below[0]] if len(below) else 0.0
    tmax = min(tp.max(), ET + t0)
    g = np.arange(0.0, tmax - t0, 0.5)

    def I(tn, vn):
        t, v = raw(tn, vn); return np.interp(g, t - t0, v)
    if v162 and v163:
        Pavg = 0.5 * (I("tPT162", "PT162") + I("tPT163", "PT163"))
    else:
        Pavg = I("tPT162", "PT162") if v162 else I("tPT163", "PT163")
    W = I("tWeight", "Weight")
    fluid = np.vstack([I("tTT", c) for c in FLUID])
    iwall = np.vstack([I("tTT", c) for c in IWALL])
    return g, Pavg, W, fluid, iwall


def run_model(P0_bar, T0_C, nozzle_mm, riser, ET):
    T0 = T0_C + 273.15
    d = {
        "vessel": {"length": 1.0, "diameter": 0.273, "thickness": 0.0254, "heat_capacity": 500,
                   "density": 7950, "orientation": "vertical", "type": "Flat-end"},
        "initial": {"temperature": T0, "pressure": P0_bar * 1e5, "fluid": "CO2"},
        "calculation": {"type": "energybalance", "time_step": 0.1, "end_time": ET,
                        "non_equilibrium": True, "h_gas_liquid": "calc_two_sided"},
        "valve": {"flow": "discharge", "type": "none", "back_pressure": 101325.},
        "release": {"type": "liquid" if riser else "gas", "diameter": nozzle_mm / 1000., "discharge_coef": 1.0,
                    "liquid_nonequilibrium": 0.0, "back_pressure": 101325., "atm_pressure": 101325., "eos": "CoolProp",
                    "solid_in_vessel": True, "solid_h_inner": "cooper", "solid_h_gas_wall": "churchill",
                    "solid_h_gas_liquid": 0.0, "solid_h_gas_solid": 3.0,
                    "discharge_location": 0.009 if riser else 0.0},   # riser inlet 9 mm above floor
        "heat_transfer": {"type": "specified_h", "temp_ambient": T0, "h_outer": 0.0, "h_inner": "churchill"},
    }
    for dt in (0.1, 0.02):                 # dt=0.1 fast; fall back to 0.02 on a near-triple flash fail
        d["calculation"]["time_step"] = dt
        try:
            hd = HydDown(d); hd.run(disable_pbar=True); return hd
        except Exception:
            if dt == 0.02:
                raise
    return hd


def main():
    z = zipfile.ZipFile(os.path.join(REPO, "background", "19589510.zip"))
    for name, P0, T0, noz, riser in TESTS:
        try:
            ET = endtime(name, noz)
            g, Pavg, W, fluid, iwall = measured(z, name, ET)
            hd = run_model(P0, T0, noz, riser, ET)
            t = np.asarray(hd.time_array)
            kind = "riser / liquid draw" if riser else "no-riser / gas draw"
            fig, ax = plt.subplots(2, 2, figsize=(12, 8))
            fig.suptitle("Munkejord %s  (%.0f bar, %.1f C, %.1f mm, %s) - model (CoolProp) vs 1 Hz data"
                         % (name, P0, T0, noz, kind), color=NAVY, fontsize=12)
            a = ax[0, 0]; a.plot(g, Pavg, color=SLATE, lw=2, label="measured vessel P")
            a.plot(t, hd.P / 1e5, color=RED, lw=1.7, ls="--", label="model")
            a.set_ylabel("pressure [bar]"); a.set_title("Pressure"); a.legend(fontsize=8); a.grid(alpha=.3)
            a = ax[0, 1]; a.plot(g, W, color=SLATE, lw=2, label="measured weight")
            a.plot(t, hd.mass_fluid, color=RED, lw=1.7, ls="--", label="model total")
            a.plot(t, hd.m_liquid, color=NAVY, lw=1, ls="-.", label="liquid"); a.plot(t, hd.m_gas, color=AMBER, lw=1, label="gas")
            a.plot(t, hd.m_solid, color="k", lw=1, ls=":", label="dry ice")
            a.set_ylabel("mass [kg]"); a.set_title("Inventory"); a.legend(fontsize=7.5); a.grid(alpha=.3)
            a = ax[1, 0]
            a.fill_between(g, fluid.min(0), fluid.max(0), color=SLATE, alpha=0.3, label="measured band (6 heights)")
            a.plot(g, np.median(fluid, 0), color=SLATE, lw=1.4, label="measured midline")
            alive = np.asarray(hd.mass_fluid) > 0.02
            cond = np.asarray(hd.m_liquid) + np.asarray(hd.m_solid)
            a.plot(t, np.where(alive, np.asarray(hd.T_gas) - 273.15, np.nan), color=RED, lw=1.6, ls="--", label="model gas")
            a.plot(t, np.where(alive & (cond > 0.1), np.asarray(hd.T_liquid) - 273.15, np.nan), color=NAVY, lw=1.6, ls="-.", label="model liquid/solid")
            a.set_ylabel("fluid T [C]"); a.set_xlabel("time [s]"); a.set_title("Fluid temperature"); a.legend(fontsize=7.5); a.grid(alpha=.3)
            a = ax[1, 1]
            a.fill_between(g, iwall.min(0), iwall.max(0), color=GREEN, alpha=0.25, label="measured band (inside wall)")
            a.plot(g, np.median(iwall, 0), color=GREEN, lw=1.4, label="measured midline")
            a.plot(t, hd.T_vessel - 273.15, color=RED, lw=1.6, ls="--", label="model gas-wall")
            a.plot(t, hd.T_vessel_wetted - 273.15, color=NAVY, lw=1.6, ls="-.", label="model wetted wall")
            a.set_ylabel("inside-wall T [C]"); a.set_xlabel("time [s]"); a.set_title("Wall temperature"); a.legend(fontsize=7.5); a.grid(alpha=.3)
            for row in ax:
                for a in row:
                    a.set_xlim(0, ET)
            fig.tight_layout(rect=[0, 0, 1, 0.96])
            for e in ("pdf", "png"):
                fig.savefig(os.path.join(OUT, "munke_%s.%s" % (name, e)), dpi=150)
            plt.close(fig)
            print("%s: m0=%.1f final=%.1f dryice=%.1f kg | endP=%.2f bar  -> munke_%s.pdf"
                  % (name, hd.mass_fluid[0], hd.mass_fluid[-1], hd.m_solid[-1], hd.P[-1] / 1e5, name), flush=True)
        except Exception as ex:
            print("%s FAILED: %s" % (name, str(ex)[:120]), flush=True)
    print("ALL DONE", flush=True)


if __name__ == "__main__":
    main()
