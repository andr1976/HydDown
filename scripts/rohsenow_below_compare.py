"""CARDICE gas tests: below-triple wetted-wall HTC = hardcoded 150 vs Rohsenow ('calc').
Rohsenow above the triple point is active in BOTH cases (h_inner: calc). Report retained dry
ice and the wetted-wall temperature (measured bottom-of-vessel ~ -75 C)."""
import sys, copy, yaml, numpy as np
sys.path.insert(0, r"C:/Users/AndersAndreasen/Documents/GitHub/HydDown/src")
from hyddown.hdclass import HydDown
BASE = r"C:/Users/AndersAndreasen/Documents/GitHub/HydDown/"

def run(name, shi):
    d = copy.deepcopy(yaml.safe_load(open(BASE + f"validation/{name}.yml")))
    d["release"]["eos"] = "CoolProp"          # thermopack-free
    d["release"]["solid_h_inner"] = shi
    hd = HydDown(d)
    try:
        hd.run(disable_pbar=True)
    except TypeError:
        hd.run()
    ww = hd.T_vessel_wetted - 273.15
    # only the below-triple portion (wetted wall active) - ignore the 0-K parked pre-handoff
    ww = ww[hd.T_vessel_wetted > 100.0]
    return {"retDI": float(hd.m_solid[-1]), "ww_min": float(np.min(ww)),
            "ww_end": float(hd.T_vessel_wetted[-1] - 273.15),
            "gw_min": float(np.min(hd.T_vessel[hd.T_vessel > 100.0]) - 273.15)}

MEAS = {"CARDICE_test5": -75, "CARDICE_test7": -75, "CARDICE_test9": -75}  # measured coldest wall
print(f"{'test':>14} {'solid_h_inner':>14} | {'retDI(kg)':>9} {'wetwall_min':>11} "
      f"{'wetwall_end':>11} {'gaswall_min':>11}  (meas wall ~ -75 C)", flush=True)
for name in ("CARDICE_test5", "CARDICE_test7", "CARDICE_test9"):
    for shi in (150.0, "calc"):
        try:
            r = run(name, shi)
            print(f"{name:>14} {str(shi):>14} | {r['retDI']:9.1f} {r['ww_min']:11.1f} "
                  f"{r['ww_end']:11.1f} {r['gw_min']:11.1f}", flush=True)
        except Exception as e:
            print(f"{name:>14} {str(shi):>14} | ERROR: {str(e)[:60]}", flush=True)
print("DONE", flush=True)
