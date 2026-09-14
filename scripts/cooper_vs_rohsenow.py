"""Cooper vs Rohsenow below-triple wetted-wall HTC across all 6 CARDICE tests (CoolProp backend).
Reports retained dry ice + wetted-wall minimum temperature (measured coldest wall ~ -75 C for the
gas tests; liquid tests empty before the plateau so the wetted wall barely activates)."""
import copy, yaml, numpy as np, sys
sys.path.insert(0, r"C:/Users/AndersAndreasen/Documents/GitHub/HydDown/src")
from hyddown.hdclass import HydDown
BASE = r"C:/Users/AndersAndreasen/Documents/GitHub/HydDown/"
TESTS = [("CARDICE_test5", "T5 gas 20bar"), ("CARDICE_test7", "T7 gas 15bar"),
         ("CARDICE_test9", "T9 gas 10bar"), ("CARDICE_batch6", "T6 liq 20bar"),
         ("CARDICE_test8", "T8 liq 15bar"), ("CARDICE_test10", "T10 liq 10bar")]

def run(name, shi):
    d = copy.deepcopy(yaml.safe_load(open(BASE + f"validation/{name}.yml")))
    d["release"]["eos"] = "CoolProp"; d["release"]["solid_h_inner"] = shi
    hd = HydDown(d)
    try:
        hd.run(disable_pbar=True)
    except TypeError:
        hd.run()
    ww = hd.T_vessel_wetted - 273.15
    ww = ww[hd.T_vessel_wetted > 100.0]
    ww_min = float(np.min(ww)) if ww.size else float("nan")
    return float(hd.m_solid[-1]), ww_min

print(f"{'test':>14} | {'retDI Rohsenow':>14} {'retDI Cooper':>12} | "
      f"{'wall Rohsenow':>13} {'wall Cooper':>11} {'dWall':>6}", flush=True)
for name, title in TESTS:
    dr, wr = run(name, "calc")     # Rohsenow
    dc, wc = run(name, "cooper")   # Cooper
    print(f"{title:>14} | {dr:14.1f} {dc:12.1f} | {wr:13.1f} {wc:11.1f} {wc-wr:6.2f}", flush=True)
print("DONE", flush=True)
