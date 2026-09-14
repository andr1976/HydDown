"""Run every CARDICE validation with BOTH backends (thermopack tcPR vs CoolProp-only) and
compare key observables to confirm no regression from dropping thermopack."""
import sys, copy, yaml, numpy as np
sys.path.insert(0, r"C:/Users/AndersAndreasen/Documents/GitHub/HydDown/src")
from hyddown.hdclass import HydDown

YMLS = [
    "validation/CARDICE_test5.yml", "validation/CARDICE_test7.yml", "validation/CARDICE_test9.yml",
    "validation/CARDICE_batch6.yml", "validation/CARDICE_test8.yml", "validation/CARDICE_test10.yml",
    "validation/reconciled/CARDICE_test5.yml", "validation/reconciled/CARDICE_test6.yml",
    "validation/reconciled/CARDICE_test7.yml", "validation/reconciled/CARDICE_test8.yml",
    "validation/reconciled/CARDICE_test9.yml", "validation/reconciled/CARDICE_test10.yml",
]
BASE = r"C:/Users/AndersAndreasen/Documents/GitHub/HydDown/"

def metrics(hd):
    m = {"m0": float(hd.mass_fluid[0]), "retained_solid": float(hd.m_solid[-1]),
         "max_solid": float(np.max(hd.m_solid)), "P_end": float(hd.P[-1] / 1e5),
         "cum_dryice": float(hd.m_dryice_cum[-1]) if hasattr(hd, "m_dryice_cum") else 0.0}
    m["Tgas_min"] = float(np.min(hd.T_gas)) - 273.15 if hasattr(hd, "T_gas") else float("nan")
    m["Tliq_min"] = float(np.min(hd.T_liquid)) - 273.15 if hasattr(hd, "T_liquid") else float("nan")
    return m

def run(path, eos):
    d = yaml.safe_load(open(BASE + path))
    d = copy.deepcopy(d)
    d["release"]["eos"] = eos
    hd = HydDown(d)
    try:
        hd.run(disable_pbar=True)
    except TypeError:
        hd.run()
    return metrics(hd), type(hd.release_model).__name__

print(f"{'test':>34} {'backend':>16} | {'m0':>7} {'retDI':>7} {'maxDI':>7} {'Pend':>6} "
      f"{'Tg_min':>7} {'Tl_min':>7} {'cumDI':>7}", flush=True)
for y in YMLS:
    row = {}
    for eos in ("tcPR", "CoolProp"):
        try:
            m, backend = run(y, eos)
            row[eos] = m
            tag = y.replace("validation/", "").replace(".yml", "")
            print(f"{tag:>34} {backend:>16} | {m['m0']:7.1f} {m['retained_solid']:7.1f} "
                  f"{m['max_solid']:7.1f} {m['P_end']:6.2f} {m['Tgas_min']:7.1f} "
                  f"{m['Tliq_min']:7.1f} {m['cum_dryice']:7.1f}", flush=True)
        except Exception as e:
            print(f"{y:>34} {eos:>16} | ERROR: {str(e)[:70]}", flush=True)
    # regression delta on retained dry ice
    if "tcPR" in row and "CoolProp" in row:
        a, b = row["tcPR"]["retained_solid"], row["CoolProp"]["retained_solid"]
        dP = abs(row["tcPR"]["P_end"] - row["CoolProp"]["P_end"])
        print(f"{'  -> delta retained DI':>34} {'':>16} | tcPR {a:.1f} vs CP {b:.1f} kg "
              f"(d={b-a:+.1f} kg), dPend={dP:.2f} bar", flush=True)
print("DONE", flush=True)
