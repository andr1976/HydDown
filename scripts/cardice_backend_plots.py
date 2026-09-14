"""Overlay plots: CARDICE validations, thermopack tcPR (solid) vs CoolProp-only (dashed),
to visually confirm the thermopack-free backend matches across the full blowdown."""
import sys, copy, yaml, numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, r"C:/Users/AndersAndreasen/Documents/GitHub/HydDown/src")
from hyddown.hdclass import HydDown

BASE = r"C:/Users/AndersAndreasen/Documents/GitHub/HydDown/"
TESTS = [("CARDICE_test5", "T5 gas 20 bar"), ("CARDICE_test7", "T7 gas 15 bar"),
         ("CARDICE_test9", "T9 gas 10 bar"), ("CARDICE_batch6", "T6 liquid 20 bar"),
         ("CARDICE_test8", "T8 liquid 15 bar"), ("CARDICE_test10", "T10 liquid 10 bar")]
navy, red, amber, slate = "#002D40", "#D61F39", "#E6A740", "#82979F"

def run(name, eos):
    d = copy.deepcopy(yaml.safe_load(open(BASE + f"validation/{name}.yml")))
    d["release"]["eos"] = eos
    hd = HydDown(d)
    try:
        hd.run(disable_pbar=True)
    except TypeError:
        hd.run()
    t = hd.time_array / 3600.0  # hours
    Tg = hd.T_gas - 273.15 if hasattr(hd, "T_gas") else hd.T_fluid - 273.15
    return t, hd.P / 1e5, hd.m_solid.copy(), Tg

fig, axes = plt.subplots(2, 3, figsize=(16, 8.5), dpi=150)
for ax, (name, title) in zip(axes.flat, TESTS):
    t1, P1, S1, Tg1 = run(name, "tcPR")
    t2, P2, S2, Tg2 = run(name, "CoolProp")
    # pressure (left axis)
    ax.plot(t1, P1, color=navy, lw=2.2, label="P  tcPR")
    ax.plot(t2, P2, color=red, ls="--", lw=1.6, label="P  CoolProp")
    ax.set_ylabel("Pressure (bar)", color=navy); ax.tick_params(axis="y", labelcolor=navy)
    ax.set_xlabel("time (h)"); ax.set_title(title); ax.grid(alpha=0.3)
    # retained dry ice (right axis)
    axr = ax.twinx()
    axr.plot(t1, S1, color=slate, lw=2.2, label="dry ice  tcPR")
    axr.plot(t2, S2, color=amber, ls="--", lw=1.6, label="dry ice  CoolProp")
    axr.set_ylabel("in-vessel dry ice (kg)", color=slate); axr.tick_params(axis="y", labelcolor=slate)
    if ax is axes.flat[0]:
        ax.legend(loc="upper right", fontsize=7); axr.legend(loc="center right", fontsize=7)
fig.suptitle("CARDICE validations: thermopack tcPR (solid) vs CoolProp-only (dashed) — "
             "pressure + in-vessel dry ice", fontsize=13)
fig.tight_layout(rect=[0, 0, 1, 0.97])
out = BASE + "validation/CARDICE_backend_compare.pdf"
fig.savefig(out); fig.savefig(out.replace(".pdf", ".png"))
print("wrote", out, flush=True)
print("DONE", flush=True)
