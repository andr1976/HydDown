"""Plot the below-triple wetted-wall descent for CARDICE gas tests: hardcoded solid_h_inner=150
vs Rohsenow ('calc'). Rohsenow is active above the triple point in both cases. 4th panel shows
the Rohsenow HTC vs wall superheat (the mechanism)."""
import sys, copy, yaml, numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, r"C:/Users/AndersAndreasen/Documents/GitHub/HydDown/src")
from hyddown.hdclass import HydDown
import CoolProp as CP
from CoolProp.CoolProp import AbstractState
import hyddown.transport as tp
from hyddown.co2_solid import P_TRIPLE
BASE = r"C:/Users/AndersAndreasen/Documents/GitHub/HydDown/"
navy, red, amber, slate = "#002D40", "#D61F39", "#E6A740", "#82979F"

def run(name, shi):
    d = copy.deepcopy(yaml.safe_load(open(BASE + f"validation/{name}.yml")))
    d["release"]["eos"] = "CoolProp"; d["release"]["solid_h_inner"] = shi
    hd = HydDown(d)
    try:
        hd.run(disable_pbar=True)
    except TypeError:
        hd.run()
    t = hd.time_array / 3600.0
    ww = hd.T_vessel_wetted - 273.15
    mask = hd.T_vessel_wetted > 100.0  # below-triple wetted wall active
    return t[mask], ww[mask]

TESTS = [("CARDICE_test5", "T5 gas 20 bar"), ("CARDICE_test7", "T7 gas 15 bar"),
         ("CARDICE_test9", "T9 gas 10 bar")]
fig, axes = plt.subplots(2, 2, figsize=(13, 9), dpi=150)
for ax, (name, title) in zip(axes.flat[:3], TESTS):
    t1, w1 = run(name, 150.0)
    t2, w2 = run(name, "calc")
    ax.plot(t1, w1, color=navy, lw=2.4, label="hardcoded 150 W/m$^2$K")
    ax.plot(t2, w2, color=red, ls="--", lw=1.8, label="Rohsenow ('calc')")
    ax.axhline(-75, color=slate, ls=":", lw=1.2, label="measured wall ~ -75 C")
    ax.set_xlabel("time (h)"); ax.set_ylabel("wetted-wall temperature (C)")
    ax.set_title(title); ax.grid(alpha=0.3); ax.legend(fontsize=8)

# 4th panel: Rohsenow HTC vs wall superheat (mechanism)
ax = axes.flat[3]
liq = AbstractState("HEOS", "CO2"); wet = AbstractState("HEOS", "CO2")
liq.update(CP.PQ_INPUTS, P_TRIPLE, 0.0); wet.update(CP.PQ_INPUTS, P_TRIPLE, 0.0)
Te = np.linspace(0.2, 16, 60)
h = [tp.h_inside_wetted(1.58, 273.15 - 60 + te, 273.15 - 60, wet, liq) for te in Te]
ax.plot(Te, h, color=red, lw=2.4, label="Rohsenow h(superheat)")
ax.axhline(150, color=navy, ls="--", lw=1.8, label="hardcoded 150")
ax.axhline(3000, color=slate, ls=":", lw=1.0, label="cap 3000")
ax.set_xlabel("wall superheat  T$_w$ - T$_{cold}$ (K)")
ax.set_ylabel("boiling HTC (W/m$^2$K)")
ax.set_title("Rohsenow HTC vs superheat (why it matches 150)")
ax.grid(alpha=0.3); ax.legend(fontsize=8)

fig.suptitle("CARDICE gas tests: below-triple wetted wall — hardcoded 150 vs Rohsenow "
             "('calc'). Rohsenow active above triple in both.", fontsize=12)
fig.tight_layout(rect=[0, 0, 1, 0.97])
out = BASE + "validation/CARDICE_rohsenow_below.pdf"
fig.savefig(out); fig.savefig(out.replace(".pdf", ".png"))
print("wrote", out, flush=True); print("DONE", flush=True)
