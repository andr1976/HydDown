"""Below-triple wetted-wall descent for CARDICE gas tests: fixed 150 vs Rohsenow ('calc') vs
Cooper ('cooper'). 4th panel: the two correlations' HTC vs wall superheat (both capped 3000)."""
import copy, yaml, numpy as np, sys, math
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, r"C:/Users/AndersAndreasen/Documents/GitHub/HydDown/src")
from hyddown.hdclass import HydDown
import CoolProp as CP
from CoolProp.CoolProp import AbstractState, PropsSI
import hyddown.transport as tp, ht
from hyddown.co2_solid import P_TRIPLE
BASE = r"C:/Users/AndersAndreasen/Documents/GitHub/HydDown/"
navy, red, amber, slate = "#002D40", "#D61F39", "#E6A740", "#82979F"

def run(name, shi):
    d = copy.deepcopy(yaml.safe_load(open(BASE + f"validation/{name}.yml")))
    d["release"]["eos"] = "CoolProp"; d["release"]["solid_h_inner"] = shi
    hd = HydDown(d)
    try: hd.run(disable_pbar=True)
    except TypeError: hd.run()
    t = hd.time_array / 3600.0; ww = hd.T_vessel_wetted - 273.15
    m = hd.T_vessel_wetted > 100.0
    return t[m], ww[m]

TESTS = [("CARDICE_test5", "T5 gas 20 bar"), ("CARDICE_test7", "T7 gas 15 bar"),
         ("CARDICE_test9", "T9 gas 10 bar")]
fig, axes = plt.subplots(2, 2, figsize=(13, 9), dpi=150)
for ax, (name, title) in zip(axes.flat[:3], TESTS):
    for shi, c, ls, lab in [(150.0, slate, "-", "fixed 150"), ("calc", red, "--", "Rohsenow"),
                            ("cooper", navy, "-.", "Cooper")]:
        t, w = run(name, shi)
        ax.plot(t, w, color=c, ls=ls, lw=2.0, label=lab)
    ax.axhline(-75, color="k", ls=":", lw=1.0, label="measured ~ -75 C")
    ax.set_xlabel("time (h)"); ax.set_ylabel("wetted-wall T (C)")
    ax.set_title(title); ax.grid(alpha=0.3); ax.legend(fontsize=8)

# mechanism panel: Cooper vs Rohsenow HTC vs superheat at the triple point (both capped 3000)
ax = axes.flat[3]
liq = AbstractState("HEOS", "CO2"); wet = AbstractState("HEOS", "CO2")
liq.update(CP.PQ_INPUTS, P_TRIPLE, 0.0); wet.update(CP.PQ_INPUTS, P_TRIPLE, 0.0)
Pc = PropsSI("Pcrit", "CO2"); MW = PropsSI("molar_mass", "CO2") * 1000
Te = np.linspace(0.2, 16, 60)
hR = [tp.h_inside_wetted(1.58, 273.15 - 60 + te, 273.15 - 60, wet, liq) for te in Te]
hC = [min(ht.Cooper(P=P_TRIPLE, Pc=Pc, MW=MW, Te=te, Rp=1e-6), 3000.0) for te in Te]
ax.plot(Te, hR, color=red, ls="--", lw=2.2, label="Rohsenow (cap 3000)")
ax.plot(Te, hC, color=navy, ls="-.", lw=2.2, label="Cooper (cap 3000)")
ax.axhline(150, color=slate, ls=":", lw=1.4, label="fixed 150")
ax.set_xlabel("wall superheat  T$_w$-T$_{cold}$ (K)"); ax.set_ylabel("HTC (W/m$^2$K)")
ax.set_title("Cooper vs Rohsenow HTC vs superheat (at triple point)")
ax.grid(alpha=0.3); ax.legend(fontsize=8)

fig.suptitle("CARDICE gas tests: below-triple wetted wall - fixed 150 vs Rohsenow vs Cooper "
             "(identical dry ice; wall within 0.4 C)", fontsize=12)
fig.tight_layout(rect=[0, 0, 1, 0.97])
out = BASE + "validation/CARDICE_cooper_vs_rohsenow.pdf"
fig.savefig(out); fig.savefig(out.replace(".pdf", ".png"))
print("wrote", out, flush=True); print("DONE", flush=True)
