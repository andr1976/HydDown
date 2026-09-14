"""Full Test 71 overlay (CoolProp + Rohsenow, NOMINAL 213 kg thermal mass) vs measured 1 Hz data:
pressure, inventory + dry ice, fluid temperature, discharge rate."""
import os, sys, subprocess, csv
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__))
navy, red, amber, slate, grey = "#002D40", "#D61F39", "#E6A740", "#82979F", "#4C4D4E"

# run the trial with the nominal thermal mass + Rohsenow (saves trial71_out.npz)
env = dict(os.environ, EOS="CP", MSTEEL="213", ROHSENOW="1")
subprocess.run([sys.executable, os.path.join(HERE, "trial_test71.py")], env=env, check=True)
d = np.load(os.path.join(HERE, "trial71_out.npz"))
tm_, P_, Tf_, m_, ms_ = d["t"], d["P"], d["Tf"], d["m"], d["ms"]
mdot_model = np.gradient(-m_, tm_)  # kg/s

# measured
rows = list(csv.reader(open(os.path.join(HERE, "exp71_blowdown.csv"))))
h = rows[0]; data = np.array([[float(x) for x in r] for r in rows[1:]]); c = {n: i for i, n in enumerate(h)}
tM = data[:, c["t"]]; PM = data[:, c["PT163"]]; WM = data[:, c["Weight"]]
TbM = data[:, c["TT154"]]; TtM = data[:, c["TT114"]]
# measured discharge rate: -dW/dt (61s rolling slope)
rM = np.full_like(tM, np.nan)
for i in range(len(tM)):
    k = (tM >= tM[i] - 30) & (tM <= tM[i] + 30)
    if k.sum() >= 10:
        rM[i] = -np.polyfit(tM[k], WM[k], 1)[0]

fig, ax = plt.subplots(2, 2, figsize=(13, 8.5), dpi=150)
# 1 pressure
ax[0, 0].plot(tM, PM, color=slate, lw=2.2, label="measured PT163")
ax[0, 0].plot(tm_, P_, color=red, ls="--", lw=2, label="model")
ax[0, 0].axhline(5.18, color="k", ls=":", lw=0.8, label="triple 5.18 bar")
ax[0, 0].set_ylabel("pressure (bar)"); ax[0, 0].set_title("Pressure"); ax[0, 0].legend(fontsize=8)
# 2 inventory + dry ice
ax[0, 1].plot(tM, WM, color=slate, lw=2.2, label="measured weight")
ax[0, 1].plot(tm_, m_, color=red, ls="--", lw=2, label="model total")
ax[0, 1].plot(tm_, ms_, color=navy, ls="-.", lw=1.8, label="model dry ice")
ax[0, 1].axhline(8.4, color=amber, ls=":", lw=1.4, label="measured retained 8.4 kg")
ax[0, 1].set_ylabel("mass (kg)"); ax[0, 1].set_title("Inventory + dry ice"); ax[0, 1].legend(fontsize=8)
# 3 fluid temperature
ax[1, 0].plot(tM, TbM, color=slate, lw=2.2, label="measured bottom TT154")
ax[1, 0].plot(tM, TtM, color=grey, lw=1.3, label="measured top TT114")
ax[1, 0].plot(tm_, Tf_, color=red, ls="--", lw=2, label="model fluid")
ax[1, 0].set_ylabel("fluid T (C)"); ax[1, 0].set_title("Fluid temperature"); ax[1, 0].legend(fontsize=8)
# 4 discharge rate
fin = np.isfinite(rM)
ax[1, 1].plot(tM[fin], rM[fin], color=slate, lw=1.0, alpha=0.8, label="measured (-dW/dt)")
ax[1, 1].plot(tm_, mdot_model, color=red, ls="--", lw=2, label="model (-dm/dt)")
ax[1, 1].set_ylabel("discharge rate (kg/s)"); ax[1, 1].set_title("Discharge rate")
ax[1, 1].set_ylim(-0.1, np.nanpercentile(rM[fin], 98) * 1.4 if fin.any() else 3); ax[1, 1].legend(fontsize=8)
for a in ax.flat:
    a.set_xlabel("time (s)"); a.set_xlim(0, 260); a.grid(alpha=0.3)
fig.suptitle("Munkejord Test 71 (8 mm, no riser, 122.6 bar / 25.2 C) - HydDown CoolProp + Rohsenow, "
             "NOMINAL 213 kg wall vs 1 Hz data", fontsize=11)
fig.tight_layout(rect=[0, 0, 1, 0.97])
out = os.path.join(HERE, "test71_nominal_full.pdf")
fig.savefig(out); fig.savefig(out.replace(".pdf", ".png"))
print("wrote", out, flush=True)
print(f"retained dry ice model {ms_[-1]:.1f} kg vs measured 8.4 kg; min Tf {Tf_.min():.1f} C", flush=True)
print("DONE", flush=True)
