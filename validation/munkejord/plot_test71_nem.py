"""Overlay the HydDown NEM Test 71 run vs measured 1 Hz data. Model starts at the two-phase
onset (~52 bar); its time is shifted so t=0 aligns with the measured 52 bar crossing."""
import os, csv
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__))
navy, red, amber, slate, grey = "#002D40", "#D61F39", "#E6A740", "#82979F", "#4C4D4E"

d = np.load(os.path.join(HERE, "test71_nem_out.npz"))
tmod = d["t"]; P = d["P"]/1e5; m = d["m"]; ms = d["ms"]
Tg = d["Tg"]-273.15; Tl = d["Tl"]-273.15
mdot = d["mdot"]

rows = list(csv.reader(open(os.path.join(HERE, "exp71_blowdown.csv"))))
h = rows[0]; a = np.array([[float(x) for x in r] for r in rows[1:]]); c = {n: i for i, n in enumerate(h)}
tM = a[:, c["t"]]; PM = a[:, c["PT163"]]; WM = a[:, c["Weight"]]; TbM = a[:, c["TT154"]]; TtM = a[:, c["TT114"]]
# align: measured time where P first crosses below 52 bar = model t=0
t0 = tM[np.argmax(PM < 52.0)]
tmod = tmod + t0
rM = np.full_like(tM, np.nan)
for i in range(len(tM)):
    k = (tM >= tM[i]-30) & (tM <= tM[i]+30)
    if k.sum() >= 10: rM[i] = -np.polyfit(tM[k], WM[k], 1)[0]

fig, ax = plt.subplots(2, 2, figsize=(13, 8.5), dpi=150)
ax[0, 0].plot(tM, PM, color=slate, lw=2.2, label="measured PT163")
ax[0, 0].plot(tmod, P, color=red, ls="--", lw=2, label="HydDown NEM")
ax[0, 0].axhline(5.18, color="k", ls=":", lw=0.8); ax[0, 0].set_ylabel("pressure (bar)")
ax[0, 0].set_title("Pressure (model starts at 52 bar flash onset)"); ax[0, 0].legend(fontsize=8)

ax[0, 1].plot(tM, WM, color=slate, lw=2.2, label="measured weight")
ax[0, 1].plot(tmod, m, color=red, ls="--", lw=2, label="model total")
ax[0, 1].plot(tmod, ms, color=navy, ls="-.", lw=1.8, label="model dry ice")
ax[0, 1].axhline(8.4, color=amber, ls=":", lw=1.4, label="measured retained 8.4 kg")
ax[0, 1].set_ylabel("mass (kg)"); ax[0, 1].set_title("Inventory + dry ice"); ax[0, 1].legend(fontsize=8)

# THE key panel: NEM gas/liquid split
ax[1, 0].plot(tM, TtM, color=grey, lw=1.6, label="measured top TT114")
ax[1, 0].plot(tM, TbM, color=slate, lw=2.0, label="measured bottom TT154")
ax[1, 0].plot(tmod, Tg, color=red, ls="--", lw=2, label="model GAS (NEM)")
ax[1, 0].plot(tmod, Tl, color=navy, ls="-.", lw=2, label="model LIQUID/solid (NEM)")
ax[1, 0].set_ylabel("temperature (C)"); ax[1, 0].set_title("Gas / liquid temperature split (NEM)"); ax[1, 0].legend(fontsize=8)

fin = np.isfinite(rM)
ax[1, 1].plot(tM[fin], rM[fin], color=slate, lw=1.0, alpha=0.8, label="measured (-dW/dt)")
ax[1, 1].plot(tmod, mdot, color=red, ls="--", lw=2, label="model mass rate")
ax[1, 1].set_ylabel("discharge rate (kg/s)"); ax[1, 1].set_title("Discharge rate")
ax[1, 1].set_ylim(-0.05, 0.7); ax[1, 1].legend(fontsize=8)

for a_ in ax.flat: a_.set_xlabel("time (s)"); a_.set_xlim(0, 240); a_.grid(alpha=0.3)
fig.suptitle("Munkejord Test 71 via HydDown NEM (CoolProp, nominal 222 kg wall, Cooper HTC) "
             "vs 1 Hz data", fontsize=11)
fig.tight_layout(rect=[0, 0, 1, 0.97])
out = os.path.join(HERE, "test71_nem_overlay.pdf")
fig.savefig(out); fig.savefig(out.replace(".pdf", ".png"))
print("wrote", out, "  t0(52bar)=", round(float(t0), 1), "s", flush=True); print("DONE", flush=True)
