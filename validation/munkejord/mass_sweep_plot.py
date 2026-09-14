"""Test 71 retained dry ice vs lumped wall mass (Rohsenow HTC), showing that the NOMINAL
thermal mass under-predicts dry ice in a lumped wall (over-coupling) - the wall needs
transient conduction to use the nominal mass correctly."""
import os, subprocess, sys, re
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__))
navy, red, amber, slate = "#002D40", "#D61F39", "#E6A740", "#82979F"

masses = [140, 160, 175, 189, 213, 250, 300, 390]
ret = []
for m in masses:
    env = dict(os.environ, EOS="CP", MSTEEL=str(m), ROHSENOW="1")
    out = subprocess.run([sys.executable, os.path.join(HERE, "trial_test71.py")],
                         capture_output=True, text=True, env=env).stdout
    r = re.search(r"retained dry ice=([0-9.]+)", out)
    ret.append(float(r.group(1)) if r else np.nan)
    print(f"m_steel={m} kg -> retained dry ice {ret[-1]:.2f} kg", flush=True)

masses = np.array(masses); ret = np.array(ret)
fig, ax = plt.subplots(figsize=(9, 5.5), dpi=150)
ax.plot(masses, ret, "o-", color=navy, lw=2, label="lumped wall + Rohsenow")
ax.axhline(8.4, color=red, ls="--", lw=1.8, label="measured (8.4 kg)")
ax.axhline(8.6, color=amber, ls=":", lw=1.5, label="paper GERG+FEM (8.6 kg)")
# markers for the physically-meaningful masses
for m, lab, dy in [(189, "side wall\n189", 8), (213, "nominal shell\n213", -18),
                   (390, "full vessel\n390", 8), (175, "earlier tuned\n175", -20)]:
    y = float(np.interp(m, masses, ret))
    ax.annotate(lab, (m, y), textcoords="offset points", xytext=(0, dy),
                ha="center", fontsize=8, color=slate)
    ax.plot(m, y, "s", color=slate, ms=6)
ax.set_xlabel("lumped wall thermal mass m$_{steel}$ (kg)")
ax.set_ylabel("retained dry ice (kg)")
ax.set_title("Munkejord Test 71: retained dry ice vs lumped wall mass (Rohsenow HTC)\n"
             "nominal mass under-predicts (lumped over-couples); conduction gradient needed",
             fontsize=10)
ax.grid(alpha=0.3); ax.legend()
fig.tight_layout()
out = os.path.join(HERE, "test71_mass_sweep.pdf")
fig.savefig(out); fig.savefig(out.replace(".pdf", ".png"))
print("wrote", out, flush=True); print("DONE", flush=True)
