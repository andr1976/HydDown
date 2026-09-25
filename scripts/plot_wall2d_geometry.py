"""Geometry + finite-volume discretisation sketch of the 2-D conjugate wall (paper appendix).

Panel (a): the axisymmetric steel domain (bottom plate + shell + flange + lid), the structured
r-z FV grid, the adiabatic outer boundary and the inner Robin boundary split at the liquid level.
Panel (b): a zoom of the cylindrical shell showing the through-thickness cells and the inner
(fluid-side) / outer (ambient-side) faces. A representative (coarse) grid is drawn for legibility.

Writes paper/figures/wall2d_geometry.pdf (+ .png).  Run from the repo root.
"""
import os, sys
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrowPatch

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "src"))
from hyddown import wall2d as w2
OUT = os.path.join(REPO, "paper", "figures"); os.makedirs(OUT, exist_ok=True)
NAVY = "#002D40"; RED = "#D61F39"; AMBER = "#E6A740"; SLATE = "#82979F"; STEEL = "#c9d3d8"

plt.rcParams.update({"font.family": "serif", "font.size": 9, "mathtext.fontset": "cm"})

g = w2.default_sintef_geometry()
# coarse representative grid for the sketch
W = w2.WallConduction2D(g, 7950., 500., 298.15,
                        dr_wall=0.0085, dr_bulk=0.045, dz_wall=0.025, dz_bulk=0.075)
re, ze = W.re, W.ze
r_in, r_out, r_lid, H = g["r_in"], g["r_out"], g["r_lid"], g["H"]
z_flt = H + g["t_flange"]; z_top = z_flt + g["t_lid"]
Llev = 0.42                                            # illustrative liquid level

fig = plt.subplots(1, 2, figsize=(9.2, 6.6), gridspec_kw={"width_ratios": [2.0, 1.0]})
fig, (axA, axB) = fig

# ---------------- panel (a): domain + grid + BCs ----------------
for ir, iz in W.steel_cells:
    axA.add_patch(Rectangle((re[ir], ze[iz]), re[ir + 1] - re[ir], ze[iz + 1] - ze[iz],
                            facecolor=STEEL, edgecolor=SLATE, lw=0.4))
# cavity outline
axA.add_patch(Rectangle((0, 0), r_in, z_flt, facecolor="none", edgecolor=NAVY, lw=1.0, ls=":"))
# liquid level
axA.plot([0, r_in], [Llev, Llev], color=NAVY, lw=1.3, ls="--")
axA.text(r_in * 0.5, Llev + 0.02, "liquid / condensate level", ha="center", va="bottom",
         fontsize=7.5, color=NAVY)
# fluid zone labels
axA.text(r_in * 0.5, (z_flt + Llev) / 2, "vapour\nzone", ha="center", va="center", fontsize=8, color=SLATE)
axA.text(r_in * 0.5, Llev / 2, "liquid /\nsolid zone", ha="center", va="center", fontsize=8, color=SLATE)
# inner Robin BC annotations
axA.annotate("Churchill-Chu\nfree convection", xy=(r_in, (z_flt + Llev) / 2 + 0.1),
             xytext=(r_lid * 0.62, 0.86), fontsize=7.5, color=RED, ha="left",
             arrowprops=dict(arrowstyle="->", color=RED, lw=1.0))
axA.annotate("Cooper pool\nboiling", xy=(r_in, Llev / 2), xytext=(r_lid * 0.62, 0.30),
             fontsize=7.5, color=RED, ha="left",
             arrowprops=dict(arrowstyle="->", color=RED, lw=1.0))
# adiabatic outer BC
for zc in (-0.02, 0.5, z_top - 0.04):
    axA.annotate("", xy=(r_out + 0.006, 0.5), xytext=(r_out + 0.045, 0.5),
                 arrowprops=dict(arrowstyle="-|>", color=AMBER, lw=1.1))
axA.text(r_out + 0.048, 0.5, "adiabatic outer\n(insulated)", fontsize=7.5, color=AMBER,
         ha="left", va="center", rotation=0)
# dimension labels
def dim(x, z, s, c=NAVY, ha="center"):
    axA.text(x, z, s, fontsize=7, color=c, ha=ha, va="center")
dim(r_out * 0.5, -0.025, "bottom plate 50 mm")
dim(r_lid * 0.5, (z_flt + z_top) / 2, "lid 80 mm, dia. 580 mm")
dim(r_lid * 0.72, (H + z_flt) / 2, "flange 83 mm", ha="left")
axA.annotate("shell wall\n25.4 mm", xy=((r_in + r_out) / 2, 0.66), xytext=(r_lid * 0.55, 0.60),
             fontsize=7, color=NAVY, ha="left", arrowprops=dict(arrowstyle="->", color=SLATE, lw=0.8))
dim(r_in * 0.5, 0.06, "ID 273 mm", c=SLATE)
axA.set_xlabel("radial position $r$ [m]"); axA.set_ylabel("axial position $z$ [m]")
axA.set_title("(a) axisymmetric domain and FV grid", color=NAVY, fontsize=9)
axA.set_xlim(-0.01, r_lid + 0.13); axA.set_ylim(-0.075, z_top + 0.03); axA.set_aspect("equal")

# ---------------- panel (b): shell through-thickness zoom ----------------
# a short z-window of the shell
zwin = (0.45, 0.62)
for ir, iz in W.steel_cells:
    r0, z0 = re[ir], ze[iz]
    if r_in - 1e-6 <= W.rc[ir] <= r_out + 1e-6 and zwin[0] <= W.zc[iz] <= zwin[1]:
        axB.add_patch(Rectangle((r0, z0), re[ir + 1] - r0, ze[iz + 1] - z0,
                                facecolor=STEEL, edgecolor=SLATE, lw=0.7))
        axB.plot(W.rc[ir], W.zc[iz], "o", ms=2.5, color=NAVY)   # cell-centre nodes
axB.axvline(r_in, color=RED, lw=1.6)
axB.axvline(r_out, color=AMBER, lw=1.6)
axB.annotate("inner face\n(Robin: fluid)", xy=(r_in, 0.58), xytext=(r_in - 0.028, 0.60),
             fontsize=7.5, color=RED, ha="right", va="center",
             arrowprops=dict(arrowstyle="->", color=RED, lw=1.0))
axB.annotate("outer face\n(adiabatic)", xy=(r_out, 0.49), xytext=(r_out + 0.006, 0.47),
             fontsize=7.5, color=AMBER, ha="left", va="center",
             arrowprops=dict(arrowstyle="->", color=AMBER, lw=1.0))
axB.text((r_in + r_out) / 2, zwin[1] + 0.006, "through-thickness cells", ha="center",
         va="bottom", fontsize=7.5, color=NAVY)
axB.set_xlabel("$r$ [m]"); axB.set_ylabel("$z$ [m]")
axB.set_title("(b) shell wall (zoom)", color=NAVY, fontsize=9)
axB.set_xlim(r_in - 0.05, r_out + 0.05); axB.set_ylim(zwin[0] - 0.01, zwin[1] + 0.03)
axB.set_aspect("equal")

fig.tight_layout()
for e in ("pdf", "png"):
    fig.savefig(os.path.join(OUT, "wall2d_geometry.%s" % e), dpi=200, bbox_inches="tight")
plt.close(fig)
print("wrote paper/figures/wall2d_geometry.pdf/png  (grid %dx%d cells)" % (W.nr, W.nz))
