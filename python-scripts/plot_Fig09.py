import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rcParams

rcParams.update({"font.size": 14})

M_CEN   = 0.57
M_EARTH = 3.0034e-6
LABEL_FONTSIZE = 22
TICK_FONTSIZE  = "x-large"
xmin, xmax = 1e7,2.7e10

def secplot(mp, a_p_init, color, tau_label, legend_label, ax1, ax2, T_p):
    """Load secular data and plot both trajectories."""
    fname = f"m0_sec_{mp:1.1f}_{a_p_init:1.2f}{tau_label}.txt"
    data = np.genfromtxt("../data/Fig09/" + fname, delimiter=",")
    t,a_m,a_p = data[:,0], data[:,1], data[:,3]
    R_H = a_p * (mp * M_EARTH / (3 * M_CEN))**(1/3)
    T_m = np.sqrt(a_m**3 / (mp * M_EARTH))
    T_p_nc = np.sqrt(a_p**3 / M_CEN)
    x = t / T_p
    y = a_m / R_H
    # plot lines for a_m/R_H and period ratio
    ax1.plot(x,y,"--", color=color, label=legend_label)
    ax1.plot(x,y+0.035,"-",  color=color)
    ax2.plot(x,T_p_nc/T_m,"--", color=color, label=legend_label)
    return x, y
    
def style_axes(ax1, ax2):
    """Unified axis formatting."""
    ax1.set_yticks(np.arange(0.1,0.5,0.05))
    ax1.set_ylim(0.1, 0.44)
    ax1.set_ylabel(r"$a_\mathrm{m}\ (R_\mathrm{H})$", fontsize=LABEL_FONTSIZE)
    ax2.set_ylim(5, 30)
    ax2.set_ylabel(r"Period Ratio $P_\mathrm{p}/P_\mathrm{m}$", fontsize=LABEL_FONTSIZE)
    ax2.set_xlabel("Normalized Planetary Orbits", fontsize=LABEL_FONTSIZE)
    for ax in (ax1, ax2):
        ax.grid(True, color="gray", alpha=0.7)
        ax.set_xlim(xmin,xmax)
        ax.minorticks_on()
        ax.tick_params(which="major", length=8, width=4, labelsize=TICK_FONTSIZE, direction="out")
        ax.tick_params(which="minor", length=4, width=4, labelsize=TICK_FONTSIZE, direction="out")
        ax.set_xscale("log")
    ax1.legend(loc="lower right")
    ax1.set_xticklabels([])
    # Reference line + panel labels
    ax1.hlines(0.374, 1e7, 27e9, colors="blue", linestyles="dashdot")
    ax1.text(0.025, 0.97, "a", transform=ax1.transAxes, fontsize=20, weight="bold", va="top")
    ax2.text(0.025, 0.97, "b", transform=ax2.transAxes, fontsize=20, weight="bold", va="top")
# ===========================================================
# Main figure
# ===========================================================

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 8), dpi=150)
T_p = np.sqrt(0.27**3 / M_CEN)
configs = [("red","_5e9","698"),("purple", "_5e9_100tau",  "100")]
xs, ys = [], []
for color, tau_label, legend_label in configs:
    x, y = secplot(2.0, 0.52, color, tau_label, legend_label, ax1, ax2, T_p)
    xs.append(x)
    ys.append(y)
# Filled region between the two trajectories
ax1.fill_between(xs[0], ys[0] + 0.035, ys[1] + 0.035, alpha=0.5, color="gray")
style_axes(ax1, ax2)
# ===========================================================
# Top axis using twiny()
# ===========================================================
ax_top = ax1.twiny()
ax_top.set_xscale("log")
ax_top.set_xlim(xmin*T_p,xmax*T_p)
ax_top.set_xlabel("Time (yr)", fontsize=LABEL_FONTSIZE, labelpad=15)
ax_top.tick_params(which="major", length=8, width=4, labelsize=TICK_FONTSIZE, direction="out")
ax_top.tick_params(which="minor", length=4, width=4, labelsize=TICK_FONTSIZE, direction="out")

fig.subplots_adjust(hspace=0.1)
fig.savefig("../Figs/Fig09.png",bbox_inches="tight", dpi=300)
