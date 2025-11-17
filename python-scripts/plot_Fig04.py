import numpy as np
import matplotlib 
matplotlib.use('Agg')  # non-interactive backend
import matplotlib.pyplot as plt
from matplotlib import rcParams

rcParams.update({'font.size':14})

# constants
Mcen = 0.44
m_earth = 3.0034e-6
r_earth = 4.26352e-5
fs = 'x-large'

def style_axis(ax):
    ax.minorticks_on()
    ax.tick_params(which='major',direction='out',length=8,width=4,labelsize=fs)
    ax.tick_params(which='minor',direction='out',length=4,width=4,labelsize=fs)
    ax.grid(True,color='gray',alpha=0.7,lw=1.5)

def plot_spin(ax, mp, a_p, T_p, color, label):
    fname = f"../data/Fig04/reb_tide[{mp:1.2f},{a_p:1.3f}]_698.txt"
    data = np.genfromtxt(fname, delimiter=',', comments='#')
    time = data[:,0]
    spin_per = data[:,6]
    ax.plot(time/T_p, spin_per, '.', ms=5, color=color, label=label)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11,8), dpi=150)

T_p = np.sqrt(0.17**3 / Mcen)

# configs: (axis, planet mass, color)
configs = [(ax1, 1.0, 'deepskyblue'), (ax1, 1.0, 'red'),
           (ax2, 2.0, 'deepskyblue'), (ax2, 2.0, 'red')]

# corresponding semi-major axes
a_vals = [0.17, 0.34]*2
labels = ['0.17 au', '0.34 au']*2

for (ax, mp, color), a_p, lbl in zip(configs, a_vals, labels):
    plot_spin(ax, mp, a_p, T_p, color, lbl)

# panel formatting
panel_fmt = [(ax1, 'a', 'a', (1e6,2e8)),(ax2, 'b', 'b', (1e6,6e8))]
for ax, xtop_label, letter, xlim in panel_fmt:
    style_axis(ax)
    ax.set_xscale('log')
    ax.set_xlim(*xlim)
    ax.text(0.025,0.95,letter,transform=ax.transAxes,fontsize=26,fontweight='bold',va='top')
    ax.legend(loc='upper right',markerscale=5)

ax1.set_ylabel("Planetary Spin Period (hrs)", fontsize=28)

fig.supxlabel("Normalized Planetary Orbits", fontsize=24)

def xtoyears(x): return x*T_p
def xtoorbits(x): return x/T_p

for ax in (ax1, ax2):
    top = ax.secondary_xaxis('top', functions=(xtoyears, xtoorbits))
    style_axis(top)
    top.set_xscale('log')

fig.suptitle("Time (yr)", fontsize=22, y=0.98)
fig.savefig("../Figs/Fig04.png", bbox_inches='tight', dpi=150)

