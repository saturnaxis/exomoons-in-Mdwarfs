import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.ticker import AutoMinorLocator
import matplotlib.cm as cm
import matplotlib.patheffects as path_effects
from scipy.interpolate import griddata
import scipy.ndimage as ndi

rcParams.update({'font.size': 28})

def get_zi(filename, xmin, xmax, reb):
    """Load data from filename and interpolate log10(Z) over a grid."""
    cols = (0, 1, -4) if reb else (0, 1, -1)
    Y, X, Z = np.genfromtxt(filename, delimiter=',', usecols=cols, unpack=True)
    Z = np.log10(Z)
    xi = np.round(np.arange(xmin, xmax + 0.001, 0.001), 3)
    yi = np.round(np.arange(0.8, 2.01, 0.01), 3)
    zi = griddata((X, Y), Z, (xi[None, :], yi[:, None]), method='linear', fill_value=0.0)
    return xi, yi, zi
def get_clabel_locations(CS_sub):
    """Choose contour label locations within a bounding box."""
    xmin, xmax = 0.17, 0.35
    ymin, ymax = 0.8, 2.0
    if not hasattr(CS_sub, 'levels') or not hasattr(CS_sub, 'allsegs'):
        print("Warning: CS_sub missing 'levels' or 'allsegs'; no label locations.")
        return []
    manual_locations = []
    for level_index, segs_for_level in enumerate(CS_sub.allsegs):
        if not segs_for_level:
            continue
        for path_vertices in segs_for_level:
            if len(path_vertices) < 2:
                continue
            mid_idx = len(path_vertices) // 2
            x, y = path_vertices[mid_idx]
            if xmin < x < xmax and ymin < y < ymax:
                manual_locations.append((x, y))
                found_point = True
                break

    # Deduplicate locations
    unique_locations = list({loc: loc for loc in manual_locations}.values())
    return unique_locations

def plot_panel(ax_idx, ax, fname, reb, cmap, vmin, vmax, sublbl, fig):
    """Plot a single panel: background pcolormesh + contours + optional colorbar."""
    xi, yi, zi = get_zi(fname, 0.17, 0.35, reb)
    CS = ax.pcolormesh(xi, yi, zi, cmap=cmap, vmin=vmin, vmax=vmax, alpha=1.0)
    ax.text(0.03, 0.95, sublbl, transform=ax.transAxes, fontsize=20, fontweight='bold', va='top')
    ax.set_xticks(np.arange(0.2, 0.4, 0.05))
    yticks = np.round(np.arange(0.8, 2.2, 0.2), 1)
    ax.set_yticks(yticks)
    if reb:
        zi = ndi.gaussian_filter(zi, sigma=2.0)
    CS2 = ax.contour(xi, yi, zi, levels=np.arange(6.0, 9.5, 0.5), colors='k', vmin=vmin, vmax=vmax)
    manual = get_clabel_locations(CS2)
    manual_arg = manual if manual else None  # None => let matplotlib auto-place
    contour_labels = ax.clabel(CS2, inline=True, fontsize=12, manual=manual_arg)
    stroke_effect = [path_effects.withStroke(linewidth=1.5, foreground='lightgray')]
    for label in contour_labels:
        label.set_path_effects(stroke_effect)
    # Colorbar on the first panel only
    if ax_idx == 0:
        cax = fig.add_axes([0.93, 0.105, 0.015, 0.775])
        cbar = plt.colorbar(CS, cax=cax)
        cbar.ax.tick_params(labelsize=16, width=5)
        cbar.ax.yaxis.set_minor_locator(AutoMinorLocator())
        cbar.ax.tick_params(which='minor', width=3)
    ax.minorticks_on()
    ax.tick_params('both', labelsize=16, length=8, width=5)
    ax.tick_params(which='minor', length=4, width=3)


# --- Figure / axes ---
fig = plt.figure(figsize=(10, 8), dpi=150)
gs1 = fig.add_gridspec(nrows=2, ncols=3, wspace=0.1)

ax1 = fig.add_subplot(gs1[0, :])
ax2 = fig.add_subplot(gs1[1, 0])
ax3 = fig.add_subplot(gs1[1, 1])
ax4 = fig.add_subplot(gs1[1, 2])

ax_list = [ax1, ax2, ax3, ax4]
data_fn = ['Tides_M2.txt', 'Barnes/M2_sec_698.txt', 'Barnes/M2_sec_100.txt', 'Barnes/M2_sec_10.txt']
sub_lbl = ['a', 'b', 'c', 'd']
rebound_flags = [True, False, False, False]

vmin, vmax = 5.0, 9.5
cmap_reb = cm.gist_rainbow
cmap_sec = cm.gist_rainbow
cmaps = [cmap_reb, cmap_sec, cmap_sec, cmap_sec]

# --- Master loop over panels ---
for i, ax in enumerate(ax_list):
    plot_panel(i, ax, data_fn[i], rebound_flags[i], cmaps[i], vmin, vmax, sub_lbl[i], fig)

# Panel annotations
ax4.text(0.99, 0.9, 'secular', transform=ax4.transAxes, fontsize=15, fontweight='bold', ha='right', color='k')
ax1.text(0.99, 0.9, 'rebound', transform=ax1.transAxes, fontsize=15, fontweight='bold', ha='right', color='k')

# Shared labels
fig.text(0.5, 0.03, r"$a_\mathrm{p}$ (au)", ha='center', size=24)
fig.text(0.02, 0.5, r"$M_\mathrm{p}$ ($M_\oplus$)",va='center', rotation='vertical', size=24)
fig.text(1.0, 0.5, r"$\log_{10}(t_\mathrm{max})$", va='center', rotation='vertical', size=24)

# Remove y-tick labels on right-hand panels
ax3.set_yticklabels([])
ax4.set_yticklabels([])

fig.savefig('Figs/Fig02.png', bbox_inches='tight', dpi=300)
plt.close()
