import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.ticker import AutoMinorLocator
import matplotlib.cm as cm
import matplotlib.colors as colors
import matplotlib.patheffects as path_effects
from scipy.interpolate import griddata
import scipy.ndimage as ndi

rcParams.update({'font.size': 28})
# ----------------------------------------------------
# Choose which star to plot: "M2", "M0", or "M4"
# ----------------------------------------------------
STAR = "M0"   # <- change to "M0" or "M4" as needed
# ----------------------------------------------------
# Shared helpers
# ----------------------------------------------------
def m0_rebound_Z_transform(Z):
    """Z treatment for M0 rebound."""
    newZ = []
    for z in Z:
        if z == 0:
            newZ.append(10 ** 8.4)
        elif z < 10:
            newZ.append(10 ** z)
        else:
            newZ.append(z)
    return np.array(newZ)

def get_zi(filename, xmin, xmax, cols, z_transform=None,y_min=0.8, y_max=2.0, y_step=0.01, y_offset=0.0):
    """Load data and interpolate log10(Z) over a grid."""
    Y, X, Z = np.genfromtxt("../"+filename, delimiter=',', usecols=cols, unpack=True)
    if z_transform is not None:
        Y_proj, X_proj, Z_proj = np.genfromtxt("../data/Fig06/Tides_M0_proj.txt", delimiter=',', usecols=cols, unpack=True)
        for i in range(0,len(Z_proj)):
            row_idx = np.where(np.logical_and(np.abs(Y_proj[i]-Y)<1e-3,np.abs(X_proj[i]-X)<1e-3))[0]
            Z[row_idx] = 10**Z_proj[i]
        #Z = z_transform(Z)
    Z = np.log10(Z)
    xi = np.round(np.arange(xmin, xmax + 0.001, 0.001), 3)
    yi = np.round(np.arange(y_min, y_max + y_step, y_step) + y_offset, 2)
    zi = griddata((X, Y), Z, (xi[None, :], yi[:, None]), method='linear', fill_value=0.0)
    return xi, yi, zi

def get_clabel_locations(CS_sub, bounds):
    """Choose contour label locations within a bounding box."""
    xmin, xmax, ymin, ymax = bounds
    if not hasattr(CS_sub, 'levels') or not hasattr(CS_sub, 'allsegs'):
        return []
    manual_locations = []
    for segs_for_level in CS_sub.allsegs:
        if not segs_for_level:
            continue
        for path_vertices in segs_for_level:
            if len(path_vertices) < 2:
                continue
            mid_idx = len(path_vertices) // 2
            x, y = path_vertices[mid_idx]
            if xmin < x < xmax and ymin < y < ymax:
                manual_locations.append((x, y))
                break  # move to next level
    # Deduplicate, preserve order
    return list(dict.fromkeys(manual_locations))
# ----------------------------------------------------
# CONFIG for M2 and M0 (multi-panel plots)
# ----------------------------------------------------

MULTI_CONFIG = {
    "M2": {
        "xmin": 0.17, "xmax": 0.35, "label_bounds": (0.17, 0.35, 0.8, 2.0),
        "data_files": ["data/Fig02/Tides_M2.txt","data/Fig02/M2_sec_698.txt","data/Fig02/M2_sec_100.txt","data/Fig02/M2_sec_10.txt"],
        "reb_cols": (0, 1, -4),    # rebound file (Y, X, Z)
        "sec_cols": (0, 1, -1),    # secular files
        "reb_z_transform": None,   # no special transform
        "extra_reb_contour_level": None,"arrow": None, "sigma_smooth": 2.,
        "y_min": 0.8,"y_max": 2.0,"y_step": 0.01,"y_offset": 0.0,
        "xticks_main": np.arange(0.2, 0.4, 0.05),"xticks_bottom": None,
        "output_file": "Figs/Fig02.png","rebound_label_pos": (0.99, 0.95),"rebound_label_va": "top"},
    "M0": {
        "xmin": 0.27,"xmax": 0.52,"label_bounds": (0.27, 0.52, 0.8, 2.0),
        "data_files": ["data/Fig06/Tides_M0_long.txt","data/Fig06/M0_sec_698.txt","data/Fig06/M0_sec_100.txt","data/Fig06/M0_sec_10.txt"],
        "reb_cols": (0, 1, -4),     # rebound file (Y, X, Z_raw)
        "sec_cols": (0, 1, -1), "reb_z_transform": "m0_rebound", "extra_reb_contour_level": 8.30103, "sigma_smooth": 4.,
        "arrow": {"xytext": (0.488, 1.45), "xy": (0.497, 1.595)},
        "y_min": 0.8, "y_max": 2.0, "y_step": 0.01, "y_offset": 0.0,
        "xticks_main": np.arange(0.3, 0.53, 0.05), "xticks_bottom": np.arange(0.3, 0.53, 0.1),
        "output_file": "Figs/Fig06.png", "rebound_label_pos": (0.99, 0.9),"rebound_label_va": "bottom"},
}
def plot_panel(ax_idx, ax, fname, reb, cmap, vmin, vmax, sublbl, fig, cfg):
    """Multi-panel case: pcolormesh + contours (+ extras for rebound)."""
    xmin = cfg["xmin"]
    xmax = cfg["xmax"]
    # pick columns + optional Z transform depending on rebound or secular
    if reb:
        cols = cfg["reb_cols"]
        if cfg["reb_z_transform"] == "m0_rebound":
            z_transform = m0_rebound_Z_transform
        else:
            z_transform = None
    else:
        cols = cfg["sec_cols"]
        z_transform = None

    xi, yi, zi = get_zi(fname,xmin, xmax,cols,z_transform=z_transform,y_min=cfg["y_min"],y_max=cfg["y_max"],y_step=cfg["y_step"],y_offset=cfg["y_offset"])
    CS = ax.pcolormesh(xi, yi, zi, cmap=cmap, vmin=vmin, vmax=vmax, alpha=1.0)
    ax.text(0.03, 0.95, sublbl, transform=ax.transAxes,fontsize=20, fontweight='bold', va='top')
    ax.set_xticks(cfg["xticks_main"])
    yticks = np.round(np.arange(0.8, 2.2, 0.2), 1)
    ax.set_yticks(yticks)

    if reb:
        zi = ndi.gaussian_filter(zi, sigma=cfg["sigma_smooth"])
    # Optional extra white contour + arrow (M0 rebound case)
    if reb and cfg["extra_reb_contour_level"] is not None:
        ax.contour(xi, yi, zi, levels=[cfg["extra_reb_contour_level"]], colors='white', vmin=vmin, vmax=vmax)
        if cfg["arrow"] is not None:
            ax.annotate('',xy=cfg["arrow"]["xy"], xytext=cfg["arrow"]["xytext"], arrowprops=dict(arrowstyle='->', color='white', linewidth=2.5))
    CS2 = ax.contour(xi, yi, zi, levels=np.arange(6.0, 9.5, 0.5), colors='k', vmin=vmin, vmax=vmax)
    manual = get_clabel_locations(CS2, cfg["label_bounds"])
    manual_arg = manual if manual else None  # None => auto placement
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

# ----------------------------------------------------
# M4 config (single panel, simpler plot)
# ----------------------------------------------------
M4_CONFIG = {
    "filename": "data/Fig05/Tides_M4.txt","xmin": 0.08,"xmax": 0.18,"cols": (0, 1, -4),  # (Y, X, t_max)
    "y_min": 0.8,"y_max": 2.0,"y_step": 0.02,"y_offset": 0.01,"output_file": "Figs/Fig05.png"}

def plot_m4(cfg):
    """Single-panel lifetime map for M4, using the shared get_zi."""
    xi, yi, zi = get_zi( cfg["filename"], cfg["xmin"], cfg["xmax"], cfg["cols"],
        z_transform=None, y_min=cfg["y_min"], y_max=cfg["y_max"], y_step=cfg["y_step"], y_offset=cfg["y_offset"])
    
    fig = plt.figure(figsize=(7, 5), dpi=300)
    ax = fig.add_subplot(111)

    cmap2 = cm.cubehelix_r
    vmin, vmax = 3.6, 6.8 

    ax.pcolormesh(xi, yi, zi, cmap=cmap2, vmin=vmin, vmax=vmax)
    ax.set_xlabel(r"$a_\mathrm{p}$ (au)", size=24)
    ax.set_ylabel(r"$M_\mathrm{p}$ ($M_{\oplus}$)", size=24)
    ax.tick_params('both', labelsize=16, width=5)
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(which='minor', width=3)
    ax.set_ylim(cfg["y_min"], cfg["y_max"])
    norm = colors.Normalize(vmin=vmin, vmax=vmax)
    cax = fig.add_axes([0.92, 0.11, 0.015, 0.77])
    cbar = plt.colorbar(cm.ScalarMappable(cmap=cmap2, norm=norm), cax=cax, orientation='vertical')
    cbar.set_label(r"$\log_{10}(t_\mathrm{max})$", size=24)
    cbar.ax.tick_params(labelsize=16, width=5)
    cbar.ax.yaxis.set_minor_locator(AutoMinorLocator())
    cbar.ax.tick_params(which='minor', width=3)
    fig.savefig(cfg["output_file"], bbox_inches='tight', dpi=300)
    plt.close()

# ----------------------------------------------------
# Main: decide which figure to make
# ----------------------------------------------------

if STAR in MULTI_CONFIG:
    cfg = MULTI_CONFIG[STAR]

    fig = plt.figure(figsize=(10, 8), dpi=150)
    gs1 = fig.add_gridspec(nrows=2, ncols=3, wspace=0.1)

    ax1 = fig.add_subplot(gs1[0, :])
    ax2 = fig.add_subplot(gs1[1, 0])
    ax3 = fig.add_subplot(gs1[1, 1])
    ax4 = fig.add_subplot(gs1[1, 2])

    ax_list = [ax1, ax2, ax3, ax4]

    data_fn = cfg["data_files"]
    sub_lbl = ['a', 'b', 'c', 'd']
    rebound_flags = [True, False, False, False]

    vmin, vmax = 5.0, 9.5
    cmap_reb = cm.gist_rainbow
    cmap_sec = cm.gist_rainbow
    cmaps = [cmap_reb, cmap_sec, cmap_sec, cmap_sec]

    # Master loop over panels
    for i, ax in enumerate(ax_list):
        plot_panel(i, ax, data_fn[i], rebound_flags[i], cmaps[i], vmin, vmax, sub_lbl[i], fig, cfg)

    # Panel annotations
    ax4.text(0.99, 0.9, 'secular', transform=ax4.transAxes,fontsize=15, fontweight='bold', ha='right', color='k')
    ax1.text(cfg["rebound_label_pos"][0], cfg["rebound_label_pos"][1],'rebound', transform=ax1.transAxes,fontsize=15, fontweight='bold', ha='right',
             va=cfg["rebound_label_va"], color='k')

    # Shared labels
    fig.text(0.5, 0.03, r"$a_\mathrm{p}$ (au)", ha='center', size=24)
    fig.text(0.02, 0.5, r"$M_\mathrm{p}$ ($M_\oplus$)", va='center', rotation='vertical', size=24)
    fig.text(1.0, 0.5, r"$\log_{10}(t_\mathrm{max})$", va='center', rotation='vertical', size=24)

    # Remove y-tick labels on right-hand panels
    ax3.set_yticklabels([])
    ax4.set_yticklabels([])

    # Optional bottom-panel xtick override (M0)
    if cfg["xticks_bottom"] is not None:
        for ax in (ax2, ax3, ax4):
            ax.set_xticks(cfg["xticks_bottom"])

    fig.savefig(cfg["output_file"], bbox_inches='tight', dpi=300)
    plt.close()
elif STAR == "M4":
    plot_m4(M4_CONFIG)
else:
    raise ValueError("STAR must be 'M2', 'M0', or 'M4'")