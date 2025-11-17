import numpy as np
import matplotlib.pyplot as plt

# --- Constants ---
T_SUN = 5780.0          # K, reference Teff for HZ polynomials
EARTH_MASS = 3.0035e-6  # Earth mass in solar masses
R_P_AU = 4.26352e-5     # Planet radius in AU for Hill-radius scaling

def hz_bound(T_eff, L, border):
    """Inner HZ distance (au) for Teff [K], L [L_sun]."""
    T = T_eff - T_SUN
    if border =='inner':
        S = (1.107 + 1.332e-4*T + 1.58e-8*T**2 - 8.308e-12*T**3 - 1.931e-15*T**4)
    elif border =='outer':
        S = (0.356 + 6.171e-5*T + 1.698e-9*T**2 - 3.198e-12*T**3 - 5.575e-16*T**4)
    return np.sqrt(L / S)

def hill_radius(a, m_star, m_planet=EARTH_MASS):
    """Hill radius (au) at semi-major axis a [au]."""
    return a * np.cbrt(m_planet / (3.0 * m_star))

# --- Stellar parameters for M0, M2, M4 ---
types  = np.array(['M0', 'M2', 'M4'])
Teff   = np.array([3850.0, 3560.0, 3210.0])
L_star = np.array([0.069, 0.029, 7.2e-3])
M_star = np.array([0.57, 0.44, 0.23])

colors = ['red', 'blue', 'black']
alphas = [1.0, 0.5, 1.0]
xpos   = np.arange(1, len(types) + 1)

# --- HZ edges and midpoints ---
hz_in  = hz_bound(Teff, L_star,'inner')
hz_out = hz_bound(Teff, L_star,'outer')
hz_mid = 0.5 * (hz_in + hz_out)

hz_err_lower = hz_mid - hz_in
hz_err_upper = hz_out - hz_mid

# --- Hill radii (scaled by planet radius) ---
rh_in  = hill_radius(hz_in,  M_star) / R_P_AU
rh_out = hill_radius(hz_out, M_star) / R_P_AU
rh_mid = hill_radius(hz_mid, M_star) / R_P_AU

rh_err_lower = rh_mid - rh_in
rh_err_upper = rh_out - rh_mid

# --- Figure / axes ---
fig = plt.figure(constrained_layout=True)
gs = fig.add_gridspec(nrows=2, ncols=3, wspace=0.15)
ax1 = fig.add_subplot(gs[0, 2])     # HZ vs spectral type
ax2 = fig.add_subplot(gs[1, 2])     # Hill radius vs spectral type
ax3 = fig.add_subplot(gs[:, :2])    # Bullseye HZ plot

axes  = [ax1, ax2, ax3]
names = ['hz', 'hill', 'bull']

# Shared angle grid for bullseye annuli
theta = np.linspace(0, 2 * np.pi, 50, endpoint=True)

# --- One master loop over axes ---
for ax, name in zip(axes, names):
    if name == 'hz':
        # HZ error bars
        for i in range(len(types)):
            ax.errorbar([xpos[i]], [hz_mid[i]], yerr=[[hz_err_lower[i]], [hz_err_upper[i]]],fmt='none', capsize=6, color=colors[i], alpha=alphas[i])
        # Overlap rectangle between M0 and M2 HZ
        m0_idx, m2_idx = 0, 1
        overlap_height = hz_out[m2_idx] - hz_in[m0_idx]
        for color, alpha in [('red', 1.0), ('blue', 0.5)]:
            ax.add_patch(plt.Rectangle((1, hz_in[m0_idx]), 1, overlap_height,color=color, alpha=alpha, fill=True))
        ax.set_ylabel("HZ (au)", fontsize=20)
        ax.set_xticks(xpos)
        ax.set_xticklabels([])
        ax.set_yticks([0.1, 0.3, 0.5])
        ax.text(0.8, 0.95, 'b', transform=ax.transAxes,fontsize=20, fontweight='bold', va='top')
    elif name == 'hill':
        # Hill-radius error bars
        for i in range(len(types)):
            ax.errorbar([xpos[i]], [rh_mid[i]], yerr=[[rh_err_lower[i]], [rh_err_upper[i]]],fmt='o', capsize=6, color=colors[i], alpha=alphas[i])
        ax.set_ylabel(r"$R_\mathrm{H}$ ($R_\mathrm{p}$)", fontsize=20)
        ax.set_xticks(xpos)
        ax.set_xticklabels(types)
        ax.text(0.8, 0.95, 'c', transform=ax.transAxes,fontsize=20, fontweight='bold', va='top')
    elif name == 'bull':
        # Bullseye HZ annuli
        for i in range(len(types)):
            radii = np.array([hz_in[i], hz_out[i]])
            x = np.outer(radii, np.cos(theta))
            y = np.outer(radii, np.sin(theta))
            x[1, :] = x[1, ::-1]
            y[1, :] = y[1, ::-1]
            ax.fill(x.ravel(), y.ravel(), color=colors[i], alpha=alphas[i], label=types[i], edgecolor=None)

        ax.plot(0, 0, '*', color='red', markersize=5)
        ax.legend(ncols=3, bbox_to_anchor=(1.02, 1.16))
        ax.set_xlim(-0.6, 0.6)
        ax.set_ylim(-0.6, 0.6)
        ax.set_xticks([-0.5, 0.0, 0.5])
        ax.set_yticks([-0.5, 0.0, 0.5])
        ax.set_aspect('equal')
        ax.set_xlabel("X (au)", fontsize=20)
        ax.set_ylabel("Y (au)", fontsize=20)
        ax.text(0.9, 0.95, 'a', transform=ax.transAxes, fontsize=20, fontweight='bold', va='top')

    # Shared formatting for all three axes
    ax.minorticks_on()
    ax.tick_params(which='major', direction='out', length=8.0, width=4.0, labelsize=18)
    ax.tick_params(which='minor', direction='out', length=4.0, width=2.0, labelsize=18)

    # Only the first two panels get x-limits in this style
    if name in ['hz', 'hill']:
        ax.set_xlim(0.7, 3.3)

fig.savefig('../Figs/Fig01.pdf', bbox_inches='tight', dpi=300)