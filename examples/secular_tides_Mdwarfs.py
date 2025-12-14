#!/usr/bin/env python3
"""
Tidal evolution integrator for planet–moon systems around M dwarfs.
Supports:
  - Single initial condition runs
  - Parallelized grid evaluation
Outputs: CSV lines with m_p, a_p, max e_m, final a_m/R_H, num orbits, lifetime.
"""

import numpy as np
from scipy.integrate import solve_ivp
import multiprocessing as mp
import os, sys

# ------------------------------------------------------------
# Physical constants
# ------------------------------------------------------------
G = 4.0 * np.pi**2              # AU^3 / (Msun * yr^2)
M_E = 3.0035e-6                 # Earth mass [Msun]
M_MOON = 3.6943e-8              # Lunar mass [Msun]
R_E = 4.2635e-5                 # Earth radius [AU]
C_P = 0.3308                    # Planet normalized MOI
KP = 0.298                      # Barnes 2017 k2p
PSI_E = 0.0                     # Obliquity (rad)
EP0 = 0.03                      # Initial planetary eccentricity
TSCALE = 5e9                    # Max integration time [yr]

# ------------------------------------------------------------
# Eccentricity-dependent f-functions from Hut (1981)
# ------------------------------------------------------------
def f_ecc(idx, e):
    if idx == 1: return 1 + 31/2*e**2 + 255/8*e**4 + 185/16*e**6 + 25/64*e**8
    if idx == 2: return 1 + 15/2*e**2 + 45/8*e**4 + 5/16*e**6
    if idx == 3: return 1 + 15/4*e**2 + 15/8*e**4 + 5/64*e**6
    if idx == 4: return 1 + 3/2*e**2 + 1/8*e**4
    if idx == 5: return 1 + 3*e**2 + 3/8*e**4

# ------------------------------------------------------------
# Tidal model: Barnes (2017)
# y = [a_m, e_m, a_p, e_p, Omega_p]
# ------------------------------------------------------------
def deriv_barnes(t, y, k2p, tau_p, M_p, R_p, M_m, M_star):
    am, em, ap, ep, Op = y
    beta_m = np.sqrt(1 - em**2)
    beta_p = np.sqrt(1 - ep**2)

    Z_pm = 3*G**2*k2p*M_m**2*(M_p+M_m)*R_p**5 / (am**9) * tau_p
    Z_ps = 3*G**2*k2p*M_star**2*(M_p+M_star)*R_p**5 / (ap**9) * tau_p

    n_m = np.sqrt(G*(M_p+M_m)/am**3)
    n_p = np.sqrt(G*(M_star+M_p)/ap**3)
    r_g = np.sqrt(C_P)

    dam = 2*am**2/(G*M_p*M_m) * Z_pm * (np.cos(PSI_E)*f_ecc(2,em)/beta_m**12 * Op/n_m - f_ecc(1,em)/beta_m**15)
    dem = 11*am*em/(2*G*M_p*M_m) * Z_pm * (np.cos(PSI_E)*f_ecc(4,em)/beta_m**10 * Op/n_m - f_ecc(3,em)/beta_m**13)

    dap = 2*ap**2/(G*M_star*M_p) * Z_ps * (np.cos(PSI_E)*f_ecc(2,ep)/beta_p**12 * Op/n_p - f_ecc(1,ep)/beta_p**15)
    dep = 11*ap*ep/(2*G*M_star*M_p) * Z_ps * (np.cos(PSI_E)*f_ecc(4,ep)/beta_p**10 * Op/n_p - f_ecc(4,ep)/beta_p**12)

    dOp = (Z_pm/(2*M_p*r_g**2*R_p**2*n_m)) * (2*np.cos(PSI_E)*f_ecc(2,em)/beta_m**12 - (1+np.cos(PSI_E)**2)*f_ecc(5,em)/beta_m**9 * Op/n_m)
    dOp += (Z_ps/(2*M_p*r_g**2*R_p**2*n_p)) * (2*np.cos(PSI_E)*f_ecc(2,ep)/beta_p**12 - (1+np.cos(PSI_E)**2)*f_ecc(5,ep)/beta_p**9 * Op/n_p)

    return [dam, dem, dap, dep, dOp]

# ------------------------------------------------------------
# Simulation wrapper
# ------------------------------------------------------------
def run_sim(m_mul, a_p, tau_sec, M_star, outpath):
    """Runs a single tidal evolution simulation."""
    M_p = m_mul * M_E
    M_m = M_MOON
    R_p = R_E
    tau_p = tau_sec / (3600*24*365.25)

    # Initial moon orbit
    rho_ratio = 5.515/3.34
    a_roche = 2.44 * R_p * rho_ratio**(1/3)
    a_m = 3*a_roche
    e_m = 1e-4

    n_m = np.sqrt(G*(M_p+M_m)/a_m**3)
    R_H = a_p * ((M_p+M_m)/(3*M_star))**(1/3)
    a_crit = 0.4031*(1 - 1.123*EP0)*R_H #Rosario-Franco et al. (2020)

    # Planet
    e_p = EP0
    Op = 2*np.pi / (5/(24*365.25))

    times = np.concatenate((np.arange(0, 1e4, 1e2),np.arange(1e4,1e6,1e4),np.arange(1e6,TSCALE+1e5,1e5))
    y = np.array([a_m, e_m, a_p, e_p, Op])
    e_max = 0.0

    for i in range(1, len(times)):
        sol = solve_ivp(deriv_barnes, [times[i-1], times[i]], y, method="RK45", atol=1e-12, rtol=1e-12, args=(KP, tau_p, M_p, R_p, M_m, M_star))
        y = sol.y[:, -1]
        am, em, ap, ep, Op = y

        e_max = max(e_max, em)

        # Escape / instability conditions
        if am*(1-em) < R_p or am*(1+em) > a_crit: break

    n_p = np.sqrt(G*(M_star+M_p+M_m)/a_p**3)
    lifetime = sol.t[-1]
    orbits = lifetime / (2*np.pi/n_p)

    with open(outpath, "a") as f:
        f.write(f"{m_mul:.2f}, {a_p:.3f}, {e_max:.3e}, {am/R_H:.3f}, {orbits:.3e}, {lifetime:.3e}\n")

# ------------------------------------------------------------
# Parallel worker
# ------------------------------------------------------------
def worker(args):
    return run_sim(*args)

# ------------------------------------------------------------
# Main execution
# ------------------------------------------------------------
def main():
    if len(sys.argv) < 4:
        print("Usage: secular_tides_Mdwarfs.py <start_idx> <tau_sec> <M0|M2|M4>")
        sys.exit(1)

    st_idx = int(sys.argv[1])
    tau_sec = int(sys.argv[2])
    st_type = sys.argv[3]

    # Stellar types
    if st_type == "M0": M_star, a_i, a_o = 0.57, 0.27, 0.521
    elif st_type == "M2": M_star, a_i, a_o = 0.44, 0.17, 0.351
    elif st_type == "M4": M_star, a_i, a_o = 0.23, 0.08, 0.181
    else:
        print("Error: st_type must be M0, M2, or M4")
        sys.exit(1)

    outpath = f"{st_type}_sec_{tau_sec}.txt"
    if not os.path.exists(outpath):
        with open(outpath, "w") as f:
            f.write("#M_p(M_E), a_p(AU), e_max, final a_m(R_H), N_orbits, lifetime(yr)\n")

    # Parameter grid
    params = [(round(m,3), round(a,3), tau_sec, M_star, outpath) for m in np.arange(0.8,2.01,0.01) for a in np.arange(a_i, a_o,0.001)]
    np.random.seed(42); np.random.shuffle(params)

    # Remove already completed runs
    done = set()
    data = np.genfromtxt(outpath, delimiter=",", comments="#")
    if data.size > 0:
        for row in np.atleast_2d(data):
            done.add((round(row[0],3), round(row[1],3)))
    params = [p for p in params if (p[0], p[1]) not in done]

    # Slice range
    params = params[st_idx: st_idx + 5000]

    # Parallel run
    with mp.Pool() as pool:
        pool.map(worker, params)

if __name__ == "__main__":
    main()
