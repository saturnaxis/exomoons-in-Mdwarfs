"""
Tidal evolution & spin–orbit dynamics using REBOUND/REBOUNDx.

Star types:
  M0 = {M_star=0.57 Msun, a_i=0.27, a_o=0.521}
  M2 = {M_star=0.44 Msun, a_i=0.17, a_o=0.351}
  M4 = {M_star=0.23 Msun, a_i=0.08, a_o=0.181}

Usage
-----
Grid mode (parallel):
    python tide_spin_moon_sim.py STAR_TYPE START_INDEX
Example:
    python tide_spin_moon_sim.py M0 0

Single initial condition:
    python tide_spin_moon_sim.py STAR_TYPE single M_EARTH A_AU
Example:
    python tide_spin_moon_sim.py M2 single 1.20 0.250
"""

from __future__ import annotations
import multiprocessing as mp
import os, sys, threading, time
from typing import List, Tuple
import numpy as np
import rebound, reboundx
from scipy.constants import year, au

# ---------- Global configuration ----------
LOCK = threading.Lock()
N_CORES_TARGET = 60
G = 4 * np.pi**2
M_SUN = 1.989e30
KM_IN_AU = 1.496978707e8
M_EARTH = 3e-6      # Msun
R_EARTH = 4.259e-5  # AU
K2_EARTH = 0.298
C_EARTH = 0.3308
RHO_EARTH = 5.515 / (M_SUN * 1e3) * (KM_IN_AU * 1e5) ** 3
M_MOON = 3.6943e-8
TAU = 698                  # s
INTEGRATOR = "trace"
R_P = 1.0                  # Earth radii
ECC_P = 0.03
INC_P = 0.0
BETA_P = np.radians(23.76)
R_PL = R_P * R_EARTH       # AU
P_ROT_HOURS = 5.0
D_ROCHE = 2.44 * R_PL * (5.515 / 3.34) ** (1.0 / 3.0)
MAX_ORBITS = 2e8
MIN_SNAPSHOT_MB = 0.5
M_MIN, M_MAX, M_STEP = 0.8, 2.0, 0.01
A_STEP = 0.001
BASE_OUTPUT_DIR = os.environ.get("TIDE_SIM_HOME", os.getcwd())
STAR_CONFIGS = {"M0": {"mass": 0.57, "a_i": 0.27, "a_o": 0.521}, "M2": {"mass": 0.44, "a_i": 0.17, "a_o": 0.351}, "M4": {"mass": 0.23, "a_i": 0.08, "a_o": 0.181}}
STAR_MASS = None
OUTDIR = None
PARAMS: List[Tuple[float, float]] = []

# ---------- I/O helpers ----------
def write_to_file(filename: str, text: str, mode: str) -> None:
    """Thread-safe append/write; mode='head' => overwrite, else append."""
    with LOCK:
        file_mode = "w" if mode == "head" else "a"
        with open(filename, file_mode) as f: f.write(text)

# ---------- Orbital & tidal utilities ----------
def cartesian_to_orbital(ps0, ps1):
    """Return (a, P, e, inc, Omega) for ps1 relative to ps0."""
    m_cen, m_p = ps0.m, ps1.m
    delx, dely, delz = ps1.x - ps0.x, ps1.y - ps0.y, ps1.z - ps0.z
    delvx, delvy, delvz = ps1.vx - ps0.vx, ps1.vy - ps0.vy, ps1.vz - ps0.vz
    temp_sim = rebound.Simulation()
    temp_sim.units = ("AU", "Msun", "years")
    temp_sim.add(m=m_cen)
    temp_sim.add(m=m_p, x=delx, y=dely, z=delz, vx=delvx, vy=delvy, vz=delvz)
    temp_sim.move_to_com()
    ps = temp_sim.particles
    return ps[1].a, ps[1].P, ps[1].e, ps[1].inc, ps[1].Omega

def eccentricity_function(f_idx: int, ecc: float) -> float:
    e2 = ecc**2
    if f_idx == 1: return 1.0 + 31.0 / 2.0 * e2 + 255.0 / 8.0 * e2**2 + 185.0 / 16.0 * e2**3 + 25.0 / 64.0 * e2**4
    if f_idx == 2: return 1.0 + 15.0 / 2.0 * e2 + 45.0 / 8.0 * e2**2 + 5.0 / 16.0 * e2**3
    if f_idx == 3: return 1.0 + 15.0 / 4.0 * e2 + 15.0 / 8.0 * e2**2 + 5.0 / 64.0 * e2**3
    if f_idx == 4: return 1.0 + 3.0 / 2.0 * e2 + 1.0 / 8.0 * e2**2
    if f_idx == 5: return 1.0 + 3.0 * e2 + 3.0 / 8.0 * e2**2
    raise ValueError(f"Unknown eccentricity function index f_idx={f_idx}")

def calculate_alpha(prot_hours: float, m_star: float, m_pl: float, a_pl: float) -> float:
    """Precession constant alpha (arcsec/yr) from Lissauer+ (2012) Appendix A."""
    w = 2 * np.pi / (prot_hours / 24.0 / 365.25)
    D = 25.0 * (3.0 * C_EARTH / 2.0 - 1.0) ** 2 / 4.0 + 1.0
    x = 30.0 * w**2 / (np.pi * D * G * RHO_EARTH)
    alpha = 0.0
    if x <= 4.0:
        root = np.sqrt(4.0 - x)
        req = (G * m_pl * D * (2.0 - root) / (10.0 * w**2)) ** (1.0 / 3.0)
        _ = req
        j2 = (2.0 - root) * (5.0 - D) / 30.0
        alpha = 1.5 * (G * (m_star + m_pl) / a_pl**3) / w * (j2 / C_EARTH) * (180.0 / np.pi) * 3600.0
    return alpha

def calculate_edot(ps0, ps1) -> float:
    """Tidal energy dissipation rate Edot in ps0 due to ps1."""
    k2 = K2_EARTH
    try: tau = ps0.params["tau"]
    except AttributeError: tau = ps1.params["tau"]
    sat_a, sat_period, e_sat, _, _ = cartesian_to_orbital(ps0, ps1)
    Z = 3 * G**2 * k2 * ps1.m**2 * (ps0.m + ps1.m) * ps0.r**5 / sat_a**9 * tau
    edot = 0.0
    if e_sat < 1.0:
        beta = np.sqrt(1.0 - e_sat**2)
        edot = Z / beta**15 * (eccentricity_function(1, e_sat) - eccentricity_function(2, e_sat) ** 2 / eccentricity_function(5, e_sat))
    return edot

# ---------- Spin-state tracking ----------
def calculate_spin_state(ps, t: float, filename: str, a_crit: float) -> bool:
    """Record planet spin/orbit state. Return True if Roche/unbound/escaped."""
    inside_roche = False
    pl_spin_axis_inv = ps["Planet"].params["Omega"]
    pl_orbit_normal = ps["Planet"].hvec
    rot_to_orbit_frame = rebound.Rotation.to_new_axes(newz=pl_orbit_normal)
    pl_spin_axis_orb = rot_to_orbit_frame * pl_spin_axis_inv
    mag, obliquity, phase = rebound.xyz_to_spherical(pl_spin_axis_orb)
    spin_period_pl_hours = 2 * np.pi / mag * 365.25 * 24.0
    lunar_a, lunar_orb_period, e_m, i_m, omg_m = cartesian_to_orbital(ps["Planet"], ps["Moon"])
    pl_a, pl_orb_period, pl_e, i_p, omg_p = cartesian_to_orbital(ps["Star"], ps["Planet"])
    i_p_deg, omg_p_deg = np.degrees(i_p), np.degrees(omg_p)
    sx, sy, sz = pl_spin_axis_inv
    alpha = calculate_alpha(spin_period_pl_hours, STAR_MASS, ps["Planet"].m, ps["Planet"].a)
    out_line = f"{t:1.5e},{lunar_a / R_PL:1.4f},{e_m:1.5f},{pl_a:1.4f},{pl_e:1.5f},{alpha:1.6f},{spin_period_pl_hours:1.6f},{np.degrees(obliquity):1.6f},{np.degrees(phase):1.6f},{i_p_deg:1.6f},{omg_p_deg:1.6f},{sx:1.8e},{sy:1.8e},{sz:1.8e}\n"
    write_to_file(filename, out_line, "foot")
    if lunar_a * (1 - e_m) < D_ROCHE: inside_roche = True
    if e_m > 1.0: inside_roche = True
    if lunar_a > a_crit: inside_roche = True
    return inside_roche

# ---------- Restart / skip logic ----------
def check_skip_or_restart(p_idx: int) -> int:
    """-1: run from scratch/restart, 0: restart from saved state, 1: skip (done/unstable)."""
    m_p, a_p = PARAMS[p_idx]
    subdir = os.path.join(OUTDIR, f"Mass[{m_p:1.2f}]/ap[{a_p:1.3f}]/")
    fname = os.path.join(subdir, f"reb_tide[{m_p:1.2f},{a_p:1.3f}]_{TAU:d}.txt")
    safname = os.path.join(subdir, f"sa[{m_p:1.2f},{a_p:1.3f}]_{TAU:d}.bin")
    par_fname = os.path.join(subdir, f"sim_params[{m_p:1.2f},{a_p:1.3f}]_{TAU:d}.txt")
    if not os.path.exists(fname): return -1
    if os.path.getsize(safname) / (1024 * 1024) < MIN_SNAPSHOT_MB:
        with open(par_fname, "r") as f: lines = f.readlines()
        lastline = lines[-1].strip().split(",")
        if lastline[-1].strip() == "Roche": return 1
        print(f"Binary file too small, restarting {m_p:1.2f}, {a_p:1.3f}")
        for path in (fname, safname, par_fname):
            if os.path.exists(path): os.remove(path)
        return -1
    with open(par_fname, "r") as f: lines = f.readlines()
    lastline = lines[-1].strip().split(",")
    if lastline[-1].strip() == "Roche": return 1
    simtime = float(lastline[0].strip())
    n_orb = MAX_ORBITS
    if abs(simtime - n_orb) > 1e3: return 0
    if simtime > n_orb: return 1
    return 1

# ---------- Main simulation routine ----------
def run_simulation(p_idx: int) -> None:
    """Run or resume a simulation for PARAMS[p_idx]."""
    m_p_earth, a_p = PARAMS[p_idx]
    subdir = os.path.join(OUTDIR, f"Mass[{m_p_earth:1.2f}]/ap[{a_p:1.3f}]/")
    os.makedirs(subdir, exist_ok=True)
    fname = os.path.join(subdir, f"reb_tide[{m_p_earth:1.2f},{a_p:1.3f}]_{TAU:d}.txt")
    safname = os.path.join(subdir, f"sa[{m_p_earth:1.2f},{a_p:1.3f}]_{TAU:d}.bin")
    par_fname = os.path.join(subdir, f"sim_params[{m_p_earth:1.2f},{a_p:1.3f}]_{TAU:d}.txt")
    restart_flag = "F"
    simtime = 0.0
    header_main = "#t (yr), a_m (R_p), e_m, a_p (AU), e_p, alpha (''/yr), spin_period (hr), obliquity (deg), phase (deg), i_p (deg), Omg_p (deg),S_x,S_y,S_z\n"
    header_params = "#sim time (yr), wall time (min), sim dt, tide heat (W/m^2), error\n"
    if not os.path.exists(fname):
        write_to_file(fname, header_main, "head")
        write_to_file(par_fname, header_params, "head")
    else:
        if os.path.getsize(safname) / (1024 * 1024) < MIN_SNAPSHOT_MB:
            print(f"Binary file too small, restarting {m_p_earth:1.2f}, {a_p:1.3f}")
            write_to_file(fname, header_main, "head")
            write_to_file(par_fname, header_params, "head")
            restart_flag = "F"
        else:
            with open(par_fname, "r") as f: lines = f.readlines()
            lastline = lines[-1].strip().split(",")
            if lastline[-1].strip() == "Roche": return
            simtime = float(lastline[0].strip())
            n_orb = MAX_ORBITS
            if simtime < n_orb: restart_flag = "T"
            else: return
    m_p = m_p_earth * M_EARTH
    r_hill = a_p * (m_p / (3 * STAR_MASS)) ** (1.0 / 3.0)
    a_crit = 0.4061 * (1 - 1.1257 * ECC_P) * r_hill
    a_moon = 3 * D_ROCHE
    p_moon = np.sqrt(a_moon**3 / m_p)
    sim = rebound.Simulation()
    sim.units = ("AU", "Msun", "yr")
    sim.integrator = INTEGRATOR
    sim.dt = p_moon / 40.0
    if INTEGRATOR == "trace":
        sim.ri_trace.peri_mode = "FULL_IAS15"
        sim.dt = p_moon / 5.0
    if restart_flag == "F":
        sim.add(m=STAR_MASS, hash="Star")
        sim.add(m=m_p, a=a_p, e=ECC_P, inc=np.radians(INC_P), r=R_PL, hash="Planet", primary=sim.particles["Star"])
        sim.add(m=M_MOON, a=a_moon, e=1e-6, hash="Moon", primary=sim.particles["Planet"])
        sim.move_to_com()
        ps = sim.particles
        rebx = reboundx.Extras(sim)
        ts = rebx.load_force("tides_spin")
        rebx.add_force(ts)
        spin_period_pl_yr = P_ROT_HOURS / 24.0 / 365.25
        spin_mag_pl = 2 * np.pi / spin_period_pl_yr
        ps["Planet"].params["Omega"] = rebound.spherical_to_xyz(magnitude=spin_mag_pl, theta=0.0, phi=BETA_P)
        sim.save_to_file(safname, interval=10000, delete_file=True)
    else:
        sim = rebound.Simulation(safname)
        ps = sim.particles
        lunar_a, moon_period, e_m, i_m, omg_m = cartesian_to_orbital(ps["Planet"], ps["Moon"])
        sim.dt = moon_period / 40.0
        if INTEGRATOR == "trace":
            sim.ri_trace.peri_mode = "FULL_IAS15"
            sim.dt = p_moon / 5.0
        rebx = reboundx.Extras(sim)
        ts = rebx.load_force("tides_spin")
        rebx.add_force(ts)
        reb_state = np.genfromtxt(fname, delimiter=",", comments="#")
        sx, sy, sz = reb_state[-1, -3:]
        ps["Planet"].params["Omega"] = [sx, sy, sz]
    ps["Planet"].params["k2"] = K2_EARTH / 2.0
    ps["Planet"].params["I"] = C_EARTH * m_p * R_PL**2
    ps["Planet"].params["tau"] = TAU / 3600.0 / 24.0 / 365.25
    sim.save_to_file(safname, interval=10000, delete_file=False)
    rebx = reboundx.Extras(sim)
    ts = rebx.load_force("tides_spin")
    ltot = np.array(sim.angular_momentum()) + np.array(rebx.spin_angular_momentum())
    rot = rebound.Rotation.to_new_axes(newz=ltot)
    rebx.rotate_simulation(rot)
    rebx.initialize_spin_ode(ts)
    e_prev = sim.energy()
    t0 = sim.t
    p_pl = MAX_ORBITS
    if restart_flag == "F":
        t_fin = 1e7
        times = np.concatenate((np.arange(0, 1000, 100), np.arange(1000, 10000, 1000), np.arange(10000, t_fin + 10000, 10000)))
    else:
        t_fin = p_pl - sim.t
        times = np.arange(0, t_fin + 10000, 10000)[1:]
    for t in times:
        current_sim_t = sim.t
        if current_sim_t >= p_pl: break
        start_wall = time.time()
        sim.integrate(t0 + t)
        stop_wall = time.time()
        wall_dt = (stop_wall - start_wall) / 60.0
        e_now = sim.energy()
        edot = calculate_edot(ps["Planet"], ps["Moon"])
        e_theory = e_prev + edot * (sim.t - current_sim_t)
        error = e_now - e_theory
        tide_heat = (edot * (sim.t - current_sim_t)) * M_SUN * au**2 / year**3 / (4 * np.pi * ps["Planet"].r**2 * au**2)
        time_line = f"{t0 + t:1.5e}, {wall_dt:1.4e}, {sim.dt:1.4e}, {tide_heat:1.4e}, {error:1.4e}\n"
        write_to_file(par_fname, time_line, "foot")
        if calculate_spin_state(ps, t0 + t, fname, a_crit):
            time_line_roche = f"{t0 + t:1.5e}, {wall_dt:1.4e}, {sim.dt:1.4e}, {tide_heat:1.4e}, Roche\n"
            write_to_file(par_fname, time_line_roche, "foot")
            break
        e_prev = e_now

# ---------- Parameter grid & main ----------
def generate_parameter_grid(a_i: float, a_o: float) -> List[Tuple[float, float]]:
    params: List[Tuple[float, float]] = []
    m_values = np.arange(M_MIN, M_MAX + M_STEP / 2.0, M_STEP)
    a_values = np.arange(a_i, a_o, A_STEP)
    for m in m_values:
        for a in a_values:
            params.append((np.round(m, 2), np.round(a, 3)))
    np.random.seed(42)
    np.random.shuffle(params)
    return params

def main() -> None:
    global STAR_MASS, OUTDIR, PARAMS
    if len(sys.argv) < 3:
        raise SystemExit("Usage (grid):   python tide_spin_moon_sim.py STAR_TYPE START_INDEX\nUsage (single): python tide_spin_moon_sim.py STAR_TYPE single M_EARTH A_AU")
    star_type = sys.argv[1]
    if star_type not in STAR_CONFIGS: raise SystemExit(f"Unknown STAR_TYPE '{star_type}'. Choose from {list(STAR_CONFIGS.keys())}")
    cfg = STAR_CONFIGS[star_type]
    STAR_MASS = cfg["mass"]
    OUTDIR = os.path.join(BASE_OUTPUT_DIR, f"{star_type}_runs/{INTEGRATOR}_runs/")
    os.makedirs(OUTDIR, exist_ok=True)
    mode = sys.argv[2]
    if mode == "single":
        if len(sys.argv) != 5: raise SystemExit("Single-mode usage: python tide_spin_moon_sim.py STAR_TYPE single M_EARTH A_AU")
        m_earth = float(sys.argv[3])
        a_au = float(sys.argv[4])
        PARAMS = [(m_earth, a_au)]
        decision = check_skip_or_restart(0)
        if decision == 1:
            print(f"Single IC {m_earth:.2f} M_earth, a={a_au:.3f} AU already completed or unstable.")
            return
        print(f"Running single IC: {m_earth:.2f} M_earth, a={a_au:.3f} AU")
        run_simulation(0)
        return
    # grid / parallel mode
    start_idx = int(mode)
    end_idx = start_idx + (300 if start_idx == 0 else 60)
    PARAMS = generate_parameter_grid(cfg["a_i"], cfg["a_o"])
    PARAMS = PARAMS[start_idx:end_idx]
    idx = 0
    while idx < len(PARAMS):
        decision = check_skip_or_restart(idx)
        if decision == 1:
            PARAMS.pop(idx)
        else:
            m_p, a_p = PARAMS[idx]
            print(f"{m_p:1.2f}, {a_p:1.3f}")
            idx += 1
    if not PARAMS:
        print("No simulations to run for this index range.")
        return
    ncores = min(N_CORES_TARGET, mp.cpu_count())
    with mp.Pool(processes=ncores) as pool:
        pool.map(run_simulation, range(len(PARAMS)))

if __name__ == "__main__":
    main()
