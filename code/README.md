# Tidal Spin–Orbit Evolution with REBOUND/REBOUNDx

This repository contains a Python script to simulate the tidal evolution and spin–orbit dynamics of a planet–moon system using [REBOUND](https://rebound.readthedocs.io/en/latest/) and [REBOUNDx](https://reboundx.readthedocs.io/en/latest/).

The script can:

- Scan a grid of planet masses and semi-major axes in parallel.
- Or run a single initial condition for high-detail investigation.
- Switch between M-dwarf host star types (M0, M2, M4) with predefined stellar masses and orbital parameter ranges.
- Save orbital elements, spin-state evolution, tidal heating, and restart files.
- Automatically restart unfinished runs and skip completed or unstable ones.

---

## 🔔 Citation Reminder

If you reuse or adapt any part of this code, **please cite the following paper**:

**Patel et al. (2025)**  
*Tidally Torn: Why the Most Common Stars May Lack Large, Habitable-Zone Moons*  
https://arxiv.org/abs/2511.03625

---

## 1. Requirements

- Python 3.8+
- `rebound`
- `reboundx`
- `numpy`
- `scipy`

Install:

```bash
pip install rebound reboundx numpy scipy
2. Script Overview
Primary script:

tide_spin_moon_sim.py

Key features:

Uses the tides_spin module from REBOUNDx.

Planet labeled "Planet".

Moon labeled "Moon".

Host star labeled "Star".

Uses binary REBOUND snapshots for restarts.

Records orbital elements, spin-axis components, and tidal heating.

3. Supported Star Types
The host star type is selected using the first command-line argument: M0, M2, or M4.

Star Type	Stellar Mass (Msun)	a_i (AU)	a_o (AU)
M0	0.57	0.27	0.521
M2	0.44	0.17	0.351
M4	0.23	0.08	0.181

Planet–mass grid:

0.8 → 2.0 Earth masses

Step 0.01

Semi-major axis range set by star type

Step 0.001 AU

Randomized using seeded shuffle

4. Output Structure
If you set:

bash
Copy code
export TIDE_SIM_HOME=/path/to/storage
the script writes outputs inside this directory.
Otherwise, it uses the current working directory.

Example structure:

css
Copy code
<BASE_OUTPUT_DIR>/
  M0_runs/trace_runs/
    Mass[1.20]/ap[0.250]/
      reb_tide[1.20,0.250]_698.txt
      sa[1.20,0.250]_698.bin
      sim_params[1.20,0.250]_698.txt
Output Files
reb_tide[…] .txt
Spin state & orbital evolution:

lua
Copy code
t, a_m/R_p, e_m, a_p, e_p, alpha(''/yr), spin_period(hr),
obliquity(deg), phase(deg), i_p(deg), Omg_p(deg),
Sx, Sy, Sz
sim_params[…].txt
Integration log:

scss
Copy code
sim_time(yr), wall_time(min), dt, tide_heat(W/m^2), error
If the moon becomes unstable (Roche limit, unbound, or escapes the Hill sphere):

Copy code
..., Roche
sa[…].bin
REBOUND snapshot file (used for restarts).

5. Usage
Grid Mode (Parallel Execution)
Run a block of simulations:

bash
Copy code
python tide_spin_moon_sim.py STAR_TYPE START_INDEX
Examples:

bash
Copy code
python tide_spin_moon_sim.py M0 0
python tide_spin_moon_sim.py M0 300
python tide_spin_moon_sim.py M2 0
Chunk rules:

If START_INDEX == 0 → run 300 cases

If START_INDEX > 0 → run 60 cases

The script:

Generates parameter grid

Slices by index

Skips completed or unstable runs

Restarts incomplete ones

Runs remaining jobs in parallel (up to 60 cores)

Single Initial Condition Mode
For debugging or detailed inspection:

bash
Copy code
python tide_spin_moon_sim.py STAR_TYPE single M_EARTH A_AU
Example:

bash
Copy code
python tide_spin_moon_sim.py M2 single 1.20 0.250
This runs exactly one simulation, with restart support.

6. SLURM Parallel Example
bash
Copy code
#!/bin/bash
#SBATCH --job-name=tides_M0
#SBATCH --output=logs/tides_M0_%A_%a.out
#SBATCH --error=logs/tides_M0_%A_%a.err
#SBATCH --array=0-9
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=60
#SBATCH --time=48:00:00

module load python/3.10
source /path/to/venv/bin/activate

export TIDE_SIM_HOME=/scratch/$USER/tides_output

CHUNK_START=$(( SLURM_ARRAY_TASK_ID * 60 ))
python tide_spin_moon_sim.py M0 $CHUNK_START
Adjust CPU count, array size, or wall time to match your cluster’s policies.

7. Customization Options
You may wish to modify:

N_CORES_TARGET — number of processes for parallel grid runs

P_ROT_HOURS — initial planetary rotation period

M_MIN, M_MAX, M_STEP, A_STEP — grid dimensions

INTEGRATOR — "trace" or "ias15"

For more sophisticated interior or tidal models, extend the relevant functions and document changes.

8. Additional Citations
If you use REBOUND or REBOUNDx, please cite:

Rein & Liu (2012), A&A 537, A128 — REBOUND

Tamayo et al. (2020), MNRAS 491, 2885 — REBOUNDx

If you reuse or adapt this repository, please also cite:

Patel et al. (2025) — https://arxiv.org/abs/2511.03625