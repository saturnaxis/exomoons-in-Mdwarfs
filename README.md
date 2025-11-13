# exomoons-in-Mdwarfs
This repository contains processed data and Python scripts used in the article:

**“Tidally Torn: Why the Most Common Stars May Lack Large, Habitable-Zone Moons”**  
(Patel, Quarles, Weinberg & Cuntz 2025)

The paper investigates whether Earth-mass planets in the habitable zones of M dwarfs can stably host large moons, and how tidal evolution driven by the host star influences long-term moon survival.

---

## Repository Structure

```
Mdwarf-tides/
│
├── data/                # Processed simulation outputs used by the plotting scripts
├── python-scripts/      # Figure-generation scripts used in the paper
└── Figs/                # Output directory for generated figures
```

All Python scripts assume the directory structure produced by cloning the repository:

```
git clone https://github.com/<user>/Mdwarf-tides.git
```

They expect the `data` directory to be located at the same level as `python-scripts` and will write figures to `Figs/` automatically.

---

## Python Scripts

The **`python-scripts/`** folder contains a number of standalone Python files that reproduce the figures from the paper.

Each script follows the naming convention:

```
plot_FigXX.py
```

where **XX** corresponds to the figure number in the manuscript.

### Usage

These scripts require **no command-line arguments**. Run them from inside the `python-scripts` directory:

```bash
python plot_FigXX.py
```

Each script will:

- Load the appropriate data from `../data/`
- Produce the corresponding figure
- Save the output to `../Figs/`

### Dependencies

These scripts rely on the following Python packages:

- `numpy`
- `scipy` (version ≥ 1.4)
- `matplotlib`

You can install all dependencies with:

```bash
pip install numpy scipy matplotlib
```

or via conda:

```bash
conda install numpy scipy matplotlib
```

---

## Attribution

If you use this repository, data products, or analysis tools in your research, please cite:

```
@ARTICLE{Patel2025,
       author = {{Patel}, Shaan D. and {Quarles}, Billy and {Weinberg}, Nevin N. and {Cuntz}, Manfred},
        title = "{Tidally Torn: Why the Most Common Stars May Lack Large, Habitable-Zone Moons}",
      journal = {arXiv e-prints},
     keywords = {Earth and Planetary Astrophysics, Solar and Stellar Astrophysics},
         year = 2025,
        month = nov,
          eid = {arXiv:2511.03625},
        pages = {arXiv:2511.03625},
          doi = {10.48550/arXiv.2511.03625},
 archivePrefix = {arXiv},
        eprint = {2511.03625},
  primaryClass = {astro-ph.EP},
        adsurl = {https://ui.adsabs.harvard.edu/abs/2025arXiv251103625P},
       adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
```

---
