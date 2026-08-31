# MDCheck

[![CI](https://github.com/amonreal/mdcheck/actions/workflows/test.yml/badge.svg)](https://github.com/amonreal/mdcheck/actions)
[![PyPI version](https://img.shields.io/pypi/v/mdcheck.svg?color=blue)](https://pypi.org/project/mdcheck/)
[![Python versions](https://img.shields.io/pypi/pyversions/mdcheck.svg)](https://pypi.org/project/mdcheck/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1234567.svg)](https://doi.org/10.5281/zenodo.1234567)

> **Automated Convergence, Statistical Inefficiency, and Reproducibility Assessment for Molecular Dynamics Simulations.**

---

## Overview

**MDCheck** is an open-source scientific toolkit that solves a universal methodological need in biomolecular and materials simulations: **certifying whether a molecular dynamics trajectory has converged, reached equilibrium, and accumulated sufficient statistically independent observations for publication.**

Instead of manually inspecting plots or guessing equilibration cutoffs, `mdcheck` analyzes raw timeseries (`.xvg`, `.csv`, `.dat`, `.log`) with a single command and delivers:

- 🎯 **Automated Equilibration Detection ($t_{\text{eq}}$)** via statistical inefficiency minimization.
- ⏱️ **Integrated Autocorrelation Time ($\tau_{\text{int}}$)** & **Statistical Inefficiency ($g$)** using Madras-Sokal self-consistent windowing.
- 🔢 **Effective Sample Size ($N_{\text{eff}} = N / g$)** to ensure statistically rigorous error estimation.
- 🔄 **Multi-Replica Reproducibility Matrix ($R_1 \text{ vs } R_2 \text{ vs } R_3$)** based on Jensen-Shannon Divergence (JSD) and Essential Subspace Overlap (RMSIP).
- 📉 **Linear & CUSUM Systematic Drift Diagnostics** and Flyvbjerg-Petersen block averaging.
- 🚦 **Quality Certification Badges (`PASS` / `WARNING` / `FAIL`)** with unambiguous diagnostic messages.
- 📑 **Publication-Ready Outputs**: Interactive self-contained `report.html`, vector plots (SVG/PDF/PNG 300 DPI), LaTeX summary tables (`.tex`), and a draft **Methods & Supporting Information** text snippet with automated **BibTeX citations**.

```
  Simulations (.xvg, .csv, .dat)
               │
               ▼
  ┌───────────────────────────────────────────────────────────┐
  │                         MDCheck                           │
  │  ├── Auto Equilibration (max N_eff)                       │
  │  ├── Autocorrelation & Inefficiency (tau_int, g)          │
  │  ├── Multi-Replica Overlap (Jensen-Shannon, RMSIP)        │
  │  └── Drift Detection & Block Averaging                    │
  └───────────────────────────────────────────────────────────┘
               │
               ▼
  ┌───────────────────────────────────────────────────────────┐
  │                   Publication Deliverables                │
  │  ├── report.html (Interactive Dashboard & Badges)         │
  │  ├── mdcheck_convergence_overview.pdf/svg/png             │
  │  ├── mdcheck_summary_table.tex / .csv                     │
  │  ├── methods_snippet.txt (Ready for Manuscript)           │
  │  └── citation.bib (BibTeX Reference)                      │
  └───────────────────────────────────────────────────────────┘
```

---

## Installation

### From PyPI
```bash
pip install mdcheck
```

### From Source (Development Mode)
```bash
git clone https://github.com/amonreal/mdcheck.git
cd mdcheck
pip install -e .[dev]
```

---

## Quickstart (CLI)

### 1. Test Demo Mode (Instant Synthetic Multi-Replica Simulation)
```bash
mdcheck demo -o my_demo_results/
```
Open `my_demo_results/report.html` in any web browser to see the interactive report!

### 2. Assess GROMACS XVG Trajectory
```bash
mdcheck assess -i rmsd.xvg gyrate.xvg energy.xvg -o md_quality_report/
```

### 3. Assess Multi-Replica Convergence (R1, R2, R3)
```bash
mdcheck assess -i rep1_rmsd.xvg -r rep2_rmsd.xvg rep3_rmsd.xvg -o replica_assessment/
```

---

## Python API Usage

```python
import numpy as np
from mdcheck import assess_trajectory_quality
from mdcheck.reporters import generate_publication_figures, generate_manuscript_assets, generate_html_report

# Load or define your timeseries (e.g. Backbone RMSD over 100 ns)
time_coords = np.linspace(0, 100, 2000)  # ns
rmsd_series = ... # 1D numpy array

# Assess simulation quality
report = assess_trajectory_quality(
    timeseries_dict={"Backbone_RMSD": rmsd_series},
    time_coords=time_coords
)

print(f"Overall Quality Status: {report.overall_status}")
print(f"Equilibration Time: {report.observables['Backbone_RMSD'].t_eq_time:.2f} ns")
print(f"Effective Sample Size (N_eff): {report.observables['Backbone_RMSD'].n_eff:.0f}")

# Export publication figures and LaTeX tables
generate_publication_figures({"Backbone_RMSD": rmsd_series}, time_coords, report, "output_dir/")
generate_manuscript_assets(report, "output_dir/")
generate_html_report(report, "output_dir/report.html")
```

---

## Scientific Foundations & Methodology

### 1. Automated Equilibration Detection ($t_{\text{eq}}$)
The initial non-equilibrium transient phase is automatically identified by maximizing the total effective sample size in the subsequent production interval:
$$\hat{t}_{\text{eq}} = \arg\max_{t_0} N_{\text{eff}}(t_0) = \arg\max_{t_0} \frac{N - t_0}{g(t_0)}$$

### 2. Autocorrelation & Statistical Inefficiency ($g$)
MD frames are temporally correlated. MDCheck computes the integrated autocorrelation time $\tau_{\text{int}}$ using the Madras-Sokal self-consistent cutoff window:
$$\tau_{\text{int}} = \frac{1}{2} + \sum_{k=1}^{M} C(k), \quad M \ge 6 \tau_{\text{int}}$$
$$g = 1 + 2\tau_{\text{int}}, \quad N_{\text{eff}} = \frac{N_{\text{prod}}}{g}$$

### 3. Multi-Replica Conformational Overlap
Conformational consistency between independent trajectories ($R_1, R_2, R_3$) is evaluated via the square-root of Jensen-Shannon Divergence ($\mathrm{JSD} \in [0, 1]$) and Root-Mean-Square Inner Product (RMSIP) across essential PCA subspaces:
$$\mathrm{JSD}(P \parallel Q) = \frac{1}{2} D_{\text{KL}}(P \parallel M) + \frac{1}{2} D_{\text{KL}}(Q \parallel M)$$

---

## Citation

If you use MDCheck to evaluate trajectory convergence, equilibration, statistical inefficiency, or replica reproducibility in your research, please cite:

```bibtex
@software{monreal2026mdcheck,
  author = {Monreal-Hern{\'a}ndez, Andre},
  title = {{MDCheck: Automated Convergence, Statistical Inefficiency, and Reproducibility Assessment for Molecular Dynamics Simulations}},
  year = {2026},
  version = {1.0.0},
  publisher = {Zenodo},
  url = {https://github.com/amonreal/mdcheck}
}
```

---

## License

This project is licensed under the **MIT License** - see the [LICENSE](LICENSE) file for details.
