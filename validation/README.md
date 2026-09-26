# MDCheck validation

This directory reproduces every table and figure of the MDCheck 1.1.0 validation manuscript.

The benchmarks use two kinds of data where the right answer is known:

1. **AR(1) processes.** The statistical inefficiency g = (1+φ)/(1−φ) and the mean (0) are exact.
2. **OpenMM simulations.** Lennard-Jones argon started from an fcc lattice, and TIP3P water in the NPT ensemble. Long independent runs define the reference means.

pymbar 4.0.3 and pyblock serve as reference implementations.

## Environment

```bash
python -m venv venv && . venv/bin/activate      # Windows: venv\Scripts\activate
pip install -e ..  pymbar==4.0.3 pyblock openmm pandas matplotlib
```

The OpenMM runs use the fastest available platform (CUDA → OpenCL → CPU). The results committed here were produced with OpenCL, OpenMM 8.6.1, Python 3.12 and NumPy 2.5.

## Pipeline (all seeds fixed)

| Step | Command | Output | Runtime (OpenCL GPU) |
|---|---|---|---|
| AR(1) benchmark | `python validate_ar1.py --reps 100` | `results/ar1_*.csv`, `results/equilibration_*.csv` | ~60 min (CPU) |
| LJ argon | `python validate_md_openmm.py --system lj --n-short 48` | `results/md_lj_*` (incl. raw series `.npz`) | ~4 min |
| TIP3P, 100 ps replicas | `python validate_md_openmm.py --system water` | `results/md_water_*` | ~25 min |
| TIP3P, 300 ps replicas | `python validate_md_openmm.py --system water --short-ps 300 --long-ps 1200 --ref-discard-ps 100 --seed 4321 --out results_water300` | `results_water300/` | ~60 min |
| Final MD analysis | `python analyze_md.py` and `python analyze_md.py --results results_water300` | `md_final_*.csv` | seconds |
| Estimator comparison | `python compare_estimators.py` | `results/estimator_comparison*.csv` | ~5 min |
| Tables | `python make_tables.py` | `tables/*.tex` | seconds |
| Figures | `python make_figures.py` | `figures/fig1–fig4` | ~1 min |

## What each check measures

- **Coverage.** The fraction of replicas whose 95% confidence interval contains the true mean (AR(1)) or the reference mean (MD). For MD, the uncertainty of the reference is added in quadrature.
- **Equilibration.** The t_eq chosen by MDCheck compared with `pymbar.timeseries.detect_equilibration`, plus the bias of the mean left after discarding it.
- **Replica false alarms.** Triples of replicas drawn from the same ensemble should pass. The Cochran Q test has a nominal false-alarm rate of 5%. The Jensen-Shannon criterion of version 1.0.0 is reported for comparison.

## Findings that changed the code (1.0.0 → 1.1.0)

- **Confidence interval.** The moving-block bootstrap with blocks of length 2·τ_int was ~25% too narrow (0.82–0.90 coverage). It was replaced by the t interval on N_eff − 1 degrees of freedom.
- **Replica check.** The Jensen-Shannon threshold of 0.15 raised an alarm for 100% of same-ensemble replica groups. It was replaced by Cochran's Q test.
