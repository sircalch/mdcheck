# MDCheck

**Convergence, equilibration and statistical-inefficiency checks for molecular dynamics simulations.**

An average from a molecular dynamics (MD) trajectory is only as good as the evidence that the trajectory was equilibrated and long enough. MDCheck reads time series produced by your simulation (RMSD, radius of gyration, energy, density, temperature and so on) and checks the points that are usually decided by eye:

- **Equilibration time**, chosen to maximise the number of statistically independent samples of the production region.
- **Statistical inefficiency and effective sample size** from the integrated autocorrelation time.
- **Drift** in the production phase.
- **Confidence intervals** for the mean that account for autocorrelation.
- **Reproducibility across replicas**, with a test that compares the replicas with each other.

Every decision is recorded, and results are written as a self-contained HTML report, figures, LaTeX and CSV tables and a draft methods paragraph.

```bash
pip install mdcheck
mdcheck demo -o my_first_report/
```

Open `my_first_report/report.html` in a browser to see a report built from a synthetic multi-replica dataset.

## Where to go next

- [Installation](installation.md)
- [Checking one trajectory](tutorials/single-trajectory.md)
- [Worked example: three replicas](tutorials/worked-example.md)
- [How the checks work](quality-control.md): what each status means and how to change the thresholds
- [Validation](validation.md): benchmarks with known answers and comparison with pymbar
- [API reference](api.md)

## Scope

MDCheck works on **time series that you extract from a trajectory** (for example with GROMACS tools, LAMMPS or OpenMM logs). It does not run simulations, read coordinate trajectories or validate a force field. A `PASS` means the checks implemented here found no problem; it does not mean the simulation protocol is correct.

## Getting help

Questions and bug reports are welcome in the [issue tracker](https://github.com/sircalch/mdcheck/issues). See [Support and contributing](support.md).
