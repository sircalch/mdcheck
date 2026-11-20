# Checking one trajectory

## Input files

MDCheck reads a time series from a file, detecting the format by extension and content:

| Format | Notes |
|---|---|
| GROMACS `.xvg` | Column names are read from the legends |
| `.csv`, `.tsv`, `.dat`, `.txt` | First column is time; the others are observables |
| LAMMPS and OpenMM log files | Thermodynamic output is read as columns |

The repository ships a GROMACS example, `sample_data/rmsd.xvg`:

```text
# GROMACS XVG generated for MDCheck
@ title "Backbone RMSD"
@ xaxis label "Time (ns)"
@ yaxis label "RMSD (nm)"
```

## Command line

```bash
mdcheck assess -i sample_data/rmsd.xvg -o rmsd_report/
```

The summary shows, for each observable, the equilibration time, the effective sample size, the mean with its 95 % confidence interval and the status:

```text
 [RESULT] Overall Trajectory Certification: WARNING
 [SCORE]  SIMULATION QUALITY = BORDERLINE / ACCEPTABLE WITH CAUTION
 * Backbone            : t_eq = 14.11 | N_eff =    105 | Mean = 0.218 [0.215, 0.221] | Status: WARNING
```

Here MDCheck discards the first 14.11 time units as equilibration. The remaining samples are strongly autocorrelated, so the 1,500 frames carry the information of about 105 independent ones, which is below the 500 needed for `PASS` and above the 100 below which the status would be `FAIL`. Open `rmsd_report/report.html` to see the diagnostic message for each check.

## Comparing replicas from the command line

With several files, the first one is the primary trajectory, and the others are matched to it **by column name**: an observable present in more than one file is compared across them.

```bash
mdcheck assess -i run1_rg.xvg -r run2_rg.xvg run3_rg.xvg -o rg_report/
```

For three replicas stored as columns of a single file, use the Python API shown in the [worked example](worked-example.md).

## Python API

```python
import numpy as np
from mdcheck import assess_trajectory_quality
from mdcheck.parsers.generic import load_timeseries_file

time, series, meta = load_timeseries_file("sample_data/rmsd.xvg")
report = assess_trajectory_quality(timeseries_dict=series, time_coords=time)

print(report.overall_status)
for name, obs in report.observables.items():
    print(name, obs.status, round(obs.t_eq_time, 2), round(obs.n_eff))
```

See [API reference](../api.md) for all arguments and [How the checks work](../quality-control.md) for the meaning of each status.
