# Validation

The `validation/` folder of the repository contains the scripts and the result files behind these benchmarks, so they can be repeated. They use two kinds of data where the right answer is known: AR(1) processes, for which the statistical inefficiency \(g = (1+\varphi)/(1-\varphi)\) and the mean are exact, and OpenMM simulations of Lennard-Jones argon and TIP3P water, where long independent runs define the reference means. pymbar 4.0.3 and pyblock serve as reference implementations.

## Findings that changed the code (1.0.0 to 1.1.0)

- **Confidence interval.** The moving-block bootstrap with blocks of length \(2\tau_{\text{int}}\) was about 25 % too narrow, with a coverage of 0.82 to 0.90 for a nominal 95 %. It was replaced by the \(t\) interval on \(N_{\text{eff}} - 1\) degrees of freedom, which reaches 0.94 to 0.99 on the AR(1) benchmark with \(\varphi\) from 0.8 to 0.99 (`results/ar1_inefficiency_coverage.csv`).
- **Replica check.** The Jensen–Shannon threshold of 0.15 raised an alarm for 100 % of replica groups drawn from the same ensemble. It was replaced by Cochran's \(Q\) test, which has a controlled false-alarm rate.

## What each benchmark measures

| Quantity | Reference | Files |
|---|---|---|
| Statistical inefficiency \(g\) | Exact value of AR(1) processes; pymbar; pyblock | `compare_estimators.py`, `results/ar1_inefficiency_coverage.csv`, `results/estimator_comparison.csv` |
| Coverage of the 95 % interval | AR(1) exact mean; MD reference means | `results/ar1_*.csv`, `results/md_final_summary.csv` |
| Equilibration time | `pymbar.timeseries.detect_equilibration`, and the bias left after discarding | `results/equilibration_summary.csv`, `results/md_final_summary.csv` |
| Replica false alarms | Triples of replicas from the same ensemble | `compare_replica_tests.py`, `results/replica_tests_summary.csv` |
| Runtime | Scaling with series length | `benchmark_runtime.py`, `results/runtime_benchmark.csv` |

## Reading the results with care

The benchmarks show where the checks work and where they are limited. For example, in the TIP3P water benchmark the coverage of the nominal 95 % interval was 0.71 to 0.79 (`results/md_final_summary.csv`), below the AR(1) result: when a series is only a few correlation times long, any estimate of \(\tau_{\text{int}}\) is uncertain and the intervals are optimistic. This is why a low \(N_{\text{eff}}\) is reported as a warning in the first place.

The full protocol, seeds, software versions and runtimes are in `validation/README.md`.
