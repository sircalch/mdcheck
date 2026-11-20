# How the checks work

Each observable gets a status (`PASS`, `WARNING` or `FAIL`) with a short diagnostic; the replica comparison gets its own. The overall status is the worst of them. A `WARNING` asks you to look at a result; a `FAIL` means the average should not be reported as it is.

## Equilibration time

The first part of a simulation reflects the starting structure, not the equilibrium ensemble. MDCheck tries candidate start times \(t_0\) and, for each one, estimates the statistical inefficiency \(g(t_0)\) of the remaining series. It keeps the \(t_0\) that **maximises the number of independent samples** in the production region:

\[
t_{\text{eq}} = \arg\max_{t_0} N_{\text{eff}}(t_0) = \arg\max_{t_0} \frac{N - t_0}{g(t_0)}
\]

This is the same criterion used by `pymbar.timeseries.detect_equilibration`; [Validation](validation.md) compares both.

Function: [`detect_equilibration`](api.md).

## Autocorrelation, inefficiency and effective sample size

Consecutive MD frames are correlated, so \(N\) frames carry less information than \(N\) independent samples. MDCheck computes the autocorrelation function by FFT and the integrated autocorrelation time \(\tau_{\text{int}}\) with the self-consistent Madras–Sokal window (default `c_window = 6`):

\[
\tau_{\text{int}} = \tfrac12 + \sum_{k=1}^{M} C(k), \qquad g = 1 + 2\tau_{\text{int}}, \qquad N_{\text{eff}} = N / g
\]

Functions: [`integrated_autocorrelation_time`](api.md), [`statistical_inefficiency`](api.md), [`effective_sample_size`](api.md).

## Confidence interval of the mean

The interval for the mean of the production region uses a *t* distribution with \(N_{\text{eff}} - 1\) degrees of freedom:

\[
\bar{x} \pm t_{N_{\text{eff}}-1}\, s \sqrt{g / N}
\]

In the AR(1) benchmarks this interval reaches its nominal coverage, whereas a moving-block bootstrap with blocks of about \(2\tau_{\text{int}}\) was too narrow; see [Validation](validation.md). A moving-block bootstrap ([`bootstrap_ci`](api.md)) and Flyvbjerg–Petersen block averaging ([`block_averaging`](api.md)) are also available.

Function: [`inefficiency_corrected_ci`](api.md).

## Drift and stationarity of the production region

A linear fit of the production region gives its relative change across the run.

| Status | Condition (defaults) |
|---|---|
| `PASS` | slope not significant (\(p > 0.05\)), or relative change \(\le 5\) % |
| `WARNING` | relative change between 5 % and 15 % |
| `FAIL` | relative change above 15 % |

The thresholds are arguments of [`assess_drift`](api.md). A Geweke diagnostic ([`geweke_diagnostic`](api.md)) additionally compares the first 10 % and the last 50 % of the production region; if it rejects stationarity, the observable is at least a `WARNING`.

## Status of each observable

| Condition (defaults) | Status |
|---|---|
| \(N_{\text{eff}} < 100\) | `FAIL` |
| \(100 \le N_{\text{eff}} < 500\) | `WARNING` (marginal sampling) |
| More than 50 % of the run discarded as equilibration | `WARNING` |
| Drift `FAIL` | `FAIL` |
| Drift `WARNING` or Geweke rejects stationarity | `WARNING` |

The thresholds `min_neff_pass`, `min_neff_warn` and `max_discard_warn` are arguments of [`assess_trajectory_quality`](api.md). They are defaults for typical MD observables; adjust them to the precision your study needs.

## Reproducibility across replicas

When replicas are provided, MDCheck compares their means with Cochran's \(Q\) test for heterogeneity. Each replica mean has a standard error based on its own effective sample size, so the test accounts for autocorrelation.

| Status | Condition (defaults) |
|---|---|
| `PASS` | \(p \ge 0.05\) |
| `WARNING` | \(0.01 \le p < 0.05\) |
| `FAIL` | \(p < 0.01\) |

The report also gives the pairwise Jensen–Shannon distances, Kolmogorov–Smirnov tests and Wasserstein distances between replicas, and \(I^2\). The Jensen–Shannon distance is reported with a level (`LOW`, `MODERATE`, `HIGH`) but **does not set the status**: histogram distances between finite correlated samples are rarely zero even for replicas from the same ensemble (see [Validation](validation.md)).

Functions: [`assess_replica_consistency`](api.md), [`jensen_shannon_distance`](api.md), [`subspace_overlap_rmsip`](api.md).
