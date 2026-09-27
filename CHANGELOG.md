# Changelog

## 1.1.0 (unreleased)

Changes prompted by the validation in `validation/`, which compares MDCheck with pymbar 4.0.3, pyblock and
ArviZ on AR(1) processes with known answers and on OpenMM simulations of Lennard-Jones argon and TIP3P water.

### Fixed
- **Confidence interval of the mean.** The moving-block bootstrap with blocks of length 2·τ_int was about
  25% too narrow for correlated data, covering the true mean in only 82–90% of AR(1) cases. The reported
  interval is now `mean ± t(N_eff − 1) · s · sqrt(g/N)` (`inefficiency_corrected_ci`), which gives 0.90–0.99
  coverage. `bootstrap_ci` is kept for API compatibility.
- **Replica consistency.** A histogram Jensen–Shannon distance above 0.15 set the FAIL status. On replicas of
  the same ensemble this raised an alarm for 100% of short-replica groups. The status now comes from
  Cochran's Q test on autocorrelation-corrected standard errors (χ²₍ₖ₋₁₎; PASS for p ≥ 0.05, WARNING for
  0.01 ≤ p < 0.05, FAIL for p < 0.01). JSD, KS and Wasserstein distances remain as descriptors. The result
  now also includes `cochran_q`, `heterogeneity_p_value`, `i_squared` and `replica_standard_errors`.
- **AMBER parser.** `parse_amber_dat` failed with pandas ≥ 3 (`delim_whitespace` was removed) and discarded
  the `#Frame ...` header written by cpptraj.
- **README.** `g` was written as `1 + 2·τ_int`; the code computes, correctly, `g = 2·τ_int = 1 + 2·Σ C(k)`.
  The README also advertised a CUSUM drift test that did not exist.

### Added
- LAMMPS `log.lammps` parser (`parse_lammps_log`), which merges consecutive thermo blocks and is detected
  automatically.
- A reproducible validation suite in `validation/`: scripts, raw MD time series, results, tables and figures.

### Changed
- `detect_equilibration` searches a grid of stride `max(10, N // 200)` and refines it hierarchically down to
  one sample (about 250 evaluations regardless of N). It returns the same t_eq as the exhaustive stride-10
  search on all 120 OpenMM replicas of the validation, and runs about 4× faster at N = 10⁴.
  Pass `step_search=10` to reproduce the old search.

## 1.0.0 (2026-08-31)

Initial release.
