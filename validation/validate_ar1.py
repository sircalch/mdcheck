"""
Known-answer validation of MDCheck against pymbar and pyblock.

Benchmarks (all on AR(1) processes x_t = phi * x_{t-1} + e_t, where the statistical
inefficiency g = (1 + phi) / (1 - phi) and the mean (0) are known exactly):

  1. Statistical inefficiency g: MDCheck vs pymbar.timeseries vs analytic.
  2. 95% confidence-interval coverage of the true mean:
     naive (s/sqrt(N)), MDCheck v1.0.0 moving-block bootstrap, MDCheck g-corrected t interval, g-corrected SE (MDCheck g and
     pymbar g), and pyblock (Flyvbjerg-Petersen with Wolff/Lee optimal block).
  3. Equilibration detection on AR(1) plus an exponentially decaying initial transient:
     t_eq and residual bias of MDCheck vs pymbar.timeseries.detect_equilibration.

Usage:
    python validation/validate_ar1.py [--reps 200] [--out validation/results]

Outputs CSV tables (and a JSON with the software versions) under --out.
"""
import argparse
import json
import os
import platform
import sys
import time

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import mdcheck  # noqa: E402
from mdcheck.core.autocorrelation import integrated_autocorrelation_time  # noqa: E402
from mdcheck.core.bootstrap import bootstrap_ci, inefficiency_corrected_ci  # noqa: E402
from mdcheck.core.equilibration import detect_equilibration  # noqa: E402

import pymbar  # noqa: E402
from pymbar import timeseries  # noqa: E402
import pyblock  # noqa: E402

Z95 = 1.959963984540054


def ar1(phi, n, rng, burn=None):
    """Stationary AR(1) with unit marginal variance."""
    sigma_e = np.sqrt(1.0 - phi * phi)
    x = np.empty(n)
    x[0] = rng.normal()  # draw from the stationary distribution
    e = rng.normal(0.0, sigma_e, n)
    for t in range(1, n):
        x[t] = phi * x[t - 1] + e[t]
    return x


def g_true(phi):
    return (1.0 + phi) / (1.0 - phi)


def pyblock_se(x):
    """Standard error of the mean at pyblock's optimal block size (Wolff/Lee criterion)."""
    reblock = pyblock.blocking.reblock(x)
    opt = pyblock.blocking.find_optimal_block(len(x), reblock)
    if not opt or opt[0] is None or (isinstance(opt[0], float) and np.isnan(opt[0])):
        return np.nan
    return float(reblock[opt[0]].std_err)


def bench_inefficiency_and_coverage(phis, n, reps, rng):
    rows = []
    for phi in phis:
        for r in range(reps):
            x = ar1(phi, n, rng)
            m = x.mean()
            s = x.std(ddof=1)

            _, g_md, _ = integrated_autocorrelation_time(x)
            g_pm = float(timeseries.statistical_inefficiency(x))

            _, lo_bs, hi_bs = bootstrap_ci(x, random_state=int(rng.integers(1 << 31)))
            _, lo_ci, hi_ci = inefficiency_corrected_ci(x, g=g_md)
            se_naive = s / np.sqrt(n)
            se_g_md = s * np.sqrt(g_md / n)
            se_g_pm = s * np.sqrt(g_pm / n)
            se_pb = pyblock_se(x)

            rows.append({
                "phi": phi, "rep": r, "n": n, "g_true": g_true(phi),
                "g_mdcheck": g_md, "g_pymbar": g_pm,
                "cover_naive": abs(m) <= Z95 * se_naive,
                "cover_mdcheck_bootstrap": lo_bs <= 0.0 <= hi_bs,
                "cover_mdcheck_ci": lo_ci <= 0.0 <= hi_ci,
                "cover_g_mdcheck": abs(m) <= Z95 * se_g_md,
                "cover_g_pymbar": abs(m) <= Z95 * se_g_pm,
                "cover_pyblock": (abs(m) <= Z95 * se_pb) if np.isfinite(se_pb) else np.nan,
                "halfwidth_bootstrap_over_true": (hi_bs - lo_bs) / 2.0 / (Z95 * np.sqrt(g_true(phi) / n)),
            })
    return pd.DataFrame(rows)


def bench_equilibration(phi, n, amp, tau_relax, reps, rng):
    rows = []
    t = np.arange(n)
    transient = amp * np.exp(-t / tau_relax)
    for r in range(reps):
        x = ar1(phi, n, rng) + transient
        md = detect_equilibration(x)
        t0_pm, g_pm, neff_pm = timeseries.detect_equilibration(x)
        rows.append({
            "phi": phi, "n": n, "amp": amp, "tau_relax": tau_relax, "rep": r,
            "t_eq_mdcheck": md["t_eq_index"], "t_eq_pymbar": int(t0_pm),
            "neff_mdcheck": md["n_eff"], "neff_pymbar": float(neff_pm),
            "bias_no_discard": float(x.mean()),
            "bias_mdcheck": float(x[md["t_eq_index"]:].mean()),
            "bias_pymbar": float(x[int(t0_pm):].mean()),
        })
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=200)
    ap.add_argument("--n", type=int, default=20000)
    ap.add_argument("--seed", type=int, default=20260924)
    ap.add_argument("--out", default=os.path.join(os.path.dirname(__file__), "results"))
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    t0 = time.time()
    phis = [0.0, 0.5, 0.8, 0.9, 0.95, 0.99]
    df = bench_inefficiency_and_coverage(phis, args.n, args.reps, rng)
    df.to_csv(os.path.join(args.out, "ar1_raw.csv"), index=False)

    summ = df.groupby("phi").agg(
        g_true=("g_true", "first"),
        g_mdcheck_mean=("g_mdcheck", "mean"), g_mdcheck_sd=("g_mdcheck", "std"),
        g_pymbar_mean=("g_pymbar", "mean"), g_pymbar_sd=("g_pymbar", "std"),
        cov_naive=("cover_naive", "mean"),
        cov_mdcheck_bootstrap_v100=("cover_mdcheck_bootstrap", "mean"),
        cov_mdcheck_ci=("cover_mdcheck_ci", "mean"),
        cov_g_mdcheck=("cover_g_mdcheck", "mean"),
        cov_g_pymbar=("cover_g_pymbar", "mean"),
        cov_pyblock=("cover_pyblock", "mean"),
        bootstrap_halfwidth_ratio=("halfwidth_bootstrap_over_true", "mean"),
    ).reset_index()
    summ.to_csv(os.path.join(args.out, "ar1_inefficiency_coverage.csv"), index=False)
    print(summ.to_string(index=False, float_format=lambda v: f"{v:.3f}"))

    eq_reps = max(20, args.reps // 4)
    eq = pd.concat([
        bench_equilibration(0.9, 10000, 5.0, 300.0, eq_reps, rng),
        bench_equilibration(0.95, 10000, 5.0, 1000.0, eq_reps, rng),
    ])
    eq.to_csv(os.path.join(args.out, "equilibration_raw.csv"), index=False)
    eq_s = eq.groupby(["phi", "tau_relax"]).agg(
        t_eq_mdcheck=("t_eq_mdcheck", "median"), t_eq_pymbar=("t_eq_pymbar", "median"),
        abs_diff_t_eq=("t_eq_mdcheck", lambda s: np.nan),
        bias_no_discard=("bias_no_discard", "mean"),
        bias_mdcheck=("bias_mdcheck", "mean"), bias_pymbar=("bias_pymbar", "mean"),
    ).reset_index()
    eq_s["median_abs_diff_t_eq"] = eq.assign(d=(eq.t_eq_mdcheck - eq.t_eq_pymbar).abs()) \
        .groupby(["phi", "tau_relax"])["d"].median().values
    eq_s = eq_s.drop(columns=["abs_diff_t_eq"])
    eq_s.to_csv(os.path.join(args.out, "equilibration_summary.csv"), index=False)
    print()
    print(eq_s.to_string(index=False, float_format=lambda v: f"{v:.3f}"))

    meta = {
        "mdcheck": mdcheck.__version__, "pymbar": pymbar.__version__,
        "pyblock": getattr(pyblock, "__version__", "unknown"), "numpy": np.__version__,
        "python": platform.python_version(), "seed": args.seed, "reps": args.reps, "n": args.n,
        "runtime_s": round(time.time() - t0, 1),
    }
    with open(os.path.join(args.out, "versions.json"), "w") as fh:
        json.dump(meta, fh, indent=2)
    print("\n", meta)


if __name__ == "__main__":
    main()
