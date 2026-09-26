"""
Compares statistical-inefficiency estimators on identical data: CI coverage of the reference mean
on the OpenMM trajectories (validation/results/md_*_series.npz) and on AR(1) processes.

Estimators (all use the same ACF and the same MDCheck equilibration cutoff):
  sokal_c6        tau = 1/2 + sum rho(k), stop at first rho < 0 or k >= 6 tau   (MDCheck 1.0.0)
  first_negative  g = 1 + 2 sum rho(k) up to the first rho(k) <= 0             (no window)
  geyer_imse      Geyer (1992) initial monotone positive sequence of rho(2k)+rho(2k+1)
  pymbar          pymbar.timeseries.statistical_inefficiency

    python validation/compare_estimators.py [--results validation/results]
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import t as student_t

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from mdcheck.core.autocorrelation import compute_autocorrelation  # noqa: E402
from mdcheck.core.equilibration import detect_equilibration  # noqa: E402
from pymbar import timeseries  # noqa: E402


def g_sokal(rho, c=6.0):
    tau = 0.5
    for k in range(1, len(rho)):
        if rho[k] < 0:
            break
        tau += rho[k]
        if k >= c * tau:
            break
    return max(1.0, 2.0 * tau)


def g_first_negative(rho):
    s = 0.0
    for k in range(1, len(rho)):
        if rho[k] <= 0:
            break
        s += rho[k]
    return max(1.0, 1.0 + 2.0 * s)


def g_geyer(rho):
    n_pairs = (len(rho) - 1) // 2
    gam = rho[0:2 * n_pairs:2] + rho[1:2 * n_pairs + 1:2]
    total, prev = 0.0, np.inf
    for gk in gam:
        if gk <= 0:
            break
        gk = min(gk, prev)
        total += gk
        prev = gk
    return max(1.0, 2.0 * total - 1.0)


ESTIMATORS = {
    "sokal_c6": lambda x, rho: g_sokal(rho),
    "first_negative": lambda x, rho: g_first_negative(rho),
    "geyer_imse": lambda x, rho: g_geyer(rho),
    "pymbar": lambda x, rho: float(timeseries.statistical_inefficiency(x)),
}


def covers(prod, g, ref):
    n = len(prod)
    se = prod.std(ddof=1) * np.sqrt(g / n)
    tc = student_t.ppf(0.975, max(1.0, n / g - 1.0))
    return abs(prod.mean() - ref) <= tc * se


def eval_md(npz_path):
    d = np.load(npz_path)
    disc = int(d["ref_discard_samples"])
    rows = []
    for key in [k for k in d.files if k.startswith("short_")]:
        obs = key[len("short_"):]
        ref = float(np.mean([run[disc:].mean() for run in d["long_" + obs]]))
        for r, x in enumerate(d[key]):
            t0 = detect_equilibration(x)["t_eq_index"]
            prod = x[t0:]
            rho = compute_autocorrelation(prod)
            for name, fn in ESTIMATORS.items():
                g = fn(prod, rho)
                rows.append({"data": os.path.basename(npz_path).split("_")[1], "observable": obs,
                             "replica": r, "estimator": name, "g": g, "covers": covers(prod, g, ref)})
    return pd.DataFrame(rows)


def eval_ar1(phis=(0.5, 0.9, 0.99), n=20000, reps=300, seed=7):
    rng = np.random.default_rng(seed)
    rows = []
    for phi in phis:
        for r in range(reps):
            x = np.empty(n)
            x[0] = rng.normal()
            e = rng.normal(0, np.sqrt(1 - phi * phi), n)
            for i in range(1, n):
                x[i] = phi * x[i - 1] + e[i]
            rho = compute_autocorrelation(x)
            for name, fn in ESTIMATORS.items():
                g = fn(x, rho)
                rows.append({"data": f"ar1_phi{phi}", "observable": "x", "replica": r,
                             "estimator": name, "g": g, "covers": covers(x, g, 0.0)})
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", default=os.path.join(os.path.dirname(__file__), "results"))
    ap.add_argument("--ar1-reps", type=int, default=300)
    args = ap.parse_args()
    frames = [eval_ar1(reps=args.ar1_reps)]
    for s in ("lj", "water"):
        p = os.path.join(args.results, f"md_{s}_series.npz")
        if os.path.exists(p):
            frames.append(eval_md(p))
    df = pd.concat(frames, ignore_index=True)
    df.to_csv(os.path.join(args.results, "estimator_comparison_raw.csv"), index=False)
    summ = df.groupby(["data", "observable", "estimator"]).agg(
        n=("covers", "size"), coverage=("covers", "mean"), g_median=("g", "median")).reset_index()
    summ.to_csv(os.path.join(args.results, "estimator_comparison.csv"), index=False)
    print(summ.pivot_table(index=["data", "observable"], columns="estimator", values="coverage")
          .to_string(float_format=lambda v: f"{v:.3f}"))
    print()
    print(summ.pivot_table(index=["data", "observable"], columns="estimator", values="g_median")
          .to_string(float_format=lambda v: f"{v:.2f}"))


if __name__ == "__main__":
    main()
