"""
Coverage of the MDCheck confidence interval as a function of series length N and of the
Madras-Sokal window constant c, on AR(1) processes (exact g, mean 0).

Intervals compared on the same series:
  naive      mean +/- 1.96 s / sqrt(N)
  z_g        mean +/- 1.96 s sqrt(g/N)                      (normal quantile)
  t_g        mean +/- t_{N/g - 1} s sqrt(g/N)               (MDCheck 1.1.0)
with g estimated with c in {4, 6, 8} (default 6).

    python validation/validate_ar1_length.py [--reps 500]
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import t as student_t

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from mdcheck.core.autocorrelation import integrated_autocorrelation_time  # noqa: E402


def ar1(phi, n, rng):
    x = np.empty(n)
    x[0] = rng.normal()
    e = rng.normal(0.0, np.sqrt(1.0 - phi * phi), n)
    for i in range(1, n):
        x[i] = phi * x[i - 1] + e[i]
    return x


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=500)
    ap.add_argument("--seed", type=int, default=2718)
    ap.add_argument("--out", default=os.path.join(os.path.dirname(__file__), "results"))
    args = ap.parse_args()
    rng = np.random.default_rng(args.seed)
    rows = []
    for phi in (0.5, 0.9, 0.95):
        g_true = (1 + phi) / (1 - phi)
        for n in (250, 500, 1000, 2000, 5000, 20000):
            for r in range(args.reps):
                x = ar1(phi, n, rng)
                m, s = x.mean(), x.std(ddof=1)
                row = {"phi": phi, "g_true": g_true, "n": n, "neff_true": n / g_true, "rep": r,
                       "naive": abs(m) <= 1.96 * s / np.sqrt(n)}
                for c in (4.0, 6.0, 8.0):
                    _, g, _ = integrated_autocorrelation_time(x, c_window=c)
                    se = s * np.sqrt(g / n)
                    row[f"g_c{c:g}"] = g
                    row[f"z_c{c:g}"] = abs(m) <= 1.96 * se
                    row[f"t_c{c:g}"] = abs(m) <= student_t.ppf(0.975, max(1.0, n / g - 1.0)) * se
                rows.append(row)
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(args.out, "ar1_length_raw.csv"), index=False)
    agg = {k: (k, "mean") for k in ["naive", "z_c6", "t_c4", "t_c6", "t_c8"]}
    agg.update({"g_c6_median": ("g_c6", "median"), "neff_true": ("neff_true", "first")})
    summ = df.groupby(["phi", "n"]).agg(**agg).reset_index()
    summ.to_csv(os.path.join(args.out, "ar1_length_coverage.csv"), index=False)
    pd.set_option("display.width", 200)
    print(summ.to_string(index=False, float_format=lambda v: f"{v:.3f}"))


if __name__ == "__main__":
    main()
