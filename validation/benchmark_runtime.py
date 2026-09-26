"""
Runtime and agreement of equilibration detection: MDCheck 1.1 (hierarchical cutoff search),
MDCheck 1.0 (exhaustive stride-10 grid), pymbar.timeseries.detect_equilibration (exact and
fast=True with nskip = N // 200). AR(1), phi = 0.9, with an exponential initial transient.

    python validation/benchmark_runtime.py [--out validation/results]
"""
import argparse
import os
import platform
import sys
import time

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from mdcheck.core.equilibration import detect_equilibration  # noqa: E402
from pymbar import timeseries  # noqa: E402


def series(n, rng):
    x = np.empty(n)
    x[0] = 0.0
    e = rng.normal(0.0, np.sqrt(1 - 0.81), n)
    for i in range(1, n):
        x[i] = 0.9 * x[i - 1] + e[i]
    return x + 3.0 * np.exp(-np.arange(n) / (n / 20))


def timed(fn):
    t = time.perf_counter()
    out = fn()
    return time.perf_counter() - t, out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.join(os.path.dirname(__file__), "results"))
    args = ap.parse_args()
    rng = np.random.default_rng(0)
    rows = []
    for n in (1_000, 10_000, 100_000, 1_000_000):
        x = series(n, rng)
        row = {"n": n}
        row["mdcheck_s"], md = timed(lambda: detect_equilibration(x))
        row["mdcheck_t0"] = md["t_eq_index"]
        if n <= 10_000:
            row["mdcheck10_grid_s"], mo = timed(lambda: detect_equilibration(x, step_search=10))
            row["mdcheck10_grid_t0"] = mo["t_eq_index"]
            row["pymbar_s"], pm = timed(lambda: timeseries.detect_equilibration(x))
            row["pymbar_t0"] = int(pm[0])
        row["pymbar_fast_s"], pf = timed(lambda: timeseries.detect_equilibration(x, fast=True, nskip=max(1, n // 200)))
        row["pymbar_fast_t0"] = int(pf[0])
        rows.append(row)
        print(row, flush=True)
    df = pd.DataFrame(rows)
    df["machine"] = f"{platform.processor()} / Python {platform.python_version()}"
    df.to_csv(os.path.join(args.out, "runtime_benchmark.csv"), index=False)


if __name__ == "__main__":
    main()
