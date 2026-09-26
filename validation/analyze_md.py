"""
Final analysis of the OpenMM validation trajectories (validation/results/md_*_series.npz,
written by validate_md_openmm.py). Produces the MD tables of the manuscript.

For every short replica and observable:
  * MDCheck t_eq (max-N_eff) and pymbar.timeseries.detect_equilibration t_eq;
  * coverage of the reference mean by the MDCheck 95% CI, testing
    |m - ref| <= t * sqrt(SE_replica^2 + SE_ref^2), where SE_ref is the between-run standard error
    of the long reference runs;
  * the same test with the naive SE (s/sqrt(n)), without equilibration discard, and with a
    conservative extra discard of `--extra-discard-ps`.
Replica consistency: non-overlapping groups of 3 replicas of the same ensemble are passed to
assess_replica_consistency; the fraction not returning PASS is the false-alarm rate
(nominal 5% for WARNING or FAIL, 1% for FAIL).

    python validation/analyze_md.py [--results validation/results]
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import t as student_t

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from mdcheck.core.autocorrelation import integrated_autocorrelation_time  # noqa: E402
from mdcheck.core.equilibration import detect_equilibration  # noqa: E402
from mdcheck.core.replicas import assess_replica_consistency  # noqa: E402
from pymbar import timeseries  # noqa: E402


def se_g(x):
    _, g, _ = integrated_autocorrelation_time(x)
    return float(np.std(x, ddof=1) * np.sqrt(g / len(x))), float(g)


def covers(x, ref, se_ref, g=None, naive=False):
    n = len(x)
    if naive:
        se, dof = float(np.std(x, ddof=1) / np.sqrt(n)), n - 1
    else:
        se, g_est = se_g(x) if g is None else (float(np.std(x, ddof=1) * np.sqrt(g / n)), g)
        dof = max(1.0, n / g_est - 1.0)
    return bool(abs(x.mean() - ref) <= student_t.ppf(0.975, dof) * np.hypot(se, se_ref))


def analyze(npz_path, extra_ps):
    d = np.load(npz_path)
    system = os.path.basename(npz_path).split("_")[1]
    dt = float(d["sample_interval_ps"])
    disc = int(d["ref_discard_samples"])
    extra = int(round(extra_ps / dt))
    rows, rep_rows, ref_rows = [], [], []
    for key in [k for k in d.files if k.startswith("short_")]:
        obs = key[len("short_"):]
        long_means = np.array([run[disc:].mean() for run in d["long_" + obs]])
        ref = float(long_means.mean())
        se_ref = float(long_means.std(ddof=1) / np.sqrt(len(long_means)))
        ref_rows.append({"system": system, "observable": obs, "reference_mean": ref,
                         "reference_se": se_ref, "n_long_runs": len(long_means),
                         "long_run_ps": d["long_" + obs].shape[1] * dt, "discard_ps": disc * dt})
        prods = []
        for r, x in enumerate(d[key]):
            md = detect_equilibration(x)
            t0 = md["t_eq_index"]
            t0_pm, _, _ = timeseries.detect_equilibration(x)
            prod = x[t0:]
            prods.append(prod)
            m, g = prod.mean(), md["g"]
            se_prod = float(np.std(prod, ddof=1) * np.sqrt(g / len(prod)))
            se_full, g_full = se_g(x)
            rows.append({
                "system": system, "observable": obs, "replica": r, "n_samples": len(x),
                "bias_no_discard_over_se_prod": (x.mean() - ref) / se_prod,
                "halfwidth_ratio_no_discard": se_full / se_prod,
                "t_eq_ps_mdcheck": t0 * dt, "t_eq_ps_pymbar": int(t0_pm) * dt, "g_mdcheck": g,
                "n_eff": md["n_eff"], "bias": m - ref, "bias_over_se": (m - ref) / np.hypot(se_prod, se_ref),
                "cover_mdcheck": covers(prod, ref, se_ref, g=g),
                "cover_naive": covers(prod, ref, se_ref, naive=True),
                "cover_no_discard": covers(x, ref, se_ref),
                "cover_extra_discard": covers(x[min(len(x) - 20, t0 + extra):], ref, se_ref),
            })
        for gi in range(len(prods) // 3):
            res = assess_replica_consistency(prods[3 * gi:3 * gi + 3])
            rep_rows.append({"system": system, "observable": obs, "group": gi, "status": res["status"],
                             "p_value": res["heterogeneity_p_value"], "max_jsd": res["max_jsd"],
                             "jsd_flag": res["jsd_flag"]})
    return pd.DataFrame(rows), pd.DataFrame(rep_rows), pd.DataFrame(ref_rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", default=os.path.join(os.path.dirname(__file__), "results"))
    ap.add_argument("--extra-discard-ps", type=float, default=20.0)
    args = ap.parse_args()
    parts = [analyze(os.path.join(args.results, f"md_{s}_series.npz"), args.extra_discard_ps)
             for s in ("lj", "water") if os.path.exists(os.path.join(args.results, f"md_{s}_series.npz"))]
    df = pd.concat([p[0] for p in parts], ignore_index=True)
    rep = pd.concat([p[1] for p in parts], ignore_index=True)
    refs = pd.concat([p[2] for p in parts], ignore_index=True)
    df.to_csv(os.path.join(args.results, "md_final_raw.csv"), index=False)
    rep.to_csv(os.path.join(args.results, "md_final_replica_groups.csv"), index=False)
    refs.to_csv(os.path.join(args.results, "md_reference_means.csv"), index=False)

    summ = df.groupby(["system", "observable"]).agg(
        n_replicas=("replica", "count"),
        t_eq_ps_mdcheck=("t_eq_ps_mdcheck", "median"), t_eq_ps_pymbar=("t_eq_ps_pymbar", "median"),
        g_median=("g_mdcheck", "median"), n_eff_median=("n_eff", "median"),
        mean_bias_over_se=("bias_over_se", "mean"),
        coverage_mdcheck=("cover_mdcheck", "mean"), coverage_naive=("cover_naive", "mean"),
        coverage_no_discard=("cover_no_discard", "mean"),
        coverage_extra_discard=("cover_extra_discard", "mean"),
        bias_no_discard_over_se=("bias_no_discard_over_se_prod", "mean"),
        ci_width_ratio_no_discard=("halfwidth_ratio_no_discard", "median"),
    ).reset_index()
    teq = df.assign(d=(df.t_eq_ps_mdcheck - df.t_eq_ps_pymbar).abs()).groupby(["system", "observable"])["d"].median()
    summ["median_abs_diff_t_eq_ps"] = teq.values
    rs = rep.groupby(["system", "observable"]).agg(
        n_groups=("group", "count"),
        replica_alarm_rate=("status", lambda s: float((s != "PASS").mean())),
        replica_fail_rate=("status", lambda s: float((s == "FAIL").mean())),
        legacy_jsd_alarm_rate=("max_jsd", lambda s: float((s > 0.15).mean())),
    ).reset_index()
    summ = summ.merge(rs, on=["system", "observable"])
    summ.to_csv(os.path.join(args.results, "md_final_summary.csv"), index=False)
    pd.set_option("display.width", 250)
    print(refs.to_string(index=False))
    print()
    print(summ.T.to_string())


if __name__ == "__main__":
    main()
