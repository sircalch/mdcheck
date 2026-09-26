"""
Replica-consistency tests on the OpenMM trajectories: MDCheck's Cochran Q test versus the
Gelman-Rubin potential scale reduction factor (rank-normalised split R-hat, ArviZ) and the
Jensen-Shannon criterion of MDCheck 1.0.0.

Data: validation/results/md_{lj,water}_series.npz and validation/results_water300/md_water_series.npz.
Production segments start at the MDCheck equilibration time of each replica.

1. False-alarm rate: random triples of replicas of the same ensemble (no shift).
2. Power: a shift of k standard errors of the replica mean (k = 1..6) is added to one replica of
   each triple, using the replica's own autocorrelation-corrected SE, so the noise is real MD
   noise and only the mean differs.

Alarm definitions: Q test status != PASS (p < 0.05); R-hat > 1.01 (Vehtari et al. 2021);
R-hat > 1.1 (Gelman and Rubin 1992 convention); max JSD > 0.15 (MDCheck 1.0.0).

    python validation/compare_replica_tests.py [--triples 400]
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
import arviz as az

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from mdcheck.core.autocorrelation import integrated_autocorrelation_time  # noqa: E402
from mdcheck.core.equilibration import detect_equilibration  # noqa: E402
from mdcheck.core.replicas import assess_replica_consistency  # noqa: E402

HERE = os.path.dirname(__file__)
DATASETS = [
    ("LJ U (150 ps)", os.path.join(HERE, "results", "md_lj_series.npz"), "potential_energy_kJmol"),
    ("TIP3P rho (100 ps)", os.path.join(HERE, "results", "md_water_series.npz"), "density_gcm3"),
    ("TIP3P U (100 ps)", os.path.join(HERE, "results", "md_water_series.npz"), "potential_energy_kJmol"),
    ("TIP3P rho (300 ps)", os.path.join(HERE, "results_water300", "md_water_series.npz"), "density_gcm3"),
    ("TIP3P U (300 ps)", os.path.join(HERE, "results_water300", "md_water_series.npz"), "potential_energy_kJmol"),
]


def productions(npz_path, obs):
    d = np.load(npz_path)
    prods, ses = [], []
    for x in d["short_" + obs]:
        t0 = detect_equilibration(x)["t_eq_index"]
        p = x[t0:]
        _, g, _ = integrated_autocorrelation_time(p)
        prods.append(p)
        ses.append(float(np.std(p, ddof=1) * np.sqrt(g / len(p))))
    return prods, np.array(ses)


def alarms(triple):
    res = assess_replica_consistency(triple)
    m = min(len(t) for t in triple)
    arr = np.array([t[-m:] for t in triple])
    rhat = float(az.rhat(arr))
    return {"q_alarm": res["status"] != "PASS", "q_fail": res["status"] == "FAIL",
            "rhat": rhat, "rhat_101": rhat > 1.01, "rhat_11": rhat > 1.1,
            "jsd_alarm": res["max_jsd"] > 0.15}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--triples", type=int, default=400)
    ap.add_argument("--seed", type=int, default=31)
    ap.add_argument("--out", default=os.path.join(HERE, "results"))
    args = ap.parse_args()
    rng = np.random.default_rng(args.seed)
    rows = []
    for label, path, obs in DATASETS:
        prods, ses = productions(path, obs)
        n = len(prods)
        for _ in range(args.triples):
            idx = rng.choice(n, 3, replace=False)
            for k in (0, 1, 2, 3, 4, 6):
                trip = [prods[idx[0]], prods[idx[1]], prods[idx[2]] + k * ses[idx[2]]]
                rows.append({"dataset": label, "shift_se": k, **alarms(trip)})
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(args.out, "replica_tests_raw.csv"), index=False)
    summ = df.groupby(["dataset", "shift_se"]).agg(
        n=("q_alarm", "size"), q_test=("q_alarm", "mean"), q_fail=("q_fail", "mean"),
        rhat_gt_1_01=("rhat_101", "mean"), rhat_gt_1_1=("rhat_11", "mean"),
        jsd_gt_0_15=("jsd_alarm", "mean"), rhat_median=("rhat", "median")).reset_index()
    summ.to_csv(os.path.join(args.out, "replica_tests_summary.csv"), index=False)
    pd.set_option("display.width", 200)
    print(summ.to_string(index=False, float_format=lambda v: f"{v:.3f}"))


if __name__ == "__main__":
    main()
