"""
Worked example: three replicas of the radius of gyration, from a CSV file to a report.

Uses sample_data/multi_replica_gyrate.csv (time in ns plus one column per replica).

Run from the repository root:

    python examples/replicas_to_report.py
"""

import os

from mdcheck import assess_trajectory_quality
from mdcheck.parsers.generic import load_timeseries_file
from mdcheck.reporters import (
    generate_html_report,
    generate_manuscript_assets,
    generate_publication_figures,
)

DATA_FILE = os.path.join("sample_data", "multi_replica_gyrate.csv")
OUTPUT_DIR = "replicas_example_output"


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # 1. Read the file: one column per replica.
    time, columns, meta = load_timeseries_file(DATA_FILE)
    replicas = [columns[name] for name in columns]
    print(f"Replicas: {list(columns)}  ({len(time)} frames each)")

    # 2. Assess the first replica as the primary trajectory and compare all three.
    primary = {"Radius_of_Gyration": replicas[0]}
    report = assess_trajectory_quality(
        timeseries_dict=primary,
        time_coords=time,
        replica_dict={"Radius_of_Gyration": replicas},
    )

    obs = report.observables["Radius_of_Gyration"]
    print(f"\nRadius of gyration, replica 1: t_eq = {obs.t_eq_time:.2f} ns, N_eff = {obs.n_eff:.0f}, "
          f"mean = {obs.mean_prod:.3f} [{obs.ci_lower_95:.3f}, {obs.ci_upper_95:.3f}], status {obs.status}")

    rep = report.replica_assessment
    print(f"Replica comparison: status {rep['status']}, mean JSD = {rep['mean_jsd']:.3f}")
    print(f"Overall status: {report.overall_status}")

    # 3. Write the report, figures, tables and methods paragraph.
    generate_publication_figures(primary, time, report, OUTPUT_DIR,
                                 replica_dict={"Radius_of_Gyration": replicas})
    generate_manuscript_assets(report, OUTPUT_DIR)
    generate_html_report(report, os.path.join(OUTPUT_DIR, "report.html"))
    print(f"\nReport written to {os.path.join(OUTPUT_DIR, 'report.html')}")


if __name__ == "__main__":
    main()
