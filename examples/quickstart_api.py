"""
Quickstart API Example for MDCheck.
"""

import os
import numpy as np
from mdcheck import assess_trajectory_quality
from mdcheck.reporters import (
    generate_publication_figures,
    generate_manuscript_assets,
    generate_html_report
)


def main():
    print("Running MDCheck Python API quickstart example...")
    output_dir = "quickstart_output"
    os.makedirs(output_dir, exist_ok=True)
    
    # 1. Synthetic 50 ns trajectory data
    n_points = 1000
    t = np.linspace(0, 50, n_points)
    
    # Backbone RMSD (Angstroms)
    rmsd = 1.0 * (1.0 - np.exp(-t / 5.0)) + 1.2 + np.random.normal(0, 0.08, n_points)
    # Radius of Gyration (nm)
    rg = 1.50 + np.random.normal(0, 0.02, n_points)
    
    timeseries = {
        "Backbone_RMSD": rmsd,
        "Radius_of_Gyration": rg
    }
    
    # 2. Run Assessment
    report = assess_trajectory_quality(timeseries_dict=timeseries, time_coords=t)
    
    print(f"\nOverall Simulation Status: {report.overall_status}")
    print(f"Summary: {report.score_summary}")
    for k, v in report.observables.items():
        print(f" - {k}: t_eq = {v.t_eq_time:.2f} ns | N_eff = {v.n_eff:.0f} | Status: {v.status}")
        
    # 3. Export all publication assets
    generate_publication_figures(timeseries, t, report, output_dir)
    assets = generate_manuscript_assets(report, output_dir)
    
    with open(assets["methods_text"], "r", encoding="utf-8") as f:
        methods_txt = f.read()
    with open(assets["citation_bib"], "r", encoding="utf-8") as f:
        bib_txt = f.read()
        
    html_p = os.path.join(output_dir, "report.html")
    generate_html_report(report, html_p, methods_text=methods_txt, citation_bib=bib_txt)
    
    print(f"\nCompleted! Check out: {os.path.abspath(html_p)}")


if __name__ == "__main__":
    main()
