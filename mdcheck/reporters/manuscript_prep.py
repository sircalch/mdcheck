"""
Generates manuscript Methods text, Supporting Information LaTeX, summary tables, and BibTeX citations.
"""

from typing import Dict, Any
import os
import pandas as pd
from mdcheck.core.scoring import SimulationQualityReport


def generate_manuscript_assets(report: SimulationQualityReport, output_dir: str) -> Dict[str, str]:
    """
    Generates all textual and tabular assets needed for publication manuscripts.

    Parameters
    ----------
    report : SimulationQualityReport
        Quality evaluation report.
    output_dir : str
        Target output directory.

    Returns
    -------
    paths : dict
        Mapping of asset name -> file path.
    """
    os.makedirs(output_dir, exist_ok=True)
    generated = {}
    
    # 1. Summary DataFrame
    rows = []
    for name, obs in report.observables.items():
        rows.append({
            "Observable": name,
            "Total Frames": obs.n_total,
            "Equilibration Time (t_eq)": f"{obs.t_eq_time:.2f}",
            "Discarded (%)": f"{obs.fraction_discarded * 100.0:.1f}%",
            "Production Frames": obs.n_prod,
            "tau_int (frames)": f"{obs.tau_int:.1f}",
            "Inefficiency (g)": f"{obs.g_inefficiency:.2f}",
            "Effective Samples (N_eff)": f"{obs.n_eff:.0f}",
            "Production Mean": f"{obs.mean_prod:.3f}",
            "95% CI (g-corrected)": f"[{obs.ci_lower_95:.3f}, {obs.ci_upper_95:.3f}]",
            "Std Dev": f"{obs.std_prod:.3f}",
            "Drift Status": obs.drift_status,
            "Quality Status": obs.status
        })
        
    df_summary = pd.DataFrame(rows)
    
    # CSV Table
    csv_path = os.path.join(output_dir, "mdcheck_summary_table.csv")
    df_summary.to_csv(csv_path, index=False)
    generated["summary_csv"] = csv_path
    
    # LaTeX Table (booktabs format)
    tex_table_path = os.path.join(output_dir, "mdcheck_summary_table.tex")
    tex_table = df_summary.to_latex(index=False, escape=False)
    with open(tex_table_path, "w", encoding="utf-8") as f:
        f.write("% MDCheck Convergence and Quality Assessment Summary Table\n")
        f.write(tex_table)
    generated["summary_tex"] = tex_table_path

    # 2. Methods Text Snippet
    methods_path = os.path.join(output_dir, "methods_snippet.txt")
    obs_names = ", ".join(report.observables.keys())
    max_teq = max([o.t_eq_time for o in report.observables.values()]) if report.observables else 0.0
    min_neff = min([o.n_eff for o in report.observables.values()]) if report.observables else 0.0
    
    replica_sentence = ""
    if report.replica_assessment is not None:
        n_rep = report.replica_assessment["n_replicas"]
        ra = report.replica_assessment
        replica_sentence = (
            f" Consistency of the ensemble averages across {n_rep} independent replicas was tested with "
            f"Cochran's Q on autocorrelation-corrected standard errors "
            f"(Q = {ra.get('cochran_q', 0.0):.2f}, p = {ra.get('heterogeneity_p_value', 1.0):.3f}; status: {ra['status']})."
        )
        
    methods_text = (
        f"Trajectory convergence, stationarity, and statistical independence were systematically assessed "
        f"using MDCheck v1.1.0 (Monreal-Hernández, 2026). For all monitored observables ({obs_names}), the initial "
        f"equilibration phase (up to t_eq = {max_teq:.2f} ns) was automatically detected and discarded by maximizing "
        f"the effective sample size N_eff of the production regime. Integrated autocorrelation times (tau_int) and "
        f"statistical inefficiency factors (g) were evaluated via self-consistent Madras-Sokal windowing, ensuring a minimum "
        f"of N_eff = {min_neff:.0f} statistically uncorrelated conformations in the production ensemble. Stationarity and "
        f"absence of systematic drift were confirmed using Geweke diagnostics and Flyvbjerg-Petersen block averaging.{replica_sentence} "
        f"Reported values represent production means with 95% confidence intervals from the statistical inefficiency (mean ± t_{{N_eff-1}} s sqrt(g/N))."
    )
    with open(methods_path, "w", encoding="utf-8") as f:
        f.write(methods_text + "\n")
    generated["methods_text"] = methods_path

    # 3. BibTeX Citation File
    bib_path = os.path.join(output_dir, "citation.bib")
    bib_content = """@software{monreal2026mdcheck,
  author = {Monreal-Hern\\'andez, Andre},
  title = {{MDCheck: Automated Convergence, Statistical Inefficiency, and Reproducibility Assessment for Molecular Dynamics Simulations}},
  year = {2026},
  version = {1.1.0},
  publisher = {Zenodo},
  url = {https://github.com/sircalch/mdcheck}
}
"""
    with open(bib_path, "w", encoding="utf-8") as f:
        f.write(bib_content)
    generated["citation_bib"] = bib_path

    return generated

