"""
Command Line Interface (CLI) for MDCheck.
"""

import sys
import os
import argparse
import numpy as np

from mdcheck import __version__
from mdcheck.parsers.generic import load_timeseries_file
from mdcheck.core.scoring import assess_trajectory_quality
from mdcheck.reporters.plot_generator import generate_publication_figures
from mdcheck.reporters.manuscript_prep import generate_manuscript_assets
from mdcheck.reporters.html_report import generate_html_report


def print_banner():
    banner = rf"""
  __  __ _____   _____ _               _    
 |  \/  |  __ \ / ____| |             | |   
 | \  / | |  | | |    | |__   ___  ___| | __
 | |\/| | |  | | |    | '_ \ / _ \/ __| |/ /
 | |  | | |__| | |____| | | |  __/ (__|   < 
 |_|  |_|_____/ \_____|_| |_|\___|\___|_|\_\ v{__version__}

 Molecular Dynamics Convergence & Statistical Inefficiency Toolkit
 Monreal-Hernández et al., 2026
"""
    print(banner)


def run_demo(output_dir: str = "mdcheck_demo_output"):
    """
    Generates a realistic 3-replica molecular dynamics dataset (RMSD, Radius of Gyration, Potential Energy)
    and executes the full MDCheck assessment pipeline.
    """
    print(f"\n[MDCheck] Running demonstration mode...")
    os.makedirs(output_dir, exist_ok=True)
    
    n_frames = 2000
    time_coords = np.linspace(0.0, 100.0, n_frames)  # 100 ns simulation
    
    rng = np.random.default_rng(42)
    
    # Simulate Replica 1, 2, 3
    # Observable 1: RMSD (Equilibrating from 0.8 to 1.8 A over first 20 ns, then stable fluctuations)
    def make_rmsd(seed):
        r = np.random.default_rng(seed)
        noise = np.zeros(n_frames)
        # AR(1) correlated noise
        phi = 0.85
        for i in range(1, n_frames):
            noise[i] = phi * noise[i-1] + np.sqrt(1 - phi**2) * r.normal(0, 0.08)
        # Equilibration rise in first 400 frames (20 ns)
        eq_curve = 1.0 * (1.0 - np.exp(-time_coords / 8.0)) + 0.8
        return eq_curve + noise

    # Observable 2: Radius of Gyration Rg (stable around 1.45 nm)
    def make_rg(seed):
        r = np.random.default_rng(seed)
        noise = np.zeros(n_frames)
        phi = 0.80
        for i in range(1, n_frames):
            noise[i] = phi * noise[i-1] + np.sqrt(1 - phi**2) * r.normal(0, 0.02)
        eq_curve = 0.05 * np.exp(-time_coords / 5.0) + 1.45
        return eq_curve + noise

    r1_rmsd = make_rmsd(101)
    r2_rmsd = make_rmsd(102)
    r3_rmsd = make_rmsd(103)
    
    r1_rg = make_rg(201)
    r2_rg = make_rg(202)
    r3_rg = make_rg(203)
    
    primary_timeseries = {
        "Backbone_RMSD": r1_rmsd,
        "Radius_of_Gyration": r1_rg
    }
    
    replica_dict = {
        "Backbone_RMSD": [r1_rmsd, r2_rmsd, r3_rmsd],
        "Radius_of_Gyration": [r1_rg, r2_rg, r3_rg]
    }
    
    print("  -> Evaluating trajectory equilibration, autocorrelation, and replica overlap...")
    report = assess_trajectory_quality(
        timeseries_dict=primary_timeseries,
        time_coords=time_coords,
        replica_dict=replica_dict
    )
    
    # Generate assets
    print("  -> Generating publication-ready vector figures (SVG/PDF/PNG)...")
    generate_publication_figures(primary_timeseries, time_coords, report, output_dir, replica_dict=replica_dict)
    
    print("  -> Formulating Methods text snippet, summary LaTeX tables, and BibTeX citations...")
    text_assets = generate_manuscript_assets(report, output_dir)
    
    with open(text_assets["methods_text"], "r", encoding="utf-8") as f:
        methods_content = f.read()
    with open(text_assets["citation_bib"], "r", encoding="utf-8") as f:
        bib_content = f.read()
        
    html_path = os.path.join(output_dir, "report.html")
    print(f"  -> Rendering interactive HTML report at {html_path}...")
    generate_html_report(report, html_path, methods_text=methods_content, citation_bib=bib_content)
    
    print("\n" + "="*70)
    print(f" [RESULT] Overall Trajectory Certification: {report.overall_status}")
    print(f" [SCORE]  {report.score_summary}")
    print("="*70)
    for name, obs in report.observables.items():
        print(f" * {name:20s}: t_eq = {obs.t_eq_time:5.2f} ns | N_eff = {obs.n_eff:6.0f} | Mean = {obs.mean_prod:.3f} [{obs.ci_lower_95:.3f}, {obs.ci_upper_95:.3f}] | Status: {obs.status}")
    if report.replica_assessment:
        print(f" * Multi-Replica Overlap : Mean JSD = {report.replica_assessment['mean_jsd']:.3f} | Status: {report.replica_assessment['status']}")
    print("="*70)
    print(f"\nAll outputs successfully saved to: {os.path.abspath(output_dir)}/")
    print(f"Open {os.path.abspath(html_path)} in your browser to inspect the full report.\n")


def run_assess(args):
    """
    Executes assessment on user-provided input files.
    """
    output_dir = args.output
    os.makedirs(output_dir, exist_ok=True)
    
    input_files = args.input
    replica_files = args.replicas if args.replicas else []
    
    if not input_files and not replica_files:
        print("[Error] Please specify input files with --input or --replicas.", file=sys.stderr)
        sys.exit(1)
        
    all_files = input_files + replica_files
    print(f"\n[MDCheck] Loading timeseries from {len(all_files)} file(s)...")
    
    loaded_data = []
    for fp in all_files:
        t, s_dict, meta = load_timeseries_file(fp)
        loaded_data.append((t, s_dict, meta, fp))
        print(f"  -> Loaded '{meta['title']}': {len(t)} frames, observables: {list(s_dict.keys())}")
        
    primary_t, primary_dict, primary_meta, _ = loaded_data[0]
    
    replica_dict = None
    if len(loaded_data) > 1:
        replica_dict = {}
        for obs_name in primary_dict.keys():
            rep_list = []
            for _, s_dict, _, _ in loaded_data:
                if obs_name in s_dict:
                    rep_list.append(s_dict[obs_name])
            if len(rep_list) > 1:
                replica_dict[obs_name] = rep_list

    print("\n[MDCheck] Performing statistical convergence, equilibration, and inefficiency analysis...")
    report = assess_trajectory_quality(
        timeseries_dict=primary_dict,
        time_coords=primary_t,
        replica_dict=replica_dict
    )
    
    print("  -> Generating publication-quality vector charts...")
    generate_publication_figures(primary_dict, primary_t, report, output_dir, replica_dict=replica_dict)
    
    print("  -> Generating manuscript text, LaTeX summary table, and BibTeX citations...")
    text_assets = generate_manuscript_assets(report, output_dir)
    
    with open(text_assets["methods_text"], "r", encoding="utf-8") as f:
        methods_content = f.read()
    with open(text_assets["citation_bib"], "r", encoding="utf-8") as f:
        bib_content = f.read()
        
    html_path = os.path.join(output_dir, "report.html")
    print(f"  -> Writing HTML quality report to {html_path}...")
    generate_html_report(report, html_path, methods_text=methods_content, citation_bib=bib_content)
    
    print("\n" + "="*70)
    print(f" [RESULT] Overall Trajectory Certification: {report.overall_status}")
    print(f" [SCORE]  {report.score_summary}")
    print("="*70)
    for name, obs in report.observables.items():
        print(f" * {name:20s}: t_eq = {obs.t_eq_time:5.2f} | N_eff = {obs.n_eff:6.0f} | Mean = {obs.mean_prod:.3f} [{obs.ci_lower_95:.3f}, {obs.ci_upper_95:.3f}] | Status: {obs.status}")
    if report.replica_assessment:
        print(f" * Multi-Replica Overlap : Mean JSD = {report.replica_assessment['mean_jsd']:.3f} | Status: {report.replica_assessment['status']}")
    print("="*70)
    print(f"\nReport ready at: {os.path.abspath(html_path)}\n")


def print_citation():
    bib = """@software{monreal2026mdcheck,
  author = {Monreal-Hern\\'andez, Andre},
  title = {{MDCheck: Automated Convergence, Statistical Inefficiency, and Reproducibility Assessment for Molecular Dynamics Simulations}},
  year = {2026},
  version = {1.0.0},
  publisher = {Zenodo},
  url = {https://github.com/sircalch/mdcheck}
}"""
    print("\nIf you use MDCheck in your publications, please cite:\n")
    print("APA Style:")
    print("Monreal-Hernández, A. (2026). MDCheck: Automated Convergence, Statistical Inefficiency, and Reproducibility Assessment for Molecular Dynamics Simulations (v1.0.0). Zenodo. https://github.com/sircalch/mdcheck\n")
    print("BibTeX:")
    print(bib)
    print()


def main():
    parser = argparse.ArgumentParser(
        prog="mdcheck",
        description="MDCheck: Automated Convergence, Statistical Inefficiency, and Reproducibility Assessment for Molecular Dynamics Simulations."
    )
    parser.add_argument("-v", "--version", action="version", version=f"mdcheck {__version__}")
    
    subparsers = parser.add_subparsers(dest="command", help="Available subcommands")
    
    # Assess command
    assess_parser = subparsers.add_parser("assess", help="Assess trajectory convergence, equilibration, and reproducibility")
    assess_parser.add_argument("-i", "--input", nargs="+", help="Primary trajectory timeseries file(s) (.xvg, .csv, .dat, .txt)")
    assess_parser.add_argument("-r", "--replicas", nargs="+", help="Multi-replica files for consistency comparison")
    assess_parser.add_argument("-o", "--output", default="mdcheck_output", help="Directory to save report, plots, and manuscript assets (default: mdcheck_output)")
    
    # Demo command
    demo_parser = subparsers.add_parser("demo", help="Run MDCheck on a synthetic multi-replica simulation dataset")
    demo_parser.add_argument("-o", "--output", default="mdcheck_demo_output", help="Output directory (default: mdcheck_demo_output)")
    
    # Cite command
    subparsers.add_parser("cite", help="Display BibTeX and APA citation details")
    
    if len(sys.argv) == 1:
        print_banner()
        parser.print_help()
        sys.exit(0)
        
    args = parser.parse_args()
    
    if args.command == "assess":
        print_banner()
        run_assess(args)
    elif args.command == "demo":
        print_banner()
        run_demo(args.output)
    elif args.command == "cite":
        print_banner()
        print_citation()
    else:
        parser.print_help()


if __name__ == "__main__":
    main()

