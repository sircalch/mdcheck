"""
Tests for scoring, reporters, and CLI workflows.
"""

import os
import tempfile
import numpy as np
import pytest
from mdcheck.core.scoring import assess_trajectory_quality
from mdcheck.reporters.plot_generator import generate_publication_figures
from mdcheck.reporters.manuscript_prep import generate_manuscript_assets
from mdcheck.reporters.html_report import generate_html_report
from mdcheck.cli import run_demo


def test_full_scoring_and_reporting_pipeline():
    rng = np.random.default_rng(42)
    n = 1000
    t = np.linspace(0, 50, n)
    
    # Observable 1: Well-behaved RMSD
    eq_phase = np.linspace(1.0, 2.0, 200)
    prod_phase = rng.normal(2.0, 0.05, 800)
    rmsd = np.concatenate([eq_phase, prod_phase])
    
    # Observable 2: Rg
    rg = rng.normal(1.5, 0.02, 1000)
    
    timeseries_dict = {"RMSD": rmsd, "Rg": rg}
    
    report = assess_trajectory_quality(timeseries_dict, time_coords=t)
    assert report.overall_status in ["PASS", "WARNING"]
    assert "RMSD" in report.observables
    assert "Rg" in report.observables
    assert report.observables["RMSD"].n_eff > 100
    
    with tempfile.TemporaryDirectory() as tmpdir:
        # Test plot generation
        saved_plots = generate_publication_figures(timeseries_dict, t, report, tmpdir, formats=["png", "svg"])
        assert len(saved_plots) > 0
        for p in saved_plots:
            assert os.path.exists(p)
            
        # Test manuscript assets
        assets = generate_manuscript_assets(report, tmpdir)
        assert os.path.exists(assets["summary_csv"])
        assert os.path.exists(assets["summary_tex"])
        assert os.path.exists(assets["methods_text"])
        assert os.path.exists(assets["citation_bib"])
        
        # Test HTML report
        html_p = os.path.join(tmpdir, "report.html")
        generate_html_report(report, html_p, methods_text="Sample methods", citation_bib="@software{}")
        assert os.path.exists(html_p)
        assert os.path.getsize(html_p) > 500


def test_cli_demo_execution():
    with tempfile.TemporaryDirectory() as tmpdir:
        run_demo(output_dir=tmpdir)
        assert os.path.exists(os.path.join(tmpdir, "report.html"))
        assert os.path.exists(os.path.join(tmpdir, "mdcheck_summary_table.csv"))
        assert os.path.exists(os.path.join(tmpdir, "citation.bib"))
