"""
Report, manuscript snippets, and publication-ready figure generators.
"""

from mdcheck.reporters.plot_generator import generate_publication_figures
from mdcheck.reporters.manuscript_prep import generate_manuscript_assets
from mdcheck.reporters.html_report import generate_html_report

__all__ = [
    "generate_publication_figures",
    "generate_manuscript_assets",
    "generate_html_report"
]
