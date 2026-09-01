"""
MDCheck: Automated Convergence, Statistical Inefficiency, and Reproducibility
Assessment for Molecular Dynamics Simulations.
"""

__version__ = "1.0.0"
__author__ = "Andres Monreal-Hernández"
__license__ = "MIT"

from mdcheck.core.autocorrelation import (
    compute_autocorrelation,
    integrated_autocorrelation_time,
    effective_sample_size,
    statistical_inefficiency
)
from mdcheck.core.equilibration import detect_equilibration, geweke_diagnostic
from mdcheck.core.replicas import (
    assess_replica_consistency,
    jensen_shannon_distance,
    subspace_overlap_rmsip
)
from mdcheck.core.drift import assess_drift, block_averaging
from mdcheck.core.bootstrap import bootstrap_ci
from mdcheck.core.scoring import assess_trajectory_quality, SimulationQualityReport

__all__ = [
    "__version__",
    "compute_autocorrelation",
    "integrated_autocorrelation_time",
    "effective_sample_size",
    "statistical_inefficiency",
    "detect_equilibration",
    "geweke_diagnostic",
    "assess_replica_consistency",
    "jensen_shannon_distance",
    "subspace_overlap_rmsip",
    "assess_drift",
    "block_averaging",
    "bootstrap_ci",
    "assess_trajectory_quality",
    "SimulationQualityReport"
]
