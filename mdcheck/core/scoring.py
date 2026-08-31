"""
Simulation Quality Assessment, Multi-criteria Decision Matrix, and Scoring Engine.
"""

from typing import List, Dict, Any, Optional
from dataclasses import dataclass, asdict
import numpy as np

from mdcheck.core.autocorrelation import integrated_autocorrelation_time, effective_sample_size
from mdcheck.core.equilibration import detect_equilibration, geweke_diagnostic
from mdcheck.core.replicas import assess_replica_consistency
from mdcheck.core.drift import assess_drift, block_averaging
from mdcheck.core.bootstrap import bootstrap_ci


@dataclass
class ObservableReport:
    name: str
    n_total: int
    n_prod: int
    t_eq_time: float
    fraction_discarded: float
    tau_int: float
    g_inefficiency: float
    n_eff: float
    mean_prod: float
    ci_lower_95: float
    ci_upper_95: float
    std_prod: float
    drift_status: str
    relative_drift_pct: float
    geweke_pvalue: float
    status: str
    diagnostic_message: str


@dataclass
class SimulationQualityReport:
    overall_status: str
    score_summary: str
    observables: Dict[str, ObservableReport]
    replica_assessment: Optional[Dict[str, Any]]
    key_metrics: Dict[str, Any]
    recommendations: List[str]
    provenance: Dict[str, Any]

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


def assess_trajectory_quality(
    timeseries_dict: Dict[str, np.ndarray],
    time_coords: Optional[np.ndarray] = None,
    replica_dict: Optional[Dict[str, List[np.ndarray]]] = None,
    min_neff_pass: float = 500.0,
    min_neff_warn: float = 100.0,
    max_discard_warn: float = 0.50
) -> SimulationQualityReport:
    """
    Evaluates full molecular dynamics trajectory quality across all provided observables.

    Parameters
    ----------
    timeseries_dict : dict
        Mapping of observable name -> 1D numpy array (e.g. {'RMSD': arr, 'Rg': arr, 'TotalEnergy': arr}).
    time_coords : np.ndarray, optional
        Array of timestamps corresponding to frames.
    replica_dict : dict, optional
        Mapping of observable name -> list of 1D arrays for replicas.
    min_neff_pass : float, default 500.0
        Minimum effective sample size for PASS.
    min_neff_warn : float, default 100.0
        Minimum effective sample size for WARNING.
    max_discard_warn : float, default 0.50
        Maximum allowed equilibration fraction before WARNING.

    Returns
    -------
    report : SimulationQualityReport
        Consolidated quality report.
    """
    observable_reports: Dict[str, ObservableReport] = {}
    recommendations: List[str] = []
    
    statuses = []
    
    for obs_name, series in timeseries_dict.items():
        y = np.asarray(series, dtype=np.float64)
        n_total = len(y)
        
        # 1. Equilibration detection
        eq_res = detect_equilibration(y, time_series=time_coords)
        t_eq_idx = eq_res["t_eq_index"]
        t_eq_time = eq_res["t_eq_time"]
        fraction_discarded = eq_res["fraction_discarded"]
        
        # Production slice
        y_prod = y[t_eq_idx:]
        t_prod = time_coords[t_eq_idx:] if time_coords is not None else None
        
        # 2. Autocorrelation & N_eff on production
        tau_int, g, _ = integrated_autocorrelation_time(y_prod)
        n_eff = float(len(y_prod)) / max(1.0, g)
        
        # 3. Drift on production
        drift_res = assess_drift(y_prod, time_series=t_prod)
        
        # 4. Geweke stationarity
        geweke_res = geweke_diagnostic(y_prod)
        
        # 5. Bootstrap confidence intervals
        mean_val, ci_low, ci_high = bootstrap_ci(y_prod, stat_func=np.mean)
        std_val = float(np.std(y_prod, ddof=1)) if len(y_prod) > 1 else 0.0
        
        # 6. Scoring rules for this observable
        obs_status = "PASS"
        diag_messages = []
        
        if n_eff < min_neff_warn:
            obs_status = "FAIL"
            diag_messages.append(f"Insufficient effective samples (N_eff = {n_eff:.1f} < {min_neff_warn:.0f}).")
        elif n_eff < min_neff_pass:
            if obs_status != "FAIL":
                obs_status = "WARNING"
            diag_messages.append(f"Marginal effective samples (N_eff = {n_eff:.1f} < {min_neff_pass:.0f}).")
            
        if fraction_discarded > max_discard_warn:
            if obs_status != "FAIL":
                obs_status = "WARNING"
            diag_messages.append(f"Long equilibration phase ({fraction_discarded*100:.1f}% discarded).")
            
        if drift_res["status"] == "FAIL":
            obs_status = "FAIL"
            diag_messages.append(f"Significant production drift ({drift_res['relative_drift']*100:.2f}% relative change).")
        elif drift_res["status"] == "WARNING":
            if obs_status != "FAIL":
                obs_status = "WARNING"
            diag_messages.append(f"Mild production drift detected.")
            
        if not geweke_res["is_stationary"]:
            if obs_status != "FAIL":
                obs_status = "WARNING"
            diag_messages.append(f"Geweke diagnostic failed stationarity (p = {geweke_res['p_value']:.4f}).")
            
        if not diag_messages:
            diag_messages.append("Observable converged; stationary production phase verified.")
            
        statuses.append(obs_status)
        
        observable_reports[obs_name] = ObservableReport(
            name=obs_name,
            n_total=n_total,
            n_prod=len(y_prod),
            t_eq_time=t_eq_time,
            fraction_discarded=fraction_discarded,
            tau_int=tau_int,
            g_inefficiency=g,
            n_eff=n_eff,
            mean_prod=mean_val,
            ci_lower_95=ci_low,
            ci_upper_95=ci_high,
            std_prod=std_val,
            drift_status=drift_res["status"],
            relative_drift_pct=float(drift_res["relative_drift"] * 100.0),
            geweke_pvalue=geweke_res["p_value"],
            status=obs_status,
            diagnostic_message=" ".join(diag_messages)
        )
        
    # Multi-replica evaluation if available
    replica_assessment = None
    if replica_dict is not None and len(replica_dict) > 0:
        primary_key = list(replica_dict.keys())[0]
        replica_assessment = assess_replica_consistency(replica_dict[primary_key])
        statuses.append(replica_assessment["status"])
        if replica_assessment["status"] != "PASS":
            recommendations.append(replica_assessment["recommendation"])

    # Global status calculation
    if "FAIL" in statuses:
        overall_status = "FAIL"
        score_summary = "SIMULATION QUALITY = INSUFFICIENT / NON-CONVERGED"
    elif "WARNING" in statuses:
        overall_status = "WARNING"
        score_summary = "SIMULATION QUALITY = BORDERLINE / ACCEPTABLE WITH CAUTION"
    else:
        overall_status = "PASS"
        score_summary = "SIMULATION QUALITY = EXCELLENT / CERTIFIED"

    # Aggregated key metrics
    total_frames = max([o.n_total for o in observable_reports.values()]) if observable_reports else 0
    min_eff = min([o.n_eff for o in observable_reports.values()]) if observable_reports else 0
    max_teq = max([o.t_eq_time for o in observable_reports.values()]) if observable_reports else 0
    
    key_metrics = {
        "overall_status": overall_status,
        "total_frames": total_frames,
        "max_equilibration_time": max_teq,
        "min_effective_sample_size": min_eff,
        "n_observables_evaluated": len(observable_reports),
        "n_replicas_evaluated": replica_assessment["n_replicas"] if replica_assessment else 1
    }

    return SimulationQualityReport(
        overall_status=overall_status,
        score_summary=score_summary,
        observables=observable_reports,
        replica_assessment=replica_assessment,
        key_metrics=key_metrics,
        recommendations=recommendations,
        provenance={
            "tool": "MDCheck",
            "version": "1.0.0",
            "citation": "Monreal-Hernández, A. (2026). MDCheck: Automated Convergence, Statistical Inefficiency, and Reproducibility Assessment for Molecular Dynamics Simulations."
        }
    )
