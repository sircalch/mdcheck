"""
Publication-ready vector figure generation for MD trajectory convergence.
"""

from typing import Dict, List, Optional
import os
import numpy as np
import matplotlib.pyplot as plt
from mdcheck.core.autocorrelation import compute_autocorrelation
from mdcheck.core.drift import block_averaging
from mdcheck.core.scoring import SimulationQualityReport

# Set clean scientific plotting style
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.size': 11,
    'axes.labelsize': 12,
    'axes.titlesize': 13,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 10,
    'figure.titlesize': 14,
    'figure.dpi': 300,
    'lines.linewidth': 1.8,
    'grid.alpha': 0.3,
    'grid.linestyle': '--'
})


def generate_publication_figures(
    timeseries_dict: Dict[str, np.ndarray],
    time_coords: np.ndarray,
    report: SimulationQualityReport,
    output_dir: str,
    replica_dict: Optional[Dict[str, List[np.ndarray]]] = None,
    formats: List[str] = ("png", "svg", "pdf")
) -> List[str]:
    """
    Generates publication-quality figures illustrating equilibration, autocorrelation,
    replica overlap, and block averaging convergence.

    Parameters
    ----------
    timeseries_dict : dict
        Mapping of observable name -> timeseries array.
    time_coords : np.ndarray
        Array of timestamps (ns/ps).
    report : SimulationQualityReport
        Quality assessment report.
    output_dir : str
        Directory to save figures.
    replica_dict : dict, optional
        Multi-replica data dictionary.
    formats : list of str
        Image formats to save.

    Returns
    -------
    saved_paths : list of str
        List of generated file paths.
    """
    os.makedirs(output_dir, exist_ok=True)
    saved_files = []
    
    # 1. Main Convergence & Equilibration Multi-panel Plot
    fig, axes = plt.subplots(len(timeseries_dict), 2, figsize=(12, 3.2 * len(timeseries_dict)), squeeze=False)
    
    for row_idx, (obs_name, y) in enumerate(timeseries_dict.items()):
        obs_rep = report.observables[obs_name]
        t_eq = obs_rep.t_eq_time
        t_eq_idx = int(np.searchsorted(time_coords, t_eq)) if time_coords is not None else 0
        
        # Subplot Left: Timeseries with Equilibration Shading
        ax_ts = axes[row_idx, 0]
        t = time_coords if time_coords is not None else np.arange(len(y))
        
        if t_eq_idx > 0:
            ax_ts.axvspan(t[0], t_eq, color="#ff9999", alpha=0.3, label=f"Equilibration ({obs_rep.fraction_discarded*100:.1f}%)")
            ax_ts.plot(t[:t_eq_idx], y[:t_eq_idx], color="#cc3333", alpha=0.7, linestyle=":")
            
        ax_ts.plot(t[t_eq_idx:], y[t_eq_idx:], color="#0066cc", label="Production Phase")
        ax_ts.axhline(obs_rep.mean_prod, color="#003366", linestyle="--", label=f"Mean: {obs_rep.mean_prod:.3f}")
        ax_ts.fill_between(t[t_eq_idx:], obs_rep.ci_lower_95, obs_rep.ci_upper_95, color="#0066cc", alpha=0.15, label="95% Bootstrap CI")
        
        ax_ts.set_xlabel("Time (ns)")
        ax_ts.set_ylabel(f"{obs_name}")
        ax_ts.set_title(f"{obs_name} Trajectory & Equilibration Cutoff (t_eq = {t_eq:.2f} ns)")
        ax_ts.grid(True)
        ax_ts.legend(loc="upper right", frameon=True, fontsize=8)
        
        # Subplot Right: Autocorrelation Function of Production
        ax_acf = axes[row_idx, 1]
        y_prod = y[t_eq_idx:]
        acf = compute_autocorrelation(y_prod, max_lag=min(len(y_prod) // 2, 500))
        lags = np.arange(len(acf))
        
        ax_acf.plot(lags, acf, color="#2e7d32", label=f"ACF (tau_int = {obs_rep.tau_int:.1f} frames)")
        ax_acf.axhline(0, color="gray", linestyle="-", linewidth=0.8)
        ax_acf.axhline(np.exp(-1), color="orange", linestyle="--", label=r"$e^{-1}$ decay")
        ax_acf.set_xlabel("Lag (frames)")
        ax_acf.set_ylabel("Autocorrelation C(lag)")
        ax_acf.set_title(f"{obs_name} Autocorrelation (N_eff = {obs_rep.n_eff:.0f})")
        ax_acf.grid(True)
        ax_acf.legend(loc="upper right", frameon=True, fontsize=8)

    plt.tight_layout()
    for fmt in formats:
        p = os.path.join(output_dir, f"mdcheck_convergence_overview.{fmt}")
        plt.savefig(p, dpi=300, bbox_inches="tight")
        saved_files.append(p)
    plt.close()

    # 2. Block Averaging Plot
    fig, ax = plt.subplots(figsize=(7, 4.5))
    for obs_name, y in timeseries_dict.items():
        obs_rep = report.observables[obs_name]
        t_eq_idx = int(np.searchsorted(time_coords, obs_rep.t_eq_time)) if time_coords is not None else 0
        y_prod = y[t_eq_idx:]
        b_sizes, std_errs = block_averaging(y_prod)
        if len(b_sizes) > 1:
            ax.plot(b_sizes, std_errs, marker="o", markersize=4, label=f"{obs_name} (Plateau SE)")
            
    ax.set_xscale("log")
    ax.set_xlabel("Block Size (frames)")
    ax.set_ylabel("Estimated Standard Error of the Mean")
    ax.set_title("Flyvbjerg-Petersen Block Averaging Convergence")
    ax.grid(True, which="both")
    ax.legend(frameon=True)
    plt.tight_layout()
    
    for fmt in formats:
        p = os.path.join(output_dir, f"mdcheck_block_averaging.{fmt}")
        plt.savefig(p, dpi=300, bbox_inches="tight")
        saved_files.append(p)
    plt.close()

    # 3. Multi-Replica Overlap Plot (if replicas provided)
    if replica_dict is not None and len(replica_dict) > 0:
        fig, axes = plt.subplots(1, len(replica_dict), figsize=(6.0 * len(replica_dict), 4.5), squeeze=False)
        for col_idx, (obs_name, rep_list) in enumerate(replica_dict.items()):
            ax = axes[0, col_idx]
            for r_idx, r_series in enumerate(rep_list):
                ax.hist(r_series, bins=30, density=True, alpha=0.4, label=f"Replica {r_idx+1}")
            ax.set_xlabel(f"{obs_name}")
            ax.set_ylabel("Probability Density")
            ax.set_title(f"Multi-Replica Overlap: {obs_name}")
            ax.grid(True)
            ax.legend(frameon=True)
            
        plt.tight_layout()
        for fmt in formats:
            p = os.path.join(output_dir, f"mdcheck_replica_distributions.{fmt}")
            plt.savefig(p, dpi=300, bbox_inches="tight")
            saved_files.append(p)
        plt.close()

    return saved_files
