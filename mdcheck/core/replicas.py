"""
Multi-replica consistency, distribution divergence, and conformational overlap.
"""

from typing import List, Dict, Any, Optional
import numpy as np
from scipy import stats
from scipy.spatial.distance import jensenshannon


def jensen_shannon_distance(
    series_a: np.ndarray,
    series_b: np.ndarray,
    n_bins: int = 50
) -> float:
    """
    Calculates the Jensen-Shannon distance (square root of JSD, bounded in [0, 1])
    between the empirical probability distributions of two replica timeseries.

    Parameters
    ----------
    series_a : np.ndarray
        Data from replica A.
    series_b : np.ndarray
        Data from replica B.
    n_bins : int, default 50
        Number of histogram bins for discrete density estimation.

    Returns
    -------
    js_dist : float
        Jensen-Shannon distance in [0, 1]. Values < 0.15 indicate excellent overlap.
    """
    a = np.asarray(series_a, dtype=np.float64)
    b = np.asarray(series_b, dtype=np.float64)
    
    val_min = min(np.min(a), np.min(b))
    val_max = max(np.max(a), np.max(b))
    
    if val_min == val_max:
        return 0.0
        
    bins = np.linspace(val_min, val_max, n_bins + 1)
    
    hist_a, _ = np.histogram(a, bins=bins, density=True)
    hist_b, _ = np.histogram(b, bins=bins, density=True)
    
    # Smooth with small epsilon to avoid numerical log(0) issues
    eps = 1e-12
    p = hist_a + eps
    q = hist_b + eps
    
    p = p / np.sum(p)
    q = q / np.sum(q)
    
    js_dist = jensenshannon(p, q, base=2.0)
    return float(np.clip(js_dist, 0.0, 1.0))


def subspace_overlap_rmsip(
    eigenvectors_a: np.ndarray,
    eigenvectors_b: np.ndarray,
    n_modes: int = 10
) -> float:
    """
    Computes the Root-Mean-Square Inner Product (RMSIP) between the principal
    component subspaces of two trajectory replicas:

        RMSIP = sqrt( (1 / s) * sum_{i=1}^s sum_{j=1}^s (v_i^A . v_j^B)^2 )

    Parameters
    ----------
    eigenvectors_a : np.ndarray
        Matrix of eigenvectors from PCA of replica A (shape: [n_features, n_modes]).
    eigenvectors_b : np.ndarray
        Matrix of eigenvectors from PCA of replica B (shape: [n_features, n_modes]).
    n_modes : int, default 10
        Number of essential modes s to compare.

    Returns
    -------
    rmsip : float
        Subspace overlap in [0, 1]. RMSIP > 0.70 represents strong essential subspace convergence.
    """
    va = np.asarray(eigenvectors_a, dtype=np.float64)
    vb = np.asarray(eigenvectors_b, dtype=np.float64)
    
    s = min(n_modes, va.shape[1], vb.shape[1])
    if s == 0:
        return 1.0
        
    va_sub = va[:, :s]
    vb_sub = vb[:, :s]
    
    # Dot products matrix between mode i and mode j
    dot_products = np.dot(va_sub.T, vb_sub)
    squared_sum = np.sum(dot_products ** 2)
    
    rmsip = np.sqrt(squared_sum / s)
    return float(np.clip(rmsip, 0.0, 1.0))


def assess_replica_consistency(
    replica_data: List[np.ndarray],
    replica_names: Optional[List[str]] = None,
    jsd_threshold_pass: float = 0.15,
    jsd_threshold_warn: float = 0.30
) -> Dict[str, Any]:
    """
    Evaluates consistency across multiple replica trajectories.

    Parameters
    ----------
    replica_data : list of np.ndarray
        List of 1D timeseries (one per replica).
    replica_names : list of str, optional
        Names or identifiers for each replica (e.g. ['R1', 'R2', 'R3']).
    jsd_threshold_pass : float, default 0.15
        JSD distance cutoff for PASS status.
    jsd_threshold_warn : float, default 0.30
        JSD distance cutoff for WARNING status.

    Returns
    -------
    result : dict
        Multi-replica consistency metrics, pairwise matrix, and status.
    """
    n_reps = len(replica_data)
    if replica_names is None:
        replica_names = [f"Replica_{i+1}" for i in range(n_reps)]
        
    if n_reps < 2:
        return {
            "n_replicas": n_reps,
            "status": "PASS",
            "mean_jsd": 0.0,
            "max_jsd": 0.0,
            "pairwise_jsd": {},
            "pairwise_ks": {},
            "wasserstein_distances": {},
            "replica_means": [float(np.mean(r)) for r in replica_data],
            "replica_stds": [float(np.std(r)) for r in replica_data],
            "recommendation": "Single replica provided. Multi-replica testing recommended for publishable rigor."
        }

    pairwise_jsd = {}
    pairwise_ks = {}
    pairwise_wasserstein = {}
    jsd_values = []
    
    for i in range(n_reps):
        for j in range(i + 1, n_reps):
            pair_key = f"{replica_names[i]} vs {replica_names[j]}"
            r_i = replica_data[i]
            r_j = replica_data[j]
            
            # JS Distance
            js_dist = jensen_shannon_distance(r_i, r_j)
            pairwise_jsd[pair_key] = js_dist
            jsd_values.append(js_dist)
            
            # KS Test
            ks_res = stats.ks_2samp(r_i, r_j)
            pairwise_ks[pair_key] = {
                "statistic": float(ks_res.statistic),
                "p_value": float(ks_res.pvalue)
            }
            
            # Wasserstein Distance
            w_dist = stats.wasserstein_distance(r_i, r_j)
            pairwise_wasserstein[pair_key] = float(w_dist)

    mean_jsd = float(np.mean(jsd_values))
    max_jsd = float(np.max(jsd_values))
    
    if max_jsd <= jsd_threshold_pass:
        status = "PASS"
        recommendation = "Replicas exhibit high distributional overlap and conformational reproducibility."
    elif max_jsd <= jsd_threshold_warn:
        status = "WARNING"
        recommendation = "Moderate divergence detected between replicas. Check individual pairwise distributions."
    else:
        status = "FAIL"
        recommendation = "Significant divergence between replicas. Simulation has not sampled identical conformational states."

    return {
        "n_replicas": n_reps,
        "status": status,
        "mean_jsd": mean_jsd,
        "max_jsd": max_jsd,
        "pairwise_jsd": pairwise_jsd,
        "pairwise_ks": pairwise_ks,
        "wasserstein_distances": pairwise_wasserstein,
        "replica_means": [float(np.mean(r)) for r in replica_data],
        "replica_stds": [float(np.std(r)) for r in replica_data],
        "recommendation": recommendation
    }
