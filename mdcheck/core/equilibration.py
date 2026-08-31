"""
Automated equilibration detection and stationarity diagnostics.
"""

from typing import Dict, Any, Optional
import numpy as np
from scipy import stats
from mdcheck.core.autocorrelation import integrated_autocorrelation_time


def detect_equilibration(
    series: np.ndarray,
    time_series: Optional[np.ndarray] = None,
    step_search: int = 10,
    c_window: float = 6.0
) -> Dict[str, Any]:
    """
    Automatically detects the equilibration time (t_eq) of a timeseries by maximizing
    the effective sample size of the production region:
        t_eq = argmax_{t_0} N_eff(t_0) = argmax_{t_0} (N - t_0) / g(t_0)

    Parameters
    ----------
    series : np.ndarray
        1D array of the observable (e.g., RMSD, Energy, Rg).
    time_series : np.ndarray, optional
        Array of corresponding timestamps (e.g. in ns or ps). If None, frame indices are used.
    step_search : int, default 10
        Step size between candidate t_0 cutoff points to speed up search.
    c_window : float, default 6.0
        Window parameter for autocorrelation time calculation.

    Returns
    -------
    result : dict
        Dictionary containing:
        - 't_eq_index': Index of first production frame
        - 't_eq_time': Timestamp of equilibration
        - 'n_total': Total number of frames
        - 'n_prod': Number of production frames
        - 'n_eff': Effective sample size in production
        - 'tau_int': Integrated autocorrelation time in production (frames)
        - 'g': Statistical inefficiency factor in production
        - 'fraction_discarded': Proportion of trajectory discarded as equilibration
        - 'neff_curve': Effective sample size as function of cutoff index
    """
    y = np.asarray(series, dtype=np.float64)
    n = len(y)
    
    if time_series is None:
        t = np.arange(n, dtype=np.float64)
    else:
        t = np.asarray(time_series, dtype=np.float64)
        
    if n < 20:
        tau, g, _ = integrated_autocorrelation_time(y, c_window=c_window)
        return {
            "t_eq_index": 0,
            "t_eq_time": float(t[0]),
            "n_total": n,
            "n_prod": n,
            "n_eff": float(n / max(1.0, g)),
            "tau_int": float(tau),
            "g": float(g),
            "fraction_discarded": 0.0,
            "neff_curve": [(0, float(n / max(1.0, g)))]
        }
        
    # Search grid for t_0 (test up to 85% of trajectory)
    max_test_idx = int(0.85 * n)
    indices = np.arange(0, max_test_idx, max(1, step_search))
    
    best_t0_idx = 0
    best_n_eff = -1.0
    best_tau = 0.5
    best_g = 1.0
    
    neff_curve = []
    
    for idx in indices:
        y_sub = y[idx:]
        n_sub = len(y_sub)
        if n_sub < 10:
            break
            
        tau, g, _ = integrated_autocorrelation_time(y_sub, c_window=c_window)
        n_eff = float(n_sub) / max(1.0, g)
        
        neff_curve.append((int(idx), float(n_eff)))
        
        if n_eff > best_n_eff:
            best_n_eff = n_eff
            best_t0_idx = idx
            best_tau = tau
            best_g = g

    # Refined local search around best_t0_idx if step_search > 1
    if step_search > 1 and best_t0_idx > 0:
        local_start = max(0, best_t0_idx - step_search)
        local_end = min(max_test_idx, best_t0_idx + step_search)
        for idx in range(local_start, local_end):
            y_sub = y[idx:]
            n_sub = len(y_sub)
            if n_sub < 10:
                break
            tau, g, _ = integrated_autocorrelation_time(y_sub, c_window=c_window)
            n_eff = float(n_sub) / max(1.0, g)
            if n_eff > best_n_eff:
                best_n_eff = n_eff
                best_t0_idx = idx
                best_tau = tau
                best_g = g

    fraction_discarded = float(best_t0_idx) / float(n)
    n_prod = n - best_t0_idx
    
    return {
        "t_eq_index": int(best_t0_idx),
        "t_eq_time": float(t[best_t0_idx]),
        "n_total": int(n),
        "n_prod": int(n_prod),
        "n_eff": float(best_n_eff),
        "tau_int": float(best_tau),
        "g": float(best_g),
        "fraction_discarded": float(fraction_discarded),
        "neff_curve": neff_curve
    }


def geweke_diagnostic(
    series: np.ndarray,
    first_fraction: float = 0.1,
    last_fraction: float = 0.5
) -> Dict[str, Any]:
    """
    Computes Geweke's stationarity diagnostic by comparing the sample mean
    of the first part of the chain with the last part of the chain.

    Parameters
    ----------
    series : np.ndarray
        1D timeseries array (typically in the production phase).
    first_fraction : float, default 0.1
        Fraction of data in the early segment.
    last_fraction : float, default 0.5
        Fraction of data in the late segment.

    Returns
    -------
    result : dict
        Contains 'z_score', 'p_value', and 'is_stationary' (p > 0.05).
    """
    x = np.asarray(series, dtype=np.float64)
    n = len(x)
    
    if n < 20:
        return {"z_score": 0.0, "p_value": 1.0, "is_stationary": True}
        
    n_first = max(5, int(n * first_fraction))
    n_last = max(10, int(n * last_fraction))
    
    x_first = x[:n_first]
    x_last = x[-n_last:]
    
    mean_first = np.mean(x_first)
    mean_last = np.mean(x_last)
    
    # Estimate spectral variance / effective variance using autocorrelation
    _, g_first, _ = integrated_autocorrelation_time(x_first)
    _, g_last, _ = integrated_autocorrelation_time(x_last)
    
    var_first = (np.var(x_first, ddof=1) * g_first) / n_first if n_first > 1 else 0.0
    var_last = (np.var(x_last, ddof=1) * g_last) / n_last if n_last > 1 else 0.0
    
    se_diff = np.sqrt(var_first + var_last)
    
    if se_diff == 0:
        z_score = 0.0
    else:
        z_score = (mean_first - mean_last) / se_diff
        
    p_value = 2.0 * (1.0 - stats.norm.cdf(abs(z_score)))
    
    # Stationarity holds if |z| < 2 (p > 0.05)
    is_stationary = bool(abs(z_score) < 1.96)
    
    return {
        "z_score": float(z_score),
        "p_value": float(p_value),
        "is_stationary": is_stationary,
        "mean_first": float(mean_first),
        "mean_last": float(mean_last)
    }
