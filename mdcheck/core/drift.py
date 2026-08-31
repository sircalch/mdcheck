"""
Energy and structural drift assessment and block averaging uncertainty estimation.
"""

from typing import Dict, Any, Optional, Tuple, List
import numpy as np
from scipy import stats


def block_averaging(
    series: np.ndarray,
    min_blocks: int = 4,
    max_blocks: Optional[int] = None
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Performs Flyvbjerg-Petersen block averaging to estimate the true standard error
    as a function of block size.

    Parameters
    ----------
    series : np.ndarray
        1D timeseries array.
    min_blocks : int, default 4
        Minimum number of blocks to keep for valid variance estimation.
    max_blocks : int, optional
        Maximum number of blocks (default: n // 2).

    Returns
    -------
    block_sizes : np.ndarray
        Array of block sizes (number of frames per block).
    std_errors : np.ndarray
        Estimated standard error of the mean for each block size.
    """
    x = np.asarray(series, dtype=np.float64)
    n = len(x)
    
    if n < min_blocks * 2:
        return np.array([1]), np.array([np.std(x, ddof=1) / np.sqrt(max(1, n))])
        
    if max_blocks is None:
        max_blocks = n // 2
        
    block_sizes_list = []
    std_errors_list = []
    
    # Test logarithmically spaced block sizes
    max_block_size = n // min_blocks
    sizes = np.unique(np.logspace(0, np.log10(max_block_size), num=30, dtype=int))
    
    for b_size in sizes:
        if b_size < 1:
            continue
        n_b = n // b_size
        if n_b < min_blocks:
            break
            
        # Truncate to exact multiple
        usable_data = x[:n_b * b_size]
        blocks = usable_data.reshape(n_b, b_size)
        block_means = np.mean(blocks, axis=1)
        
        # Standard error of the mean of block means
        se = np.std(block_means, ddof=1) / np.sqrt(n_b)
        
        block_sizes_list.append(b_size)
        std_errors_list.append(se)
        
    return np.array(block_sizes_list), np.array(std_errors_list)


def assess_drift(
    series: np.ndarray,
    time_series: Optional[np.ndarray] = None,
    relative_drift_threshold_pass: float = 0.05,
    relative_drift_threshold_warn: float = 0.15
) -> Dict[str, Any]:
    """
    Estimates systematic linear drift and tests for statistical stationarity in production.

    Parameters
    ----------
    series : np.ndarray
        1D timeseries array in production.
    time_series : np.ndarray, optional
        Corresponding time coordinates (e.g. ns).
    relative_drift_threshold_pass : float, default 0.05
        Maximum acceptable |drift| / mean ratio (5% change over trajectory length).
    relative_drift_threshold_warn : float, default 0.15
        Warning threshold (15% change).

    Returns
    -------
    result : dict
        Drift diagnostics, slope, p-value, relative change, and PASS/WARNING/FAIL status.
    """
    y = np.asarray(series, dtype=np.float64)
    n = len(y)
    
    if time_series is None:
        t = np.arange(n, dtype=np.float64)
    else:
        t = np.asarray(time_series, dtype=np.float64)
        
    if n < 4:
        return {
            "status": "PASS",
            "slope": 0.0,
            "intercept": float(np.mean(y)) if n > 0 else 0.0,
            "r_squared": 0.0,
            "p_value": 1.0,
            "total_drift": 0.0,
            "relative_drift": 0.0,
            "duration": 0.0,
            "recommendation": "Too few frames to evaluate drift."
        }

    # Linear regression
    slope, intercept, r_value, p_value, std_err = stats.linregress(t, y)
    r_squared = float(r_value ** 2)
    
    duration = float(t[-1] - t[0])
    total_drift = float(slope * duration)
    
    mean_val = float(np.mean(y))
    scale = max(abs(mean_val), float(np.std(y, ddof=1)), 1e-6)
    relative_drift = float(abs(total_drift) / scale)
    
    # Determine pass/warn/fail status
    # If slope is not statistically significant (p > 0.05) or relative drift is tiny -> PASS
    if p_value > 0.05 or relative_drift <= relative_drift_threshold_pass:
        status = "PASS"
        recommendation = "No critical systematic drift detected; timeseries is stationary."
    elif relative_drift <= relative_drift_threshold_warn:
        status = "WARNING"
        recommendation = "Mild systematic drift observed across the trajectory."
    else:
        status = "FAIL"
        recommendation = "Substantial systematic drift detected; the system has not reached stationary equilibrium."

    return {
        "status": status,
        "slope": float(slope),
        "intercept": float(intercept),
        "r_squared": float(r_squared),
        "p_value": float(p_value),
        "std_err": float(std_err),
        "total_drift": float(total_drift),
        "relative_drift": float(relative_drift),
        "duration": duration,
        "recommendation": recommendation
    }
