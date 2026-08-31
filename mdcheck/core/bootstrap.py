"""
Block bootstrap and non-parametric confidence interval estimation.
"""

from typing import Tuple, Callable, Optional
import numpy as np
from mdcheck.core.autocorrelation import integrated_autocorrelation_time


def bootstrap_ci(
    series: np.ndarray,
    stat_func: Callable[[np.ndarray], float] = np.mean,
    n_resamples: int = 2000,
    confidence_level: float = 0.95,
    block_size: Optional[int] = None,
    random_state: Optional[int] = 42
) -> Tuple[float, float, float]:
    """
    Computes moving-block bootstrap confidence intervals accounting for autocorrelation.

    Parameters
    ----------
    series : np.ndarray
        1D timeseries array.
    stat_func : callable, default np.mean
        Function computing the target statistic on a 1D array.
    n_resamples : int, default 2000
        Number of bootstrap resamples.
    confidence_level : float, default 0.95
        Confidence level (e.g. 0.95 for 95% CI).
    block_size : int, optional
        Block length. If None, automatically determined as int(2 * tau_int + 1).
    random_state : int, optional
        Random seed for reproducibility.

    Returns
    -------
    point_estimate : float
        Observed sample statistic.
    ci_lower : float
        Lower bound of bootstrap confidence interval.
    ci_upper : float
        Upper bound of bootstrap confidence interval.
    """
    x = np.asarray(series, dtype=np.float64)
    n = len(x)
    
    if n == 0:
        return 0.0, 0.0, 0.0
        
    point_estimate = float(stat_func(x))
    
    if n < 5:
        return point_estimate, point_estimate, point_estimate
        
    if block_size is None:
        tau_int, _, _ = integrated_autocorrelation_time(x)
        block_size = max(1, int(np.ceil(2.0 * tau_int)))
        
    rng = np.random.default_rng(random_state)
    
    if block_size <= 1:
        # Standard i.i.d. bootstrap
        resample_indices = rng.integers(0, n, size=(n_resamples, n))
        resamples = x[resample_indices]
        boot_stats = np.apply_along_axis(stat_func, 1, resamples)
    else:
        # Moving block bootstrap
        n_blocks = int(np.ceil(n / block_size))
        max_start = n - block_size + 1
        
        boot_stats = np.empty(n_resamples, dtype=np.float64)
        for b in range(n_resamples):
            start_indices = rng.integers(0, max_start, size=n_blocks)
            sampled_blocks = np.concatenate([x[idx : idx + block_size] for idx in start_indices])
            sampled_series = sampled_blocks[:n]
            boot_stats[b] = stat_func(sampled_series)
            
    alpha = (1.0 - confidence_level) / 2.0
    ci_lower = float(np.percentile(boot_stats, 100.0 * alpha))
    ci_upper = float(np.percentile(boot_stats, 100.0 * (1.0 - alpha)))
    
    return point_estimate, ci_lower, ci_upper
