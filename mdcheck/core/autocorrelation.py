"""
Integrated autocorrelation time, statistical inefficiency, and effective sample size.
"""

from typing import Tuple, Optional
import numpy as np


def compute_autocorrelation(series: np.ndarray, max_lag: Optional[int] = None) -> np.ndarray:
    """
    Computes the normalized autocorrelation function (ACF) of a 1D timeseries using FFT.

    Parameters
    ----------
    series : np.ndarray
        1D array of timeseries values.
    max_lag : int, optional
        Maximum lag to compute. If None, computes for all possible lags.

    Returns
    -------
    acf : np.ndarray
        Normalized autocorrelation function with acf[0] = 1.0.
    """
    x = np.asarray(series, dtype=np.float64)
    if x.ndim != 1:
        x = x.ravel()
    
    n = len(x)
    if n <= 1:
        return np.array([1.0], dtype=np.float64)
    
    x_mean = np.mean(x)
    x_zero_mean = x - x_mean
    variance = np.var(x, ddof=0)
    
    if variance == 0 or np.isnan(variance):
        return np.ones(max_lag if max_lag else n, dtype=np.float64)
    
    # Zero-padding for linear autocorrelation via FFT
    n_fft = 2 ** int(np.ceil(np.log2(2 * n - 1)))
    fx = np.fft.fft(x_zero_mean, n=n_fft)
    px = fx * np.conjugate(fx)
    autocov = np.fft.ifft(px).real[:n]
    
    # Scale by sample count at each lag for unbiased/consistent estimator
    lags = np.arange(n)
    denom = (n - lags) * variance
    acf = autocov / denom
    
    # Normalization guarantee
    if acf[0] != 0:
        acf = acf / acf[0]
    else:
        acf[0] = 1.0
        
    if max_lag is not None and max_lag < n:
        acf = acf[:max_lag + 1]
        
    return acf


def integrated_autocorrelation_time(
    series: np.ndarray,
    c_window: float = 6.0,
    max_lag: Optional[int] = None
) -> Tuple[float, float, int]:
    """
    Calculates the integrated autocorrelation time (tau_int) using a self-consistent
    Madras-Sokal windowing procedure.

    Formula:
        tau_int = 0.5 + sum_{k=1}^M C(k)
    where M is the smallest integer such that M >= c_window * tau_int.

    Parameters
    ----------
    series : np.ndarray
        1D timeseries array.
    c_window : float, default 6.0
        Madras-Sokal window parameter (typically between 4.0 and 8.0).
    max_lag : int, optional
        Maximum lag to explore.

    Returns
    -------
    tau_int : float
        Integrated autocorrelation time in units of steps/frames.
    g : float
        Statistical inefficiency factor g = 1 + 2 * tau_int.
    window_m : int
        Window cutoff index M.
    """
    x = np.asarray(series, dtype=np.float64)
    n = len(x)
    
    if n < 4:
        return 0.5, 1.0, 0
    
    if max_lag is None:
        max_lag = min(n // 2, 50000)
        
    acf = compute_autocorrelation(x, max_lag=max_lag)
    
    # Initial estimate
    tau_int = 0.5
    window_m = 0
    
    for k in range(1, len(acf)):
        # Stop if autocorrelation becomes negative or noise-dominated
        if acf[k] < 0:
            window_m = k
            break
            
        tau_int += acf[k]
        
        # Self-consistent window check: M >= c * tau_int
        if k >= c_window * tau_int:
            window_m = k
            break
    else:
        window_m = len(acf) - 1

    # Bounds guarantee
    tau_int = max(0.5, float(tau_int))
    g = max(1.0, float(1.0 + 2.0 * (tau_int - 0.5)))
    
    return tau_int, g, int(window_m)


def statistical_inefficiency(series: np.ndarray, c_window: float = 6.0) -> float:
    """
    Calculates the statistical inefficiency factor g.

    Parameters
    ----------
    series : np.ndarray
        1D timeseries array.
    c_window : float
        Window parameter.

    Returns
    -------
    g : float
        Statistical inefficiency g = 1 + 2*tau_int.
    """
    _, g, _ = integrated_autocorrelation_time(series, c_window=c_window)
    return g


def effective_sample_size(series: np.ndarray, c_window: float = 6.0) -> Tuple[float, float, float]:
    """
    Computes effective sample size N_eff = N / g.

    Parameters
    ----------
    series : np.ndarray
        1D timeseries array.
    c_window : float
        Madras-Sokal window parameter.

    Returns
    -------
    n_eff : float
        Effective number of uncorrelated samples.
    tau_int : float
        Integrated autocorrelation time.
    g : float
        Statistical inefficiency factor.
    """
    x = np.asarray(series, dtype=np.float64)
    n = len(x)
    if n == 0:
        return 0.0, 0.5, 1.0
        
    tau_int, g, _ = integrated_autocorrelation_time(x, c_window=c_window)
    n_eff = float(n) / max(1.0, g)
    return n_eff, tau_int, g
