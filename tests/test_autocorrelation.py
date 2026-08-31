"""
Tests for autocorrelation and effective sample size calculations.
"""

import numpy as np
import pytest
from mdcheck.core.autocorrelation import (
    compute_autocorrelation,
    integrated_autocorrelation_time,
    effective_sample_size,
    statistical_inefficiency
)


def test_autocorrelation_white_noise():
    rng = np.random.default_rng(42)
    x = rng.normal(0, 1, size=5000)
    
    acf = compute_autocorrelation(x, max_lag=50)
    assert np.isclose(acf[0], 1.0)
    # Lags > 0 should be near 0 for uncorrelated noise
    assert np.all(np.abs(acf[1:]) < 0.1)
    
    tau, g, _ = integrated_autocorrelation_time(x)
    assert 0.4 <= tau <= 0.8
    assert 0.8 <= g <= 1.6
    
    n_eff, _, _ = effective_sample_size(x)
    assert n_eff > 3000


def test_autocorrelation_ar1_process():
    # AR(1) process with known phi = 0.8 -> theoretical tau_int = (1 + phi) / (2 * (1 - phi)) = 1.8 / 0.4 = 4.5
    rng = np.random.default_rng(123)
    n = 10000
    phi = 0.8
    x = np.zeros(n)
    for i in range(1, n):
        x[i] = phi * x[i-1] + np.sqrt(1 - phi**2) * rng.normal(0, 1)
        
    tau, g, _ = integrated_autocorrelation_time(x)
    # Check that estimated tau_int is reasonably close to theoretical 4.5
    assert 3.0 <= tau <= 6.5
    assert g > 5.0
    
    n_eff, _, _ = effective_sample_size(x)
    assert n_eff < n / 3.0


def test_edge_cases():
    # Constant series
    const_arr = np.ones(100)
    acf = compute_autocorrelation(const_arr)
    assert len(acf) == 100
    
    # Short array
    short_arr = np.array([1.0, 2.0])
    tau, g, _ = integrated_autocorrelation_time(short_arr)
    assert tau == 0.5
    assert g == 1.0
