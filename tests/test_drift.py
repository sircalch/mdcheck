"""
Tests for drift evaluation and block averaging.
"""

import numpy as np
import pytest
from mdcheck.core.drift import assess_drift, block_averaging


def test_drift_stationary():
    rng = np.random.default_rng(42)
    t = np.linspace(0, 100, 1000)
    y = rng.normal(5.0, 0.2, 1000)
    
    res = assess_drift(y, time_series=t)
    assert res["status"] == "PASS"
    assert res["relative_drift"] < 0.05


def test_drift_significant():
    rng = np.random.default_rng(42)
    t = np.linspace(0, 100, 1000)
    # 50% change from start to end
    y = 5.0 + 0.03 * t + rng.normal(0, 0.1, 1000)
    
    res = assess_drift(y, time_series=t)
    assert res["status"] in ["WARNING", "FAIL"]
    assert res["relative_drift"] > 0.10


def test_block_averaging():
    rng = np.random.default_rng(42)
    y = rng.normal(0, 1, 1000)
    b_sizes, std_errs = block_averaging(y)
    assert len(b_sizes) > 1
    assert len(std_errs) == len(b_sizes)
    assert np.all(std_errs > 0)
