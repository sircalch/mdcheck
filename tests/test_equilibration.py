"""
Tests for equilibration detection and Geweke stationarity diagnostics.
"""

import numpy as np
import pytest
from mdcheck.core.equilibration import detect_equilibration, geweke_diagnostic


def test_equilibration_detection_synthetic():
    rng = np.random.default_rng(42)
    n = 2000
    t = np.linspace(0, 100, n)
    
    # 200 frames of initial equilibration drift, followed by stationary noise
    eq_phase = np.linspace(0.5, 2.0, 400)
    prod_phase = rng.normal(2.0, 0.1, size=1600)
    series = np.concatenate([eq_phase, prod_phase])
    
    res = detect_equilibration(series, time_series=t, step_search=10)
    
    # Should identify equilibration cutoff around frame 350-450
    assert 300 <= res["t_eq_index"] <= 500
    assert res["n_prod"] >= 1400
    assert res["fraction_discarded"] > 0.15


def test_geweke_diagnostic():
    rng = np.random.default_rng(42)
    # Stationary series
    stationary_x = rng.normal(0, 1, 1000)
    gew_stat = geweke_diagnostic(stationary_x)
    assert gew_stat["is_stationary"] is True
    assert gew_stat["p_value"] > 0.01
    
    # Non-stationary series (steep drift)
    drifting_x = np.linspace(0, 10, 1000) + rng.normal(0, 0.5, 1000)
    gew_drift = geweke_diagnostic(drifting_x)
    assert gew_drift["is_stationary"] is False
    assert abs(gew_drift["z_score"]) > 2.0
