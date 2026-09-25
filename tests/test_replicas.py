"""
Tests for multi-replica consistency, JSD, and subspace overlap (RMSIP).
"""

import numpy as np
import pytest
from mdcheck.core.replicas import (
    jensen_shannon_distance,
    subspace_overlap_rmsip,
    assess_replica_consistency
)


def test_jensen_shannon_distance():
    rng = np.random.default_rng(42)
    r1 = rng.normal(1.5, 0.2, size=1000)
    r2 = rng.normal(1.5, 0.2, size=1000)
    r3 = rng.normal(3.0, 0.2, size=1000)  # Diverged state
    
    # Identical distributions should have small JSD
    jsd_12 = jensen_shannon_distance(r1, r2)
    assert jsd_12 < 0.15
    
    # Non-overlapping distributions should have large JSD
    jsd_13 = jensen_shannon_distance(r1, r3)
    assert jsd_13 > 0.50


def test_subspace_overlap_rmsip():
    # Orthogonal identical vs rotated eigenvectors
    n_features = 20
    n_modes = 5
    
    # Perfect identity
    v_a = np.eye(n_features)[:, :n_modes]
    v_b = np.eye(n_features)[:, :n_modes]
    assert np.isclose(subspace_overlap_rmsip(v_a, v_b, n_modes=5), 1.0)
    
    # Orthogonal subspaces
    v_c = np.eye(n_features)[:, 5:10]
    assert np.isclose(subspace_overlap_rmsip(v_a, v_c, n_modes=5), 0.0)


def test_assess_replica_consistency():
    rng = np.random.default_rng(42)
    r1 = rng.normal(1.5, 0.1, 1000)
    r2 = rng.normal(1.5, 0.1, 1000)
    r3 = rng.normal(1.5, 0.1, 1000)
    
    res = assess_replica_consistency([r1, r2, r3], replica_names=["R1", "R2", "R3"])
    # Same ensemble: a calibrated test may WARN ~5% of the time (p = 0.044 for this seed) but must not FAIL.
    assert res["status"] != "FAIL"
    assert res["heterogeneity_p_value"] > 0.01
    assert res["n_replicas"] == 3
    assert res["mean_jsd"] < 0.15


def _ar1(phi, n, rng, mean=0.0):
    x = np.empty(n)
    x[0] = rng.normal()
    e = rng.normal(0.0, np.sqrt(1.0 - phi * phi), n)
    for t in range(1, n):
        x[t] = phi * x[t - 1] + e[t]
    return x + mean


def test_replica_test_false_alarm_rate_is_calibrated_for_correlated_data():
    # Three correlated replicas (g = 19) of the same ensemble: status != PASS should occur ~5%.
    rng = np.random.default_rng(5)
    reps = 400
    alarms = 0
    for _ in range(reps):
        res = assess_replica_consistency([_ar1(0.9, 3000, rng) for _ in range(3)])
        alarms += res["status"] != "PASS"
    assert 0.02 <= alarms / reps <= 0.10


def test_replica_test_detects_shifted_replica():
    # One replica displaced by ~5 standard errors of its mean.
    rng = np.random.default_rng(9)
    se = np.sqrt(19.0 / 3000)
    reps = [_ar1(0.9, 3000, rng), _ar1(0.9, 3000, rng), _ar1(0.9, 3000, rng, mean=5 * se)]
    res = assess_replica_consistency(reps)
    assert res["status"] == "FAIL"
    assert res["heterogeneity_p_value"] < 0.01
