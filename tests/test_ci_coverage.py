"""
Coverage regression test for the autocorrelation-corrected confidence interval of the mean.
"""
import numpy as np

from mdcheck.core.bootstrap import inefficiency_corrected_ci


def _ar1(phi, n, rng):
    x = np.empty(n)
    x[0] = rng.normal()
    e = rng.normal(0.0, np.sqrt(1.0 - phi * phi), n)
    for t in range(1, n):
        x[t] = phi * x[t - 1] + e[t]
    return x


def test_ci_attains_nominal_coverage_on_ar1():
    # AR(1), phi = 0.9 (g = 19), true mean 0. Nominal 95% coverage; allow sampling noise.
    rng = np.random.default_rng(11)
    hits = 0
    reps = 300
    for _ in range(reps):
        _, lo, hi = inefficiency_corrected_ci(_ar1(0.9, 5000, rng))
        hits += lo <= 0.0 <= hi
    assert 0.91 <= hits / reps <= 0.99


def test_ci_width_scales_with_inefficiency():
    rng = np.random.default_rng(3)
    x = rng.normal(size=4000)
    _, lo1, hi1 = inefficiency_corrected_ci(x, g=1.0)
    _, lo4, hi4 = inefficiency_corrected_ci(x, g=4.0)
    assert np.isclose((hi4 - lo4) / (hi1 - lo1), 2.0, rtol=0.02)
