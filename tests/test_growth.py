import math

import numpy as np
import pytest

from tissue_eng import growth


def test_exponential_doubles_after_doubling_time():
    mu = 0.03
    td = growth.doubling_time(mu)
    assert growth.exponential(td, 1e5, mu) == pytest.approx(2e5)


def test_rate_doubling_time_and_pd_are_consistent():
    # 5e5 cells grow to 3.2e6 in 72 h (worked example in notes/02).
    mu = growth.specific_growth_rate(5e5, 3.2e6, 72.0)
    assert mu == pytest.approx(math.log(6.4) / 72.0)
    assert growth.doubling_time(mu) == pytest.approx(26.9, abs=0.05)
    assert growth.population_doublings(5e5, 3.2e6) == pytest.approx(2.678, abs=1e-3)
    assert 72.0 / growth.population_doublings(5e5, 3.2e6) == pytest.approx(growth.doubling_time(mu))


def test_logistic_limits():
    t = np.array([0.0, 1e4])
    n = growth.logistic(t, n0=1e4, r=0.05, k=1e6)
    assert n[0] == pytest.approx(1e4)
    assert n[1] == pytest.approx(1e6)


def test_logistic_matches_exponential_early():
    t = 1.0
    assert growth.logistic(t, 1.0, 0.1, 1e9) == pytest.approx(growth.exponential(t, 1.0, 0.1), rel=1e-6)


def test_gompertz_limits():
    n = growth.gompertz(np.array([0.0, 1e4]), n0=10.0, r=0.1, k=1000.0)
    assert n[0] == pytest.approx(10.0)
    assert n[1] == pytest.approx(1000.0)


def test_monod_half_saturation():
    assert growth.monod(2.0, mu_max=0.04, ks=2.0) == pytest.approx(0.02)


@pytest.mark.parametrize("bad", [0.0, -1.0])
def test_rejects_non_positive(bad):
    with pytest.raises(ValueError):
        growth.doubling_time(bad)
    with pytest.raises(ValueError):
        growth.population_doublings(bad, 10.0)
