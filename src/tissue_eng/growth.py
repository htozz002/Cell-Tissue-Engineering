"""Cell population growth kinetics.

Time and rate units only need to be consistent (e.g. hours and 1/h).
Functions that take ``t`` accept scalars or array-likes and broadcast.
"""

import math

import numpy as np


def _require_positive(**values):
    for name, value in values.items():
        if np.any(np.asarray(value) <= 0):
            raise ValueError(f"{name} must be positive")


def exponential(t, n0, mu):
    """Unrestricted growth, N(t) = N0 exp(mu t)."""
    _require_positive(n0=n0)
    return n0 * np.exp(mu * np.asarray(t, dtype=float))


def logistic(t, n0, r, k):
    """Logistic growth toward carrying capacity ``k`` (e.g. confluence).

    N(t) = K / (1 + ((K - N0) / N0) exp(-r t))
    """
    _require_positive(n0=n0, k=k)
    t = np.asarray(t, dtype=float)
    return k / (1.0 + (k - n0) / n0 * np.exp(-r * t))


def gompertz(t, n0, r, k):
    """Gompertz growth: N(t) = K exp(ln(N0/K) exp(-r t)).

    Growth rate decays exponentially; often fits 3D aggregates and tumours
    better than the logistic model.
    """
    _require_positive(n0=n0, k=k)
    t = np.asarray(t, dtype=float)
    return k * np.exp(np.log(n0 / k) * np.exp(-r * t))


def specific_growth_rate(n1, n2, dt):
    """Specific growth rate mu = ln(N2/N1) / dt from two counts in log phase."""
    _require_positive(n1=n1, n2=n2, dt=dt)
    return math.log(n2 / n1) / dt


def doubling_time(mu):
    """Population doubling time t_d = ln 2 / mu."""
    _require_positive(mu=mu)
    return math.log(2.0) / mu


def population_doublings(n_initial, n_final):
    """Number of population doublings, PD = log2(N_final / N_initial).

    Sum PD over successive passages for cumulative population doublings (CPD).
    """
    _require_positive(n_initial=n_initial, n_final=n_final)
    return math.log2(n_final / n_initial)


def monod(s, mu_max, ks):
    """Monod substrate-limited growth rate, mu = mu_max S / (Ks + S)."""
    _require_positive(mu_max=mu_max, ks=ks)
    s = np.asarray(s, dtype=float)
    if np.any(s < 0):
        raise ValueError("substrate concentration must be non-negative")
    return mu_max * s / (ks + s)
