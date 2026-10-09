"""Steady-state oxygen (or nutrient) transport in cell-laden constructs.

Units are SI: concentration in mol/m^3, diffusivity in m^2/s, volumetric
consumption rate in mol/(m^3 s), lengths in m.

Geometry convention: ``r`` runs from the construct's centre (r = 0, a plane
or axis of symmetry with no flux) to its exposed surface (r = R), where the
concentration is held at the bulk value ``c0``. ``geometry`` is one of
``"slab"``, ``"cylinder"`` or ``"sphere"``. For a slab fed from one face
only, R is the full thickness; for a slab fed from both faces, R is the
half-thickness.

Governing equation (m = 0, 1, 2 for slab, cylinder, sphere)::

    D (1/r^m) d/dr (r^m dC/dr) = R(C)

with R(C) = Q for zero-order consumption or Vmax C / (Km + C) for
Michaelis-Menten consumption.
"""

import math

import numpy as np

from .constants import O2_SOLUBILITY_MEDIUM

_EXPONENT = {"slab": 0, "cylinder": 1, "sphere": 2}


def _exponent(geometry):
    try:
        return _EXPONENT[geometry]
    except KeyError:
        raise ValueError(
            f"geometry must be one of {sorted(_EXPONENT)}, got {geometry!r}"
        ) from None


def _require_positive(**values):
    for name, value in values.items():
        if value <= 0:
            raise ValueError(f"{name} must be positive")


def concentration_from_partial_pressure(p_o2_mmhg, solubility=O2_SOLUBILITY_MEDIUM):
    """Dissolved O2 concentration from partial pressure via Henry's law, C = alpha pO2."""
    return solubility * np.asarray(p_o2_mmhg, dtype=float)


def volumetric_consumption_rate(cell_density, per_cell_rate):
    """Volumetric consumption Q (mol/(m^3 s)) from cells/m^3 and mol/(cell s)."""
    return cell_density * per_cell_rate


def thiele_modulus(radius, c0, diffusivity, consumption_rate):
    """Zero-order Thiele modulus, phi = R sqrt(Q / (D C0)).

    Compares consumption to diffusive supply. The centre stays oxygenated
    while phi^2 <= 2, 4 or 6 for a slab, cylinder or sphere.
    """
    _require_positive(radius=radius, c0=c0, diffusivity=diffusivity,
                      consumption_rate=consumption_rate)
    return radius * math.sqrt(consumption_rate / (diffusivity * c0))


def critical_dimension(c0, diffusivity, consumption_rate, geometry="slab"):
    """Largest R for which zero-order consumption leaves C >= 0 everywhere.

    R_crit = sqrt(2 (m + 1) D C0 / Q): sqrt(2DC0/Q) for a slab,
    sqrt(4DC0/Q) for a cylinder and sqrt(6DC0/Q) for a sphere.
    """
    m = _exponent(geometry)
    _require_positive(c0=c0, diffusivity=diffusivity, consumption_rate=consumption_rate)
    return math.sqrt(2 * (m + 1) * diffusivity * c0 / consumption_rate)


def _depletion_shape(r, rc, m):
    """(2(m+1)D/Q) * C(r) for a zero-order profile that reaches zero at r = rc > 0."""
    if m == 0:
        return (r - rc) ** 2
    if m == 1:
        return r**2 - rc**2 - 2 * rc**2 * np.log(r / rc)
    return r**2 + 2 * rc**3 / r - 3 * rc**2


def anoxic_radius(radius, c0, diffusivity, consumption_rate, geometry="slab"):
    """Radius of the anoxic (C = 0) core under zero-order consumption.

    Returns 0 when the construct is fully oxygenated. Beyond the core the
    profile satisfies C = dC/dr = 0 at r = rc, so cells inside the core are
    predicted to be starved (and, in practice, to become necrotic).
    """
    m = _exponent(geometry)
    _require_positive(radius=radius)
    r_crit = critical_dimension(c0, diffusivity, consumption_rate, geometry)
    if radius <= r_crit:
        return 0.0
    if m == 0:
        return radius - r_crit
    # The depletion shape at r = R falls monotonically from R^2 (rc = 0)
    # to 0 (rc = R); find the rc where it equals R_crit^2.
    lo, hi = 0.0, radius
    target = r_crit**2
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if _depletion_shape(radius, mid, m) > target:
            lo = mid
        else:
            hi = mid
        if hi - lo <= 1e-12 * radius:
            break
    return 0.5 * (lo + hi)


def zero_order_profile(r, radius, c0, diffusivity, consumption_rate, geometry="slab"):
    """Analytical steady-state concentration for zero-order consumption.

    Handles an anoxic core when ``radius`` exceeds the critical dimension.
    ``r`` may be a scalar or array with 0 <= r <= radius.
    """
    m = _exponent(geometry)
    r = np.asarray(r, dtype=float)
    if np.any((r < 0) | (r > radius * (1 + 1e-12))):
        raise ValueError("r must lie within [0, radius]")
    q_over_d = consumption_rate / diffusivity / (2 * (m + 1))
    rc = anoxic_radius(radius, c0, diffusivity, consumption_rate, geometry)
    if rc == 0.0:
        return c0 - q_over_d * (radius**2 - r**2)
    outside = r > rc
    safe_r = np.where(outside, r, radius)
    return np.where(outside, q_over_d * _depletion_shape(safe_r, rc, m), 0.0)


def michaelis_menten_profile(radius, c0, diffusivity, vmax, km, geometry="slab",
                             n_points=201, tol=1e-10, max_iter=100):
    """Numerical steady-state profile for Michaelis-Menten consumption.

    Solves D lap(C) = Vmax C / (Km + C) by finite differences and Newton
    iteration. ``vmax`` is the volumetric maximum rate (cell density times
    the per-cell maximum rate), mol/(m^3 s).

    Returns ``(r, c)`` arrays of length ``n_points`` from centre to surface.
    """
    m = _exponent(geometry)
    _require_positive(radius=radius, c0=c0, diffusivity=diffusivity, vmax=vmax, km=km)
    if n_points < 3:
        raise ValueError("n_points must be at least 3")

    n = n_points - 1  # unknowns C_0 .. C_{n-1}; C_n = c0 at the surface
    h = radius / n
    r = np.linspace(0.0, radius, n_points)

    # Discrete Laplacian lap(C) = A @ C + b (finite-volume form).
    a = np.zeros((n, n))
    b = np.zeros(n)
    a[0, 0] = -2.0 * (m + 1) / h**2  # symmetry: C_{-1} = C_1
    if n > 1:
        a[0, 1] = 2.0 * (m + 1) / h**2
    else:
        b[0] = 2.0 * (m + 1) * c0 / h**2
    for i in range(1, n):
        w_minus = ((i - 0.5) / i) ** m / h**2
        w_plus = ((i + 0.5) / i) ** m / h**2
        a[i, i - 1] = w_minus
        a[i, i] = -(w_minus + w_plus)
        if i + 1 < n:
            a[i, i + 1] = w_plus
        else:
            b[i] = w_plus * c0

    def newton(c, km_stage):
        for _ in range(max_iter):
            rate = vmax * c / (km_stage + c)
            d_rate = vmax * km_stage / (km_stage + c) ** 2
            residual = diffusivity * (a @ c + b) - rate
            jacobian = diffusivity * a - np.diag(d_rate)
            step = np.linalg.solve(jacobian, residual)
            c_new = np.clip(c - step, 0.0, c0)
            if np.max(np.abs(c_new - c)) <= tol * c0:
                return c_new
            c = c_new
        raise RuntimeError("Newton iteration did not converge; try more n_points or max_iter")

    # Km << C0 makes the kinetics nearly a step function, which defeats a
    # cold Newton start. Continue from Km = C0 down to the requested Km,
    # warm-starting each stage from the previous solution.
    stages = [km]
    while stages[-1] < c0:
        stages.append(stages[-1] * 10.0)
    c = np.full(n, float(c0))
    for km_stage in reversed(stages):
        c = newton(c, km_stage)

    return r, np.append(c, c0)
