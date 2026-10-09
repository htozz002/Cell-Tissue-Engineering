import math

import numpy as np
import pytest

from tissue_eng import diffusion

C0 = 0.185      # mol/m^3, ~142 mmHg
D = 2.0e-9      # m^2/s
Q = 5.0e-3      # mol/(m^3 s)


@pytest.mark.parametrize("geometry,factor", [("slab", 2), ("cylinder", 4), ("sphere", 6)])
def test_critical_dimension(geometry, factor):
    expected = math.sqrt(factor * D * C0 / Q)
    assert diffusion.critical_dimension(C0, D, Q, geometry) == pytest.approx(expected)
    # At the critical size the centre concentration is exactly zero.
    r_crit = expected
    assert diffusion.zero_order_profile(0.0, r_crit, C0, D, Q, geometry) == pytest.approx(0.0, abs=1e-12)
    assert diffusion.thiele_modulus(r_crit, C0, D, Q) ** 2 == pytest.approx(factor)


@pytest.mark.parametrize("geometry", ["slab", "cylinder", "sphere"])
def test_zero_order_profile_boundary_conditions(geometry):
    radius = 3 * diffusion.critical_dimension(C0, D, Q, geometry)
    r = np.linspace(0, radius, 2001)
    c = diffusion.zero_order_profile(r, radius, C0, D, Q, geometry)
    assert c[-1] == pytest.approx(C0)
    assert np.all(c >= 0)
    assert np.all(np.diff(c) >= -1e-15)  # monotone toward the surface
    rc = diffusion.anoxic_radius(radius, C0, D, Q, geometry)
    assert 0 < rc < radius
    assert np.all(c[r <= rc] == 0)


def test_zero_order_profile_satisfies_ode_with_anoxic_core():
    # Check D lap(C) = Q outside the core with a second-difference stencil (sphere).
    radius = 2 * diffusion.critical_dimension(C0, D, Q, "sphere")
    rc = diffusion.anoxic_radius(radius, C0, D, Q, "sphere")
    r = np.linspace(rc + 0.1 * (radius - rc), radius * 0.99, 50)
    h = 1e-7
    c = lambda x: diffusion.zero_order_profile(x, radius, C0, D, Q, "sphere")
    lap = (c(r + h) - 2 * c(r) + c(r - h)) / h**2 + (2 / r) * (c(r + h) - c(r - h)) / (2 * h)
    np.testing.assert_allclose(D * lap, Q, rtol=1e-3)


def test_slab_anoxic_radius_closed_form():
    r_crit = diffusion.critical_dimension(C0, D, Q, "slab")
    assert diffusion.anoxic_radius(2 * r_crit, C0, D, Q, "slab") == pytest.approx(r_crit)
    assert diffusion.anoxic_radius(0.5 * r_crit, C0, D, Q, "slab") == 0.0


def test_michaelis_menten_first_order_limit_matches_cosh():
    # Km >> C0 gives first-order kinetics k = Vmax/Km: C = C0 cosh(phi r/R)/cosh(phi).
    radius, km = 500e-6, 100.0
    vmax = 1.0  # k = 0.01 1/s
    r, c = diffusion.michaelis_menten_profile(radius, C0, D, vmax, km, "slab", n_points=401)
    k = vmax / km
    phi = radius * math.sqrt(k / D)
    exact = C0 * np.cosh(phi * r / radius) / math.cosh(phi)
    np.testing.assert_allclose(c, exact, rtol=2e-3)


@pytest.mark.parametrize("geometry", ["slab", "cylinder", "sphere"])
def test_michaelis_menten_small_km_matches_zero_order(geometry):
    radius = 1.5 * diffusion.critical_dimension(C0, D, Q, geometry)
    r, c = diffusion.michaelis_menten_profile(radius, C0, D, vmax=Q, km=1e-6 * C0,
                                              geometry=geometry, n_points=801)
    exact = diffusion.zero_order_profile(r, radius, C0, D, Q, geometry)
    assert np.max(np.abs(c - exact)) < 1e-3 * C0


def test_concentration_from_partial_pressure():
    assert diffusion.concentration_from_partial_pressure(142.0) == pytest.approx(0.1846)


def test_bad_geometry():
    with pytest.raises(ValueError, match="geometry"):
        diffusion.critical_dimension(C0, D, Q, "cube")
