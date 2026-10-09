import math

import pytest

from tissue_eng import scaffold


def test_porosity():
    # PLGA (~1.3 g/cm^3) foam at 0.13 g/cm^3 is 90 % porous.
    assert scaffold.porosity_from_density(130.0, 1300.0) == pytest.approx(0.9)
    assert scaffold.porosity_from_mass(0.013e-3, 1e-7, 1300.0) == pytest.approx(0.9)
    with pytest.raises(ValueError):
        scaffold.porosity_from_density(2000.0, 1300.0)


def test_gibson_ashby():
    assert scaffold.gibson_ashby_modulus(1e9, 0.1) == pytest.approx(1e7)
    assert scaffold.gibson_ashby_collapse_strength(50e6, 0.1) == pytest.approx(0.3 * 50e6 * 0.1**1.5)


def test_specific_surface_area_scales_inverse_with_pore_size():
    a = scaffold.specific_surface_area(100e-6, 0.05)
    b = scaffold.specific_surface_area(200e-6, 0.05)
    assert a == pytest.approx(2 * b)
    assert a == pytest.approx(10.15 / 100e-6 * math.sqrt(0.05))


def test_kozeny_carman_general_form_matches_spheres():
    eps, d = 0.4, 1e-4
    # For spheres, surface per bulk volume is 6 (1 - eps) / d.
    s = 6 * (1 - eps) / d
    assert scaffold.kozeny_carman_permeability(eps, s) == pytest.approx(
        scaffold.kozeny_carman_permeability_spheres(eps, d))


def test_degradation_half_life():
    k = 0.01
    t_half = scaffold.degradation_half_life(k)
    assert scaffold.mass_remaining(t_half, k, m0=10.0) == pytest.approx(5.0)
