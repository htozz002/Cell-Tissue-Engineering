import math

import pytest

from tissue_eng import bioreactor

MU = 0.78e-3


def test_parallel_plate_round_trip():
    w, h = 0.01, 250e-6
    q = bioreactor.parallel_plate_flow_rate(1.2, MU, w, h)
    assert bioreactor.parallel_plate_shear_stress(q, MU, w, h) == pytest.approx(1.2)
    assert q == pytest.approx(1.2 * w * h**2 / (6 * MU))


def test_parallel_plate_flow_is_laminar_for_typical_settings():
    w, h = 0.01, 250e-6
    q = bioreactor.parallel_plate_flow_rate(1.2, MU, w, h)
    assert bioreactor.parallel_plate_reynolds(q, MU, w, h, 1007.0) < 100


def test_tube_shear():
    assert bioreactor.tube_wall_shear_stress(1e-6, 1e-3, 1e-3) == pytest.approx(4e-3 * 1e-6 / (math.pi * 1e-9))


def test_darcy_pressure_drop():
    # mu Q L / (k A) = 1e-3 * 1e-6 * 1e-3 / (1e-10 * 1e-5)
    assert bioreactor.darcy_pressure_drop(1e-6, 1e-3, 1e-3, 1e-5, 1e-10) == pytest.approx(1000.0)


def test_perfusion_pore_shear_is_mpa_scale():
    # 0.1 mL/min through a 6 mm diameter, 90 % porous scaffold with 300 µm pores.
    area = math.pi * (3e-3) ** 2
    tau = bioreactor.perfusion_pore_shear_stress(0.1e-6 / 60, MU, area, 0.9, 300e-6)
    assert 1e-3 < tau < 2e-3
