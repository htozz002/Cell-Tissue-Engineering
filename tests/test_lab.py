import pytest

from tissue_eng import lab


def test_hemocytometer_with_trypan_blue():
    # Mean of 50 cells per square, 1:1 with trypan blue -> 1.0e6 cells/mL.
    assert lab.hemocytometer_concentration([48, 52, 49, 51], dilution_factor=2) == pytest.approx(1.0e6)


def test_viability():
    assert lab.viability(190, 10) == pytest.approx(0.95)
    with pytest.raises(ValueError):
        lab.viability(0, 0)


def test_seeding_workflow():
    n = lab.cells_needed(10_000, lab.GROWTH_AREA_CM2["6-well"], n_vessels=6)
    assert n == pytest.approx(570_000)
    assert lab.suspension_volume(1.0e6, n) == pytest.approx(0.57)


def test_dilution():
    assert lab.dilution_volume(1.0e6, 2.5e5, 4.0) == pytest.approx(1.0)
    with pytest.raises(ValueError):
        lab.dilution_volume(1.0, 2.0, 1.0)
