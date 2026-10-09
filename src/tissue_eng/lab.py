"""Bench calculations for cell culture, in the units protocols use.

Volumes in mL, areas in cm^2, concentrations in cells/mL unless noted.
"""

import numpy as np

#: Volume over one large (1 mm x 1 mm) square of a Neubauer hemocytometer
#: with a 0.1 mm chamber depth: 0.1 µL = 1e-4 mL.
NEUBAUER_SQUARE_VOLUME_ML = 1e-4

#: Nominal growth areas (cm^2) for common Corning vessels. Other brands
#: differ slightly; check the datasheet for the vessel you actually used.
GROWTH_AREA_CM2 = {
    "96-well": 0.32,
    "48-well": 0.95,
    "24-well": 1.9,
    "12-well": 3.8,
    "6-well": 9.5,
    "35 mm dish": 9.6,
    "60 mm dish": 21.0,
    "100 mm dish": 55.0,
    "T-25": 25.0,
    "T-75": 75.0,
    "T-175": 175.0,
}


def hemocytometer_concentration(counts, dilution_factor=1.0,
                                square_volume_ml=NEUBAUER_SQUARE_VOLUME_ML):
    """Cell concentration (cells/mL) from hemocytometer square counts.

    ``counts`` is the number of cells in each large square counted.
    ``dilution_factor`` is total volume / sample volume, e.g. 2 for a
    1:1 mix with trypan blue.
    """
    counts = np.asarray(counts, dtype=float)
    if counts.size == 0:
        raise ValueError("counts must contain at least one square")
    if np.any(counts < 0):
        raise ValueError("counts must be non-negative")
    if dilution_factor < 1:
        raise ValueError("dilution_factor must be >= 1")
    return float(counts.mean() * dilution_factor / square_volume_ml)


def viability(live, dead):
    """Fraction of viable cells, live / (live + dead), e.g. from trypan blue exclusion."""
    if live < 0 or dead < 0:
        raise ValueError("counts must be non-negative")
    if live + dead == 0:
        raise ValueError("no cells counted")
    return live / (live + dead)


def dilution_volume(stock_concentration, target_concentration, target_volume):
    """Stock volume needed for a dilution, V1 = C2 V2 / C1.

    Units only need to match (same concentration unit for C1 and C2;
    V1 comes back in the unit of V2). Make up to V2 with diluent.
    """
    if stock_concentration <= 0:
        raise ValueError("stock_concentration must be positive")
    if target_concentration > stock_concentration:
        raise ValueError("target concentration exceeds stock; cannot dilute up")
    return target_concentration * target_volume / stock_concentration


def cells_needed(density_per_cm2, area_cm2, n_vessels=1):
    """Total cells to seed ``n_vessels`` at ``density_per_cm2``."""
    return density_per_cm2 * area_cm2 * n_vessels


def suspension_volume(cell_concentration, n_cells):
    """Volume of cell suspension (mL) that contains ``n_cells``."""
    if cell_concentration <= 0:
        raise ValueError("cell_concentration must be positive")
    return n_cells / cell_concentration
