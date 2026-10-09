"""Scaffold architecture, mechanics, permeability and degradation (SI units)."""

import math

import numpy as np


def _require_fraction(name, value):
    if not 0.0 < value < 1.0:
        raise ValueError(f"{name} must be between 0 and 1 (exclusive)")


def porosity_from_density(apparent_density, material_density):
    """Porosity, eps = 1 - rho_apparent / rho_material."""
    if material_density <= 0 or apparent_density <= 0:
        raise ValueError("densities must be positive")
    if apparent_density > material_density:
        raise ValueError("apparent density cannot exceed material density")
    return 1.0 - apparent_density / material_density


def porosity_from_mass(mass, bulk_volume, material_density):
    """Porosity from dry mass and bulk (envelope) volume, eps = 1 - m / (rho_s V)."""
    if bulk_volume <= 0:
        raise ValueError("bulk_volume must be positive")
    return porosity_from_density(mass / bulk_volume, material_density)


def relative_density(porosity):
    """Relative density rho*/rho_s = 1 - eps."""
    _require_fraction("porosity", porosity)
    return 1.0 - porosity


def gibson_ashby_modulus(solid_modulus, rel_density, c=1.0):
    """Elastic modulus of an open-cell foam, E* = C E_s (rho*/rho_s)^2, C ~ 1."""
    _require_fraction("rel_density", rel_density)
    return c * solid_modulus * rel_density**2


def gibson_ashby_collapse_strength(solid_yield_strength, rel_density, c=0.3):
    """Plastic collapse strength of an open-cell foam, sigma* = C sigma_ys (rho*/rho_s)^1.5."""
    _require_fraction("rel_density", rel_density)
    return c * solid_yield_strength * rel_density**1.5


def specific_surface_area(pore_diameter, rel_density):
    """Surface area per unit scaffold volume (1/m) for an open-cell foam.

    Tetrakaidecahedral unit-cell model (O'Brien et al., Biomaterials 2005):
    SA/V = (10.15 / d) sqrt(rho*/rho_s). Smaller pores give more area for
    cell attachment, which is why attachment rises as pore size falls.
    """
    if pore_diameter <= 0:
        raise ValueError("pore_diameter must be positive")
    _require_fraction("rel_density", rel_density)
    return 10.15 / pore_diameter * math.sqrt(rel_density)


def kozeny_carman_permeability(porosity, specific_surface, kozeny_constant=5.0):
    """Darcy permeability (m^2), k = eps^3 / (K S^2).

    ``specific_surface`` is the wetted surface area per unit bulk volume
    (1/m), e.g. from :func:`specific_surface_area`. K ~ 5 for random media.
    """
    _require_fraction("porosity", porosity)
    if specific_surface <= 0:
        raise ValueError("specific_surface must be positive")
    return porosity**3 / (kozeny_constant * specific_surface**2)


def kozeny_carman_permeability_spheres(porosity, particle_diameter):
    """Permeability of a packed bed of spheres, k = eps^3 d^2 / (180 (1 - eps)^2)."""
    _require_fraction("porosity", porosity)
    return porosity**3 * particle_diameter**2 / (180.0 * (1.0 - porosity) ** 2)


def mass_remaining(t, rate_constant, m0=1.0):
    """First-order (bulk hydrolytic) degradation, M(t) = M0 exp(-k t)."""
    if rate_constant < 0:
        raise ValueError("rate_constant must be non-negative")
    return m0 * np.exp(-rate_constant * np.asarray(t, dtype=float))


def degradation_half_life(rate_constant):
    """Half-life for first-order degradation, t_1/2 = ln 2 / k."""
    if rate_constant <= 0:
        raise ValueError("rate_constant must be positive")
    return math.log(2.0) / rate_constant
