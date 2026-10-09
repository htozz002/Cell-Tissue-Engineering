"""Fluid shear and pressure drop in common bioreactor geometries (SI units).

Flow rates in m^3/s, viscosity in Pa s, lengths in m, stresses in Pa.
(1 Pa = 10 dyn/cm^2; see ``constants`` for conversions.)
"""

import math


def _require_positive(**values):
    for name, value in values.items():
        if value <= 0:
            raise ValueError(f"{name} must be positive")


def parallel_plate_shear_stress(flow_rate, viscosity, width, height):
    """Wall shear stress in a parallel-plate flow chamber, tau = 6 mu Q / (w h^2).

    Assumes laminar, fully developed flow with width >> height.
    """
    _require_positive(viscosity=viscosity, width=width, height=height)
    return 6.0 * viscosity * flow_rate / (width * height**2)


def parallel_plate_flow_rate(shear_stress, viscosity, width, height):
    """Flow rate needed for a target wall shear stress, Q = tau w h^2 / (6 mu)."""
    _require_positive(viscosity=viscosity, width=width, height=height)
    return shear_stress * width * height**2 / (6.0 * viscosity)


def parallel_plate_reynolds(flow_rate, viscosity, width, height, density):
    """Reynolds number using the hydraulic diameter D_h = 2 w h / (w + h).

    Flow is laminar for Re below roughly 2000; flow chambers usually run far below.
    """
    _require_positive(viscosity=viscosity, width=width, height=height, density=density)
    mean_velocity = flow_rate / (width * height)
    hydraulic_diameter = 2.0 * width * height / (width + height)
    return density * mean_velocity * hydraulic_diameter / viscosity


def tube_wall_shear_stress(flow_rate, viscosity, radius):
    """Poiseuille wall shear stress in a tube, tau = 4 mu Q / (pi R^3)."""
    _require_positive(viscosity=viscosity, radius=radius)
    return 4.0 * viscosity * flow_rate / (math.pi * radius**3)


def darcy_pressure_drop(flow_rate, viscosity, length, area, permeability):
    """Pressure drop across a perfused porous scaffold, dP = mu Q L / (k A)."""
    _require_positive(viscosity=viscosity, length=length, area=area,
                      permeability=permeability)
    return viscosity * flow_rate * length / (permeability * area)


def perfusion_pore_shear_stress(flow_rate, viscosity, area, porosity, pore_diameter):
    """Mean wall shear stress in a perfused scaffold, cylindrical-pore model.

    Pores are treated as straight cylinders of diameter d carrying Poiseuille
    flow at the interstitial velocity u = Q / (A eps), so tau = 8 mu u / d.
    This is an order-of-magnitude estimate; local shear in real scaffolds
    varies widely (CFD on micro-CT geometry is the rigorous route).
    """
    _require_positive(viscosity=viscosity, area=area, pore_diameter=pore_diameter)
    if not 0.0 < porosity < 1.0:
        raise ValueError("porosity must be between 0 and 1 (exclusive)")
    interstitial_velocity = flow_rate / (area * porosity)
    return 8.0 * viscosity * interstitial_velocity / pore_diameter
