"""Quantitative models and bench calculations for cell and tissue engineering.

Modules
-------
growth      Cell population growth kinetics (exponential, logistic, Gompertz, Monod).
diffusion   Steady-state oxygen/nutrient transport with cellular consumption.
scaffold    Scaffold porosity, foam mechanics, surface area, permeability, degradation.
bioreactor  Shear stress and pressure drop in flow chambers, tubes and perfused scaffolds.
lab         Bench calculations: hemocytometer counts, viability, dilutions, seeding.
constants   Typical physical properties at 37 °C and unit conversions.

Every module except ``lab`` uses SI units. ``lab`` uses bench units
(mL, cm^2, cells/mL) because that is how protocols are written.
"""

from . import bioreactor, constants, diffusion, growth, lab, scaffold

__all__ = ["bioreactor", "constants", "diffusion", "growth", "lab", "scaffold"]
__version__ = "0.1.0"
