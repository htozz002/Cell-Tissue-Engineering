"""Design calculations for a flow chamber and a perfused porous scaffold.

Run from the repository root:  python examples/perfusion_design.py
"""

import math

from tissue_eng import bioreactor, scaffold
from tissue_eng import constants as k

mu, rho = k.MEDIUM_VISCOSITY, k.MEDIUM_DENSITY

# 1. Parallel-plate flow chamber for endothelial cells.
width, height = 1.0e-2, 250e-6
print("Parallel-plate chamber, 10 mm wide x 250 µm gap")
for target_dyn in (1, 5, 12, 20):
    tau = target_dyn * k.DYN_PER_CM2_TO_PA
    q = bioreactor.parallel_plate_flow_rate(tau, mu, width, height)
    re = bioreactor.parallel_plate_reynolds(q, mu, width, height, rho)
    print(f"  {target_dyn:>3} dyn/cm^2 ({tau:.1f} Pa): {q / k.ML_PER_MIN_TO_M3_PER_S:6.2f} mL/min, Re = {re:.1f}")

# 2. Perfused porous scaffold: architecture -> permeability -> pressure drop and shear.
diameter, length = 6e-3, 3e-3
area = math.pi * (diameter / 2) ** 2
print(f"\nPerfused scaffold, {diameter * 1e3:.0f} mm diameter x {length * 1e3:.0f} mm thick")
print(f"  {'pore d':>7} {'porosity':>8} {'SA/V (1/mm)':>11} {'k (m^2)':>9} {'dP (Pa)':>8} {'shear (mPa)':>11}")
flow = 0.1 * k.ML_PER_MIN_TO_M3_PER_S
for pore in (100e-6, 200e-6, 300e-6):
    for eps in (0.85, 0.95):
        sa_v = scaffold.specific_surface_area(pore, scaffold.relative_density(eps))
        perm = scaffold.kozeny_carman_permeability(eps, sa_v)
        dp = bioreactor.darcy_pressure_drop(flow, mu, length, area, perm)
        tau = bioreactor.perfusion_pore_shear_stress(flow, mu, area, eps, pore)
        print(f"  {pore * 1e6:>5.0f}µm {eps:>8.2f} {sa_v / 1e3:>11.1f} {perm:>9.1e} {dp:>8.2f} {tau * 1e3:>11.2f}")
print("  (flow rate 0.1 mL/min)")
