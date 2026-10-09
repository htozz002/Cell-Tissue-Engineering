"""How thick can a cell-laden construct be before its core runs out of oxygen?

Run from the repository root:  python examples/oxygen_limits.py
"""

from tissue_eng import constants as k
from tissue_eng import diffusion

c0 = diffusion.concentration_from_partial_pressure(k.INCUBATOR_PO2_MMHG)
D = k.O2_DIFFUSIVITY_TISSUE
per_cell_rate = 5e-17  # mol O2 / (cell s); order of magnitude for many mammalian cells

print(f"Surface O2: {k.INCUBATOR_PO2_MMHG:.0f} mmHg -> {c0 * 1e3:.0f} µM")
print(f"D = {D:.1e} m^2/s, q = {per_cell_rate:.0e} mol/(cell s)\n")
print(f"{'cells/mL':>10} | {'slab (1 face)':>13} | {'cyl. radius':>11} | {'sph. radius':>11}")
print("-" * 55)
for cells_per_ml in (1e6, 1e7, 1e8, 1e9):
    Q = diffusion.volumetric_consumption_rate(cells_per_ml * k.CELLS_PER_ML_TO_CELLS_PER_M3, per_cell_rate)
    dims = [diffusion.critical_dimension(c0, D, Q, g) * 1e6 for g in ("slab", "cylinder", "sphere")]
    print(f"{cells_per_ml:>10.0e} | " + " | ".join(f"{x:>8.0f} µm" for x in dims))

# Profile through a 600 µm slab at 1e8 cells/mL, Michaelis-Menten (Km ~ 1 µM).
Q = diffusion.volumetric_consumption_rate(1e8 * k.CELLS_PER_ML_TO_CELLS_PER_M3, per_cell_rate)
thickness = 600e-6
r, c = diffusion.michaelis_menten_profile(thickness, c0, D, vmax=Q, km=1e-3, geometry="slab")
rc = diffusion.anoxic_radius(thickness, c0, D, Q, "slab")
print(f"\n600 µm slab fed from one face, 1e8 cells/mL:")
print(f"  zero-order model predicts an anoxic zone {rc * 1e6:.0f} µm deep from the sealed face")
print("  depth from fed face (µm)   O2 (µM)  [Michaelis-Menten]")
for i in range(0, len(r), 25):
    depth = (thickness - r[i]) * 1e6
    print(f"  {depth:>22.0f}   {c[i] * 1e3:>7.1f}")
