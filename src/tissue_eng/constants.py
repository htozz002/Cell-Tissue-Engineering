"""Typical physical properties at 37 °C and unit conversions (SI).

These values are representative and fine for design estimates and homework.
For anything you report quantitatively, use a value measured for your
system or taken from a cited source.
"""

# --- Unit conversions (multiply to get SI) ---------------------------------
MMHG_TO_PA = 133.322
DYN_PER_CM2_TO_PA = 0.1
UM_TO_M = 1e-6
MM_TO_M = 1e-3
CM_TO_M = 1e-2
ML_TO_M3 = 1e-6
UL_TO_M3 = 1e-9
ML_PER_MIN_TO_M3_PER_S = 1e-6 / 60.0
CELLS_PER_ML_TO_CELLS_PER_M3 = 1e6
UM_CONC_TO_MOL_PER_M3 = 1e-3  # 1 µM = 1e-3 mol/m^3

# --- Oxygen ----------------------------------------------------------------
#: O2 diffusivity in water / culture medium at 37 °C, m^2/s.
O2_DIFFUSIVITY_MEDIUM = 3.0e-9
#: Typical O2 diffusivity in hydrogels and soft tissue, m^2/s (range ~1-2.5e-9).
O2_DIFFUSIVITY_TISSUE = 2.0e-9
#: Henry's-law O2 solubility in culture medium at 37 °C, mol/(m^3 mmHg)
#: (about 1.2-1.3 µM per mmHg).
O2_SOLUBILITY_MEDIUM = 1.3e-3
#: pO2 in a humidified 5 % CO2 incubator at sea level, mmHg
#: (0.95 x 20.9 % O2 x (760 - 47 mmHg water vapour)).
INCUBATOR_PO2_MMHG = 142.0

# --- Culture medium --------------------------------------------------------
#: Dynamic viscosity of culture medium at 37 °C, Pa s. Water is 0.69 mPa s;
#: serum-supplemented medium is typically 0.7-0.95 mPa s.
MEDIUM_VISCOSITY = 0.78e-3
#: Density of culture medium, kg/m^3.
MEDIUM_DENSITY = 1007.0
