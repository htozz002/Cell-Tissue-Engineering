# Cell and Tissue Engineering

Course notes, a small tested Python toolkit for the field's standard
quantitative models, and lab reports.

```
├── notes/         Concept notes (cells, growth, transport, scaffolds, bioreactors, translation)
├── src/tissue_eng Python package: the models used in the notes
├── examples/      Runnable design calculations
├── lab-reports/   One folder per lab: report, raw data, analysis, figures
└── tests/         pytest suite
```

## Notes

Start at [`notes/README.md`](notes/README.md). Each note covers the concepts,
key equations and a worked example, then shows the matching code.

## Lab reports

See [`lab-reports/README.md`](lab-reports/README.md). Copy `lab-reports/_template/`
for each new lab; it includes a report outline, a working `analysis.py`
(hemocytometer counts → concentration and viability), and a pre-submission checklist.

## Python toolkit

Requires Python ≥ 3.9 and NumPy.

```bash
pip install -e ".[test]"
pytest
```

| Module | What it covers |
|--------|----------------|
| `tissue_eng.growth` | Exponential, logistic, Gompertz and Monod growth; specific growth rate, doubling time, population doublings |
| `tissue_eng.diffusion` | O2 diffusion–consumption in slabs, cylinders and spheres: critical thickness, Thiele modulus, analytical zero-order profiles with anoxic cores, numerical Michaelis–Menten profiles |
| `tissue_eng.scaffold` | Porosity, Gibson–Ashby foam modulus and strength, specific surface area, Kozeny–Carman permeability, first-order degradation |
| `tissue_eng.bioreactor` | Wall shear in parallel-plate chambers and tubes, Reynolds number, Darcy pressure drop, pore shear in perfused scaffolds |
| `tissue_eng.lab` | Hemocytometer counts, trypan-blue viability, C1V1 = C2V2, seeding calculations, vessel growth areas |
| `tissue_eng.constants` | Typical values at 37 °C (O2 diffusivity and solubility, incubator pO2, medium viscosity) and unit conversions |

**Units:** every module except `lab` uses SI (m, s, Pa, mol/m³). `lab` uses
bench units (mL, cm², cells/mL). `constants` has the conversion factors.

### Quick start

```python
from tissue_eng import diffusion, bioreactor, constants as k

# How thick can a hydrogel with 1e8 cells/mL be before its core goes anoxic?
c0 = diffusion.concentration_from_partial_pressure(k.INCUBATOR_PO2_MMHG)   # mol/m^3
Q = diffusion.volumetric_consumption_rate(1e8 * k.CELLS_PER_ML_TO_CELLS_PER_M3, 5e-17)
diffusion.critical_dimension(c0, k.O2_DIFFUSIVITY_TISSUE, Q, "slab")        # ~3.8e-4 m

# Flow rate for 12 dyn/cm^2 in a 10 mm x 250 µm parallel-plate chamber
q = bioreactor.parallel_plate_flow_rate(12 * k.DYN_PER_CM2_TO_PA, k.MEDIUM_VISCOSITY, 0.01, 250e-6)
q / k.ML_PER_MIN_TO_M3_PER_S                                                # ~9.6 mL/min
```

```bash
python examples/oxygen_limits.py      # critical thickness table + O2 profile
python examples/perfusion_design.py   # flow-chamber settings, scaffold permeability and shear
```

> Constants are typical literature values for estimates and coursework. Use
> measured or cited values for anything you report quantitatively.
