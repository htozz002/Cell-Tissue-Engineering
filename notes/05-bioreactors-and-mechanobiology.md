# 05: Bioreactors and Mechanobiology

## Why bioreactors?

1. **Mass transfer**: move O2 and nutrients into, and waste out of, thick constructs.
2. **Uniform seeding**: dynamic seeding fills scaffold interiors better than pipetting on top.
3. **Mechanical conditioning**: tissues mature functionally under the loads they see in vivo.
4. **Control and monitoring**: pH, O2, temperature and sterility, reproducibly.

## Common designs

| Bioreactor | Mechanism | Used for |
|------------|-----------|----------|
| Spinner flask | Stirred medium; convection around scaffolds | Better seeding and mixing; non-uniform, partly turbulent shear |
| Rotating wall vessel | Constructs in continuous free fall; low shear | Cartilage, 3D aggregates ("simulated microgravity") |
| Perfusion | Medium forced through the scaffold pores | Bone and thick constructs; best interior mass transfer |
| Parallel-plate flow chamber | Defined laminar shear on a monolayer | Endothelial and bone-cell mechanobiology |
| Compression / hydrostatic pressure | Cyclic loading, often ≈1 Hz | Cartilage |
| Tensile / cyclic stretch | Uniaxial or biaxial strain | Tendon, ligament, muscle, heart, vessels |
| Pulsatile flow | Cyclic pressure and flow in tubes | Vascular grafts, heart valves |

## Shear stress formulas (code: `tissue_eng.bioreactor`)

All assume laminar, fully developed Newtonian flow. 1 Pa = 10 dyn/cm².

**Parallel-plate chamber** (width $w$, gap $h$, $w \gg h$):

$$\tau_w = \frac{6\mu Q}{w h^2}$$

**Tube** (Poiseuille, radius $R$):

$$\tau_w = \frac{4\mu Q}{\pi R^3}$$

**Perfused scaffold**: Darcy's law gives the pressure drop. Treating pores as
cylinders of diameter $d$ at interstitial velocity $u = Q/(A\varepsilon)$
gives an order-of-magnitude shear estimate:

$$\Delta P = \frac{\mu Q L}{k A}, \qquad \tau \approx \frac{8\mu u}{d}$$

Check laminar flow with $Re = \rho U D_h / \mu$. Flow chambers usually run at
$Re < 100$.

### Worked example

Flow chamber with $w = 10$ mm and $h = 250$ µm, medium $\mu = 0.78$ mPa·s.
To apply 12 dyn/cm² (1.2 Pa), the arterial-like shear used in many
endothelial studies:

$$Q = \frac{\tau w h^2}{6\mu} = \frac{1.2 \times 0.01 \times (250\times10^{-6})^2}{6 \times 0.78\times10^{-3}} = 1.6\times10^{-7}\ \text{m}^3/\text{s} \approx 9.6\ \text{mL/min}$$

with $Re \approx 40$ (laminar). Because $\tau \propto 1/h^2$, a 10 % error in
gap height is about a 20 % error in shear. Measure the gasket.

```python
from tissue_eng import bioreactor
q = bioreactor.parallel_plate_flow_rate(1.2, 0.78e-3, 0.01, 250e-6)   # 1.6e-7 m^3/s
```

`examples/perfusion_design.py` runs these numbers for several shear targets
and scaffold architectures.

## Physiological reference values

| Setting | Shear stress |
|---------|--------------|
| Arterial endothelium | ≈1–7 Pa (10–70 dyn/cm²) |
| Venous endothelium | ≈0.1–0.6 Pa (1–6 dyn/cm²) |
| Osteocyte lacunar–canalicular flow (estimated) | ≈0.8–3 Pa |
| Typical perfusion bioreactor (mean, in pores) | ≈0.1 to tens of mPa |

Medium viscosity is lower than blood's (≈3–4 mPa·s), so matching
physiological shear in vitro needs higher flow rates or narrower gaps.

## Mechanotransduction: how cells sense force

- **Integrins and focal adhesions** transmit force from ECM to the actin cytoskeleton.
- **Cytoskeletal tension** (RhoA/ROCK, myosin II) feeds back on adhesion and shape.
- **YAP/TAZ** move into the nucleus on stiff substrates and in spread cells,
  driving proliferation and osteogenic programmes.
- **Stretch-activated ion channels** (Piezo1/2, TRP channels) open under membrane tension.
- **Primary cilia** and the **endothelial glycocalyx** sense fluid flow.
- **Nucleus**: the LINC complex and lamin A/C transmit force to chromatin.

Examples: endothelial cells align with laminar flow and upregulate KLF2 and
eNOS. Bone adapts its mass to loading (Wolff's law; Frost's mechanostat).
Cartilage under dynamic compression makes more GAG and collagen II than
statically cultured controls.
