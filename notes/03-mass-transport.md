# 03: Mass Transport and Oxygen Limits

Cells more than about **100 to 200 µm** from a capillary become hypoxic. That
single fact limits the size of every avascular engineered tissue, and it falls
straight out of a diffusion–reaction balance.

## Governing equation

Fick's first law: $J = -D\,\dfrac{dC}{dx}$. At steady state, diffusion
balances consumption:

$$D\,\frac{1}{r^m}\frac{d}{dr}\!\left(r^m \frac{dC}{dr}\right) = R(C), \qquad m = 0,1,2 \text{ for slab, cylinder, sphere}$$

Boundary conditions: $C = C_0$ at the exposed surface ($r = R$) and no flux
at the centre ($dC/dr = 0$ at $r = 0$).

## Consumption kinetics

- **Zero order**: $R = Q = \rho_{cell}\, q_{cell}$ (cell density × per-cell rate).
  This is a good approximation while $C \gg K_m$.
- **Michaelis–Menten**: $R = V_{max} C / (K_m + C)$. For O2, $K_m$ is roughly
  1 µM (about 0.5 to 1 mmHg), far below the surface concentration, so
  zero-order is usually a good first estimate.

## Zero-order solution and critical thickness

$$C(r) = C_0 - \frac{Q}{2(m+1)D}\,(R^2 - r^2)$$

Setting $C(0) = 0$ gives the **critical dimension**:

$$R_{crit} = \sqrt{\frac{2(m+1)\,D\,C_0}{Q}} \;=\; \sqrt{\frac{2DC_0}{Q}}\ \text{(slab)},\ \sqrt{\frac{4DC_0}{Q}}\ \text{(cylinder)},\ \sqrt{\frac{6DC_0}{Q}}\ \text{(sphere)}$$

For a slab fed from one face, R is the full thickness; for one fed from both
faces, R is the half-thickness. The dimensionless form is the Thiele modulus
$\phi^2 = QR^2/(DC_0)$; the core stays oxygenated while $\phi^2 \le 2(m+1)$.

Beyond $R_{crit}$ an **anoxic core** forms. The code solves for its radius
exactly (`diffusion.anoxic_radius`).

## Typical values

| Quantity | Value |
|----------|-------|
| O2 solubility in medium, 37 °C | ≈1.3 µM/mmHg |
| Incubator pO2 (5 % CO2, humidified) | ≈142 mmHg, so $C_0 \approx 185$ µM |
| D(O2) in medium, 37 °C | ≈3 × 10⁻⁹ m²/s |
| D(O2) in hydrogel or tissue | ≈1–2.5 × 10⁻⁹ m²/s |
| Per-cell O2 consumption | ≈10⁻¹⁸ (chondrocytes) to ≈10⁻¹⁶ (hepatocytes) mol/cell/s |
| Cell density in native tissue | ≈10⁸–10⁹ cells/mL |

## Worked example

Hydrogel with $D = 2\times10^{-9}$ m²/s, surface at 142 mmHg ($C_0 = 0.185$ mol/m³),
$q = 5\times10^{-17}$ mol/cell/s:

| Cells/mL | Q (mol m⁻³ s⁻¹) | Slab, fed one face | Cylinder radius | Sphere radius |
|----------|-----------------|--------------------|-----------------|---------------|
| 10⁷ | 5 × 10⁻⁴ | 1.2 mm | 1.7 mm | 2.1 mm |
| 10⁸ | 5 × 10⁻³ | 384 µm | 543 µm | 666 µm |
| 10⁹ | 5 × 10⁻² | 122 µm | 172 µm | 210 µm |

At tissue-like densities, the limit is on the order of 100 to 200 µm, the
same as the intercapillary distance in vivo. Note that $R_{crit} \propto
1/\sqrt{\rho_{cell}}$: ten times the cells gives about one-third the thickness.

```python
from tissue_eng import diffusion
c0 = diffusion.concentration_from_partial_pressure(142)        # 0.185 mol/m^3
Q = diffusion.volumetric_consumption_rate(1e14, 5e-17)         # 1e8 cells/mL
diffusion.critical_dimension(c0, 2e-9, Q, "slab")              # 3.84e-4 m
r, c = diffusion.michaelis_menten_profile(600e-6, c0, 2e-9, vmax=Q, km=1e-3)
```

See `examples/oxygen_limits.py` for a full table and profile.

## Convection versus diffusion

The Péclet number $Pe = UL/D$ compares convective to diffusive transport.
For O2 ($D \approx 3\times10^{-9}$ m²/s) over 1 mm, even 10 µm/s of
perfusion gives $Pe \approx 3$, which is why perfusion bioreactors (note 05)
support thicker constructs than static culture.

## Strategies to beat the diffusion limit

- Perfusion culture and porous, interconnected scaffolds.
- Thin or sheet-based constructs (cell sheets stacked with vascularization between layers).
- Lower initial cell density, letting tissue mature as vessels form.
- Engineered channels (sacrificial templates, bioprinting) and prevascularization (note 06).
- O2-carrying or O2-generating materials (perfluorocarbons, calcium peroxide).
