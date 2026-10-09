# 04: Biomaterials and Scaffolds

## What a scaffold must do

- **Biocompatible**: no toxic leachables or degradation products; an acceptable host response.
- **Biodegradable** at a rate matched to new-tissue formation (or permanent, by design).
- **Porous and interconnected**, so cells can migrate in and nutrients and waste can move (note 03).
- **Mechanically appropriate**: survives handling and implantation, and presents the right stiffness to cells.
- **Bioactive surface**: adsorbs proteins and presents adhesion ligands (e.g. RGD).
- **Manufacturable and sterilizable**, reproducibly and at scale.

## Material classes

| Class | Examples | Notes |
|-------|----------|-------|
| Natural polymers | Collagen, gelatin, fibrin, hyaluronic acid, alginate, chitosan, silk fibroin, decellularized ECM | Intrinsic bioactivity (except alginate); batch variability; weak |
| Synthetic polymers | PGA, PLA, PLGA, PCL, PEG | Tunable and reproducible; bioinert unless functionalized |
| Ceramics | Hydroxyapatite (HA), β-TCP, bioactive glass (45S5) | Osteoconductive, stiff, brittle. HA resorbs slowly, β-TCP faster |
| Composites | PCL–HA, collagen–HA | Combine toughness with bioactivity |

### Biodegradable polyesters

These degrade by **bulk hydrolysis** of ester bonds into lactic and/or
glycolic acid. The acidic products can lower local pH and autocatalyse
degradation inside thick devices. Approximate full-resorption times
(Middleton & Tipton, *Biomaterials* 2000):

| Polymer | Resorption time | Comment |
|---------|-----------------|---------|
| PGA | 6–12 months | Crystalline, stiff; loses strength within weeks |
| PLGA 50:50 | 1–2 months | Amorphous, so it degrades faster than either homopolymer |
| PLGA 85:15 | 5–6 months | More lactide means slower |
| PLLA | > 24 months | Semicrystalline |
| PCL | > 24 months | Low Tm (≈60 °C), so it is easy to melt-print |

First-order mass loss is a common first model: $M(t) = M_0 e^{-kt}$, $t_{1/2} = \ln 2 / k$.

### Hydrogels

These are water-swollen polymer networks with tissue-like water content and
stiffness, and they can encapsulate cells. Crosslinking options:

- **Physical**: ionic (alginate + Ca²⁺), thermal (collagen self-assembly at 37 °C, gelatin, agarose), self-assembling peptides.
- **Chemical**: photo-crosslinked methacrylates (GelMA, PEGDA), enzymatic (thrombin → fibrin; transglutaminase), click chemistry.

## Fabrication methods

| Method | Typical feature size | Trade-offs |
|--------|----------------------|-----------|
| Solvent casting / particulate leaching | Pores set by porogen (e.g. NaCl) size | Simple; limited thickness and interconnectivity; residual solvent |
| Gas foaming (CO2) | 100 µm-scale pores | Solvent-free; poor interconnectivity |
| Freeze-drying | 20–300 µm pores | Pore size set by freezing rate; standard for collagen-GAG |
| Electrospinning | Nano- to micro-scale fibres | ECM-like fibres; small pores limit cell infiltration |
| 3D printing / bioprinting | ~100 µm and up | Designed architecture; extrusion, inkjet, stereolithography; cell-laden bioinks |
| Decellularization | Native architecture | Detergents (SDS, Triton X-100) and enzymes; aim for < 50 ng dsDNA/mg dry ECM (Crapo et al., 2011) |

## Architecture and mechanics (code: `tissue_eng.scaffold`)

**Porosity** from density or mass:

$$\varepsilon = 1 - \frac{\rho^*}{\rho_s} = 1 - \frac{m}{\rho_s V}$$

**Gibson–Ashby open-cell foam** scaling, with relative density $\rho^*/\rho_s = 1-\varepsilon$:

$$E^* \approx E_s\left(\frac{\rho^*}{\rho_s}\right)^2, \qquad \sigma^*_{pl} \approx 0.3\,\sigma_{ys}\left(\frac{\rho^*}{\rho_s}\right)^{3/2}$$

A 90 %-porous foam of a 1.5 GPa polymer gives only ≈15 MPa. Porosity costs
stiffness quadratically.

**Specific surface area** (tetrakaidecahedral cell; O'Brien et al., 2005),
with $d$ the mean pore diameter:

$$\frac{SA}{V} = \frac{10.15}{d}\sqrt{\frac{\rho^*}{\rho_s}}$$

Smaller pores give more attachment area. Larger pores help migration,
infiltration and vascularization. That conflict is why pore size is a design
trade-off.

**Permeability** (Kozeny–Carman), with $S$ = surface area per bulk volume and $K \approx 5$:

$$k = \frac{\varepsilon^3}{K S^2}$$

### Pore-size rules of thumb

- Bone: minimum ≈100 µm; > 300 µm favours new bone and capillary formation (Karageorgiou & Kaplan, 2005).
- Dermal regeneration (collagen-GAG templates): ≈20–125 µm.
- Fibroblast and osteoblast attachment rises as pore size falls, tracking SA/V.

### Stiffness of tissues (order of magnitude)

| Tissue | Modulus |
|--------|---------|
| Brain | ≈0.1–1 kPa |
| Skeletal muscle | ≈10 kPa |
| Articular cartilage (compressive, aggregate) | ≈0.5–1 MPa |
| Trabecular bone (apparent) | ≈0.05–2 GPa, strongly density-dependent |
| Tendon (tensile) | ≈1 GPa |
| Cortical bone | ≈15–20 GPa |

Substrate stiffness alone can bias MSC fate: soft (≈0.1–1 kPa) toward
neurogenic, ≈8–17 kPa toward myogenic, ≈25–40 kPa toward osteogenic
(Engler et al., *Cell* 2006).

```python
from tissue_eng import scaffold
eps = scaffold.porosity_from_density(130, 1300)             # 0.90
rd = scaffold.relative_density(eps)
scaffold.gibson_ashby_modulus(1.5e9, rd)                    # 1.5e7 Pa
sa_v = scaffold.specific_surface_area(100e-6, rd)           # 1/m
scaffold.kozeny_carman_permeability(eps, sa_v)              # m^2
```

## Sterilization

Autoclaving melts or hydrolyses most degradable polymers and denatures
proteins. Gamma irradiation causes chain scission in polyesters (faster
degradation later). Ethylene oxide leaves residues that need aeration.
Ethanol and UV disinfect but do not sterilize. Choose the method early,
because it changes the material.
