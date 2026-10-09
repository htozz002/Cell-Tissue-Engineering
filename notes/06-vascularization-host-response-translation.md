# 06: Vascularization, Host Response and Translation

## Vascularization

The O2 diffusion limit (note 03) means anything thicker than a few hundred
microns needs a blood supply quickly after implantation.

- **Angiogenesis**: sprouting of new vessels from existing ones (the host's
  route into an implant).
- **Vasculogenesis**: de novo vessel formation from endothelial progenitors.
- **Inosculation**: an engineered vascular network connecting to host vessels.

### Strategies

| Strategy | Idea |
|----------|------|
| Growth-factor delivery | VEGF (sprouting), FGF-2, PDGF-BB and angiopoietin-1 (pericyte recruitment, stabilization). Sequential VEGF then PDGF makes more mature vessels (Richardson et al., *Nat Biotechnol* 2001) |
| Co-culture prevascularization | Endothelial cells (e.g. HUVECs) with supporting cells (MSCs, fibroblasts) self-assemble into networks that can inosculate |
| Engineered channels | Sacrificial templates (e.g. carbohydrate glass, Miller et al. 2012), bioprinting, microfluidics: perfusable from day 0 |
| In vivo prevascularization | Grow the construct in a vascular bed first (arteriovenous loop, muscle or omental flap), then transfer |
| Porous architecture | Large, interconnected pores (> 300 µm in bone) allow ingrowth |

## Host response to implanted biomaterials

Typical sequence (Anderson et al., 2008):

1. **Protein adsorption** within seconds. The adsorbed layer, not the bare
   material, is what cells see (Vroman effect: proteins exchange over time).
2. **Provisional matrix** (fibrin clot) forms.
3. **Acute inflammation** (hours to days): neutrophils.
4. **Chronic inflammation**: monocytes, macrophages, lymphocytes.
5. **Granulation tissue**: fibroblasts and new capillaries.
6. **Foreign-body reaction**: macrophages fuse into foreign-body giant cells on non-degradable surfaces.
7. **Fibrous encapsulation**: an avascular collagen capsule walls the implant off.

**Macrophage polarization** is often summarized as M1 (pro-inflammatory)
and M2 (pro-remodelling). An early M1 phase followed by a timely shift to M2
is associated with constructive remodelling rather than scarring, for
example with decellularized ECM scaffolds. Allogeneic cells add an adaptive
immune response (HLA mismatch) on top of this.

## Translation and regulation

Getting from bench to bedside means scale-up under **GMP**, defined release
criteria (identity, purity, potency, sterility), stability and shelf life,
logistics (living products ship cold or cryopreserved), and cost.

### United States (FDA)

- **HCT/Ps** (human cells, tissues, and cellular and tissue-based products),
  21 CFR Part 1271:
  - **"361" products** are regulated only for disease transmission, with no
    premarket approval, if they are *minimally manipulated*, for *homologous
    use*, not combined with another article (limited exceptions), and either
    have no systemic effect or are for autologous or close-relative use.
  - **"351" products** are everything else. They are regulated as biologics
    by CBER: IND, clinical trials, then a Biologics License Application (BLA).
- **Acellular scaffolds** are usually medical devices (CDRH): 510(k), De Novo or PMA.
- **Combination products** (cells + scaffold, drug + device) are assigned to
  a lead centre by their *primary mode of action*.
- The **RMAT** designation (21st Century Cures Act, 2016) gives expedited
  review to qualifying regenerative-medicine therapies.

### European Union

Cell- and tissue-based therapies fall under the **Advanced Therapy Medicinal
Products (ATMP)** regulation (EC 1394/2007). "Tissue-engineered products"
are one ATMP class, reviewed by the EMA's Committee for Advanced Therapies.

### Approved examples (US)

| Product | Year | What it is |
|---------|------|------------|
| Integra Dermal Regeneration Template | 1996 | Acellular collagen-GAG dermal template with silicone layer (device) |
| Carticel | 1997 | Autologous cultured chondrocytes (biologic) |
| Apligraf | 1998 | Bilayered living skin: neonatal keratinocytes and fibroblasts in bovine collagen |
| Dermagraft | 2001 | Neonatal fibroblasts on a bioabsorbable polyglactin mesh |
| MACI | 2016 | Autologous chondrocytes on a porcine collagen membrane |

## Review questions

1. Why does a construct that is fine at 10⁷ cells/mL fail at 10⁹ cells/mL? Compute both critical thicknesses.
2. You need a 5 mm bone scaffold with ≥ 10 MPa modulus from a 2 GPa polymer. What is the maximum porosity (Gibson–Ashby)? Is it compatible with > 300 µm pores and good permeability?
3. Design a flow chamber (choose w and h) that applies 1.5 Pa at under 10 mL/min with Re < 100.
4. Is an autologous fat graft, mechanically processed and reinjected as fat, likely a 361 or a 351 product? What about the same cells expanded in culture and injected into a knee?
