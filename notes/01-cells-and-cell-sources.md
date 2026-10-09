# 01: Cells and Cell Sources

## Sources by donor

| Source | Donor | Pros | Cons |
|--------|-------|------|------|
| Autologous | The patient | No immune rejection | Needs biopsy, donor-site morbidity, slow expansion, variable quality (age, disease) |
| Allogeneic | Another human | Off-the-shelf, banked, quality-controlled | Immune response, disease-transmission screening |
| Xenogeneic | Another species | Abundant | Rejection, zoonosis risk |
| Syngeneic | Genetically identical twin or inbred animal | No rejection | Rarely available clinically |

## Primary cells and cell lines

- **Primary cells** are isolated directly from tissue. They are physiologically
  relevant but have finite lifespans (replicative senescence, the *Hayflick
  limit* of roughly 40 to 60 doublings for human fibroblasts) and can
  **dedifferentiate** in culture. Chondrocytes in 2D monolayer, for example,
  shift from collagen II to collagen I.
- **Cell lines** are immortalized, either spontaneously or by transformation
  (e.g. TERT or SV40 large T). They are convenient and reproducible but
  karyotypically and functionally abnormal.

## Stem cells

**Potency** runs from totipotent (zygote) to pluripotent (ESC, iPSC),
multipotent (adult or tissue stem cells) and unipotent.

| Type | Origin | Notes |
|------|--------|-------|
| Embryonic stem cells (ESC) | Inner cell mass of blastocyst | Pluripotent; ethical constraints; teratoma risk |
| Induced pluripotent (iPSC) | Somatic cells reprogrammed with Oct4, Sox2, Klf4, c-Myc (Yamanaka factors) | Patient-specific; teratoma and reprogramming-memory concerns |
| Mesenchymal stromal/stem cells (MSC) | Bone marrow, adipose, umbilical cord | Multipotent, immunomodulatory, easy to expand |
| Hematopoietic stem cells (HSC) | Bone marrow, cord blood | Give rise to blood lineages; basis of bone-marrow transplant |

**ISCT minimal criteria for MSCs** (Dominici et al., 2006):
1. Plastic-adherent in standard culture.
2. ≥95 % express CD73, CD90 and CD105; ≤2 % express CD45, CD34, CD14 or
   CD11b, CD79α or CD19, and HLA-DR.
3. Differentiate to osteoblasts, adipocytes and chondroblasts in vitro.

## Cell–matrix interactions

- The **extracellular matrix (ECM)** is collagens (I in bone, skin and
  tendon; II in cartilage; IV in basement membranes), elastin, fibronectin,
  laminin, and proteoglycans/GAGs (aggrecan, hyaluronan).
- **Integrins** (α/β heterodimers) bind ECM ligands such as the **RGD**
  (Arg-Gly-Asp) motif and link to the actin cytoskeleton through **focal
  adhesions**. This is why RGD is grafted onto inert materials like PEG and
  alginate to make them cell-adhesive.
- Cells also respond to soluble signals (growth factors, cytokines) and to
  cell–cell contacts (cadherins, gap junctions).

## Culture basics

- **Medium**: basal medium (DMEM, α-MEM, RPMI) plus serum (often 10 % FBS)
  or defined supplements, antibiotics, and lineage factors (e.g. ascorbate,
  β-glycerophosphate and dexamethasone for osteogenic differentiation; TGF-β3
  for chondrogenic).
- **Environment**: 37 °C, 5 % CO2 (with bicarbonate buffering), humidified.
  Incubator pO2 is about 142 mmHg, well above the 10 to 50 mmHg most tissues
  see, so "normoxia" in culture is hyperoxic relative to the body.
- **Passaging**: detach (trypsin/EDTA), count, re-seed. Record passage number
  and cumulative population doublings (note 02).

## Counting and seeding (code: `tissue_eng.lab`)

Neubauer hemocytometer: each large square holds 0.1 µL, so

$$\text{cells/mL} = \overline{\text{count per square}} \times \text{dilution factor} \times 10^4$$

Trypan blue exclusion: viability = live / (live + dead). A 1:1 mix with
trypan blue has a dilution factor of 2.

```python
from tissue_eng import lab
conc = lab.hemocytometer_concentration([48, 52, 49, 51], dilution_factor=2)  # 1.0e6 cells/mL
n = lab.cells_needed(10_000, lab.GROWTH_AREA_CM2["6-well"], n_vessels=6)     # 570,000 cells
vol = lab.suspension_volume(conc, n)                                         # 0.57 mL
```
