# Lab Reports

One folder per lab. Each holds the write-up, its raw data, the analysis code,
and the figures it produces, so every number in a report can be traced back
to a file.

## Starting a new report

```bash
cp -r lab-reports/_template lab-reports/lab-01-cell-counting
```

Name folders `lab-NN-short-topic` so they sort in order. Then:

1. Put instrument exports, count sheets and plate-reader files in `data/raw/`
   **unchanged**. Never edit raw data; correct or clean it in code.
2. Write analysis in `analysis.py`. It reads `data/raw/`, writes cleaned
   tables to `data/processed/` and plots to `figures/`.
3. Write the report in `report.md`, referencing figures in `figures/`.

```
lab-NN-topic/
├── report.md          # the write-up
├── analysis.py        # raw data -> processed data + figures
├── data/
│   ├── raw/           # untouched instrument/hand-count exports
│   └── processed/     # outputs of analysis.py
└── figures/           # plots used in report.md
```

Run the analysis from the repository root, after installing the package
once with `pip install -e .`:

```bash
python lab-reports/lab-01-cell-counting/analysis.py
```

## Checklist before submitting

- [ ] Every value has **units** and sensible **significant figures**.
- [ ] Replicates are stated as **n = ?** and say whether they are technical or biological.
- [ ] Error bars say what they are (SD or SEM) in the caption.
- [ ] Statistical tests are named, with α and the exact p-values.
- [ ] Each figure has a numbered caption that stands on its own.
- [ ] Methods list cell type and passage, seeding density, medium, timepoints, and instrument settings.
- [ ] Deviations from the protocol are reported.
- [ ] The discussion compares results to literature values, with citations.
- [ ] `analysis.py` regenerates every figure from `data/raw/`.

## Useful package functions

| Task | Function |
|------|----------|
| Hemocytometer count → cells/mL | `tissue_eng.lab.hemocytometer_concentration` |
| Trypan blue viability | `tissue_eng.lab.viability` |
| C1V1 = C2V2 | `tissue_eng.lab.dilution_volume` |
| Cells and volume to seed | `tissue_eng.lab.cells_needed`, `tissue_eng.lab.suspension_volume` |
| Vessel growth areas | `tissue_eng.lab.GROWTH_AREA_CM2` |
| Growth rate, doubling time, PD | `tissue_eng.growth` |
| Scaffold porosity, modulus | `tissue_eng.scaffold` |
| Flow-chamber shear | `tissue_eng.bioreactor` |
