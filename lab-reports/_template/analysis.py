"""Analysis for Lab NN.

Reads data/raw/, writes data/processed/ and figures/. Run from the repo root:

    python lab-reports/lab-NN-topic/analysis.py

The example below processes hemocytometer counts. Replace it with this lab's
analysis, and delete data/raw/example_counts.csv once you add real data.
"""

import csv
from pathlib import Path

from tissue_eng import lab

HERE = Path(__file__).resolve().parent
RAW = HERE / "data" / "raw"
PROCESSED = HERE / "data" / "processed"
FIGURES = HERE / "figures"


def analyze_counts(path):
    """One row per sample: live/dead counts in four squares, trypan blue dilution."""
    results = []
    with open(path, newline="") as f:
        for row in csv.DictReader(f):
            live = [int(row[f"live_sq{i}"]) for i in range(1, 5)]
            dead = [int(row[f"dead_sq{i}"]) for i in range(1, 5)]
            dilution = float(row["dilution_factor"])
            results.append({
                "sample": row["sample"],
                "live_cells_per_mL": lab.hemocytometer_concentration(live, dilution),
                "viability_percent": 100 * lab.viability(sum(live), sum(dead)),
            })
    return results


def main():
    PROCESSED.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)

    results = analyze_counts(RAW / "example_counts.csv")
    out = PROCESSED / "cell_counts.csv"
    with open(out, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(results[0]))
        writer.writeheader()
        writer.writerows(results)

    for r in results:
        print(f"{r['sample']:>10}: {r['live_cells_per_mL']:.2e} live cells/mL, "
              f"{r['viability_percent']:.1f} % viable")
    print(f"Wrote {out.relative_to(HERE)}")


if __name__ == "__main__":
    main()
