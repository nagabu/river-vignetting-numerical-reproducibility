# RIVER — numerical reproducibility package

**Scope:** RIVER's published numerical segmentation results, Figure 3 numeric
source and paired-seed statistical comparison checks; not end-to-end model
retraining. The manuscript associated with this package is the November 2,
2026 version: *Vignetting-Associated Differences in Dermoscopic Lesion
Segmentation Performance*, by Takashi Nagaoka.

## Contents

- `release_data/139_seed_level_6level_primary5_metrics.csv` — 180 records
  (10 repeated seeds × 6 training exposures × 3 evaluation subgroups), five
  metrics (Dice, IoU, Boundary F1, HD95, ASSD).
- `release_data/139_exposure_curve_6level_primary5_summary.csv` — 90 Figure 3
  records, means/SDs/95% confidence intervals.
- `release_data/133_all_primary_statistical_comparisons.csv` — 65 prespecified
  seed-level statistical comparison records.
- `release_data/140_upper_bound_50_to_64p3_comparison.csv` — 15 additional
  comparisons of 50% with 64.3% exposure.
- `tools/reproduce_figure3_data.py` and
  `tools/reproduce_primary_statistics.py` — independent numerical recalculation
  and equality checks against the released source tables.
- `RELEASE_MANIFEST.json` — checksums and file inventory.

**Not included:** original dermoscopic images, segmentation masks, image-level
predictions, original image identifiers, train/validation/test memberships,
patient or lesion identifiers, original immutable split-lock materials, model
checkpoints, or a validated end-to-end training pipeline. This release does not
independently reproduce the original image preprocessing or model training.

## Recalculate the numerical results

With Python and SciPy installed, from the extracted package root:

```bash
python -m pip install -r requirements-analysis.txt
python tools/reproduce_figure3_data.py --output figure3_recalculated.csv
python tools/reproduce_primary_statistics.py
```

The first script regenerates Figure 3's 90 numerical rows and checks them
against the provided table. The second reconstructs 65 primary comparison
records and 15 upper-bound comparisons from the 180 seed-level records and
checks the paired means, intervals, exact Wilcoxon statistics and Holm results.
Both operate only on numerical CSVs and require no GPU or medical images.

## Data origin and interpretation

The study used a total of 13,909 image–mask pairs from PH² (200), ISIC 2018
Task 1 (3,694) and HAM10000 (10,015); the eligible analysis cohort contains
13,279 images. The segmentation comparisons are based on repeated, group-aware
splits with overlapping test-image appearances across seeds. The 180-row table
contains **seed-level group summaries**, not 180 independent patient samples.

- PH²: https://www.fc.up.pt/addi/ph2%20database.html
- ISIC 2018 Task 1: https://challenge.isic-archive.com/data/#2018
- HAM10000 and source lesion masks:
  https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/DBW86T

Users must obtain original images and masks directly from the relevant source
custodians under their respective terms. This archive does not redistribute
any third-party image or mask, nor does it claim a validated image-to-JPEG
conversion procedure.

## Licensing

- Python code and associated code documentation: MIT License (`LICENSE`).
- Four original numerical CSV datasets: CC BY 4.0 (`LICENSE-DATA.md`).
- Third-party dermoscopic datasets: **not included or relicensed**.

For citation, use `CITATION.cff` and, when a DOI is available, cite the
published article and archived package. No DOI has been assigned in this
prepublication candidate.
