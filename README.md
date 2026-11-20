# MitoRedox Atlas analysis code

This repository contains the analysis scripts accompanying *A multi-tissue atlas of mitochondrial cysteine redox modifications during ageing*. This code snapshot corresponds to the 2026-11-20 manuscript review release.

| Step | Analysis |
| --- | --- |
| 01 | Prepare mitochondrial protein and modification measurements |
| 02 | Fit tissue and ageing contrasts |
| 03 | Test pathway enrichment and reproduce functional network partitions |
| 04 | Reproduce protein–PTM and shared-site coupling analyses |
| 05 | Fit sequence, structural and external SSG models |
| 06 | Reproduce protein sequence alignments and variant evidence summaries |

## Reproduce the main analyses

Use R 4.3.2 with data.table 1.17.0, limma 3.58.1, sandwich 3.1-1 and igraph 2.2.1. The variant step uses Python 3.12 and its standard library.

Obtain the matching review archive `mitoptm-review-2026-11-20.zip` through the manuscript submission materials. Place its `inputs/` and `reference/` directories next to `run_all.R`. These directories contain the frozen feature-level measurements, model inputs and expected results used for verification.

From the repository root, run:

```sh
Rscript run_all.R
python -B analysis/06_map_variants.py
```

The scripts run offline, write tables under `results/` and compare the reproduced results against the supplied reference tables. The workflow starts with feature-level mass-spectrometry quantification and uses frozen structure-derived predictors and external variant evidence. It re-runs the main numerical analyses and sequence alignments. Upstream spectra searches, structural feature extraction, database acquisition and final figure assembly are documented in the matching review archive.
