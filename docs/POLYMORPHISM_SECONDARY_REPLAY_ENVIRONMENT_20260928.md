# Secondary-analysis replay environment — 2026-09-28

## BIO5 / BIO14 / SRAD analyses

The original successful 2026-09-25 GitHub Actions runs used:

- Python 3.11.16
- NumPy 2.4.6
- pandas 3.0.6
- SciPy 1.17.1
- statsmodels 0.15.0
- rasterio 1.4.4

The complete pip environment recovered from the successful job log is frozen in:

`requirements-np-bio5-replay-20260928.txt`

Original successful runs:

- third-cohort white-environment mechanism: `36128744659`
- legacy BIO5 transport: `36128744674`

Original result artifacts:

- `10860659652` — `white-environment-mechanism-v1`
- `10860821810` — `legacy-white-bio5-replication-v1`

## H3 analyses

The original H3 result JSONs record their own runtime:

H3a:
- R 4.6.1
- ape 5.8.1
- phytools 2.5.2
- jsonlite 2.0.0

H3b:
- R 4.6.1
- ape 5.8.1
- phylolm 2.6.5
- jsonlite 2.0.0

The replay workflow installs these versions explicitly before rerunning the archived scripts against the permanent tree/covariate inputs.

## Role

This file freezes the software environment only. It does not alter any statistical definition, threshold, result, or claim.
