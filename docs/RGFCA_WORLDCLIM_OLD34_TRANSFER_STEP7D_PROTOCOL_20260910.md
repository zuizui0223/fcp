# RGFCA Step 7D — historical WorldClim C/S hypothesis transfer

Date: 2026-09-10 JST
Status: prospectively frozen before opening Step 7D environmental outcomes.

## Aim

Transfer the historical 34-species WorldClim-based C/S hypothesis family to the photo-derived RGFCA organization system while preserving the historical environmental estimands as closely as possible.

This is a **transferred test**, not a confirmatory replication: historical C/S are literature-supported biological organization states, whereas RGFCA C*/S* are photo-derived descriptors.

## Source freeze

- WorldClim: **WorldClim 2.1, 10 arc-minute BIO rasters**.
- Exact source URL: `https://geodata.ucdavis.edu/climate/worldclim/2_1/base/wc2.1_10m_bio.zip`.
- Retained BIO variables, matching the historical pipeline: **BIO1, BIO4, BIO5, BIO6, BIO7, BIO12, BIO14, BIO15, BIO17**.
- Historical implementation reference: `scripts/data/compute_climatic_niche_metrics.py` and `scripts/run_jbi_space_time_hypothesis_screen.py` on the preserved Chapter-1 line.

No alternative climate product, resolution, or BIO-variable set may be substituted after outcomes are opened.

## RGFCA coordinates used for predictors

Use the frozen measured-photo tables:

- discovery: `data/derived/global_monte_carlo_measured_photos_v1.csv`;
- reserve: `data/derived/rgfca_reserve_replication_measured_photos_v1.csv`.

For species admitted to Step 7A, environmental predictor construction uses **all frozen raw candidate photo rows with valid coordinates**, irrespective of `global_classifiable` or morph. This keeps climate construction colour-blind and avoids conditioning environmental predictors on colour-classification success.

Within each species, duplicate identical nine-variable WorldClim vectors are removed exactly as in the historical climate pipeline. Species require >=20 unique occupied climate vectors for H1/H3/H0 metrics.

## PCA and tranche independence

Discovery and reserve climate spaces are built **independently**. Within each tranche, standardize the nine BIO variables globally over all Step-7A species climate cells and fit PCA with three components. Do not pool the two tranches to estimate PCA axes.

This preserves independent replication at the cost of allowing arbitrary PC sign/orientation differences; all H3/H0 metrics are rotation/sign invariant in the retained 3D space.

## Primary response contrast

The primary historical analogue is restricted to pure photo states:

- `local_cooccurrence_only` = C* only;
- `spatial_segregation_only` = S* only.

Exclude `cooccurrence_and_segregation` and `unresolved_photo_state` from the primary C*/S* contrast. This is fixed because the historical 34-species flagship comparison was a binary coexistence-versus-geographic-structure contrast.

No pooling across discovery and reserve is allowed to rescue significance.

## Historical metrics retained

### H1a — temperature seasonality

Metric: species mean occupied **BIO4**.
Prediction: higher BIO4 favours C* relative to S*.
Test statistic: median metric(S* only) - median metric(C* only), predicted **negative**.

### H1b — precipitation seasonality

Metric: species mean occupied **BIO15**.
Prediction: higher BIO15 favours C* relative to S*.
Test statistic: S* - C* median difference, predicted **negative**.

### H3a — spatial environmental turnover

Historical statistic retained:

1. calculate all within-species geographic and PC1-PC3 environmental pair distances (sample 25,000 pairs with fixed seed if total pairs >20,000);
2. identify lower and upper quartiles of geographic distance;
3. `spatial_niche_turnover = (median environmental distance in upper geographic quartile - median environmental distance in lower geographic quartile) / median environmental distance`.

Prediction: larger turnover in S*. S* - C* median difference predicted **positive**.

### H3b — regional niche displacement

For species with >=40 climate cells, split coordinates into two geographic sectors with spherical-coordinate KMeans (`n_clusters=2`, `n_init=20`, fixed seed), requiring >=15 cells per sector.

`regional_centroid_sep = distance between PC1-PC3 sector centroids / within-species RMS PC distance`.

Prediction: larger separation in S*. S* - C* median difference predicted **positive**.

### H3c — regional niche overlap

Use the same two geographic sectors and the historical regularized Gaussian Bhattacharyya affinity in PC1-PC3.

Prediction: lower overlap in S*. S* - C* median difference predicted **negative**.

### H0 — total niche-size negative control

Rarefy each eligible species to 20 occupied PC1-PC3 climate cells for **99 fixed-seed draws** and take median 3D convex-hull volume.

H0 has no directional biological prediction. Report the two-sided C*/S* permutation contrast as a negative-control diagnostic. A non-significant H0 is not interpreted as proof of equality.

## Testing and multiplicity

For H1a, H1b, H3a, H3b, H3c:

- statistic: S*-only median minus C*-only median;
- 20,000 label permutations preserving the observed pure-state counts;
- one-sided p in the predeclared direction;
- Holm correction across the five directional hypotheses within each tranche.

A historical environmental hypothesis is `two_tranche_recurrent` only if:

1. effect direction matches prediction in discovery and reserve; and
2. Holm-adjusted p <=0.05 independently in both.

H0 is reported separately with a two-sided 20,000-permutation p and does not enter the five-hypothesis Holm family.

## H8 boundary

H8 (pigment-mediated abiotic trade-offs) is **not** represented by an unspecified combination of the nine WorldClim variables. It remains unopened until drought/UV/soil/extreme-stress predictors are prospectively defined and sourced.

## Interpretation boundary

- BIO4/BIO15 are climatological seasonality, not year-specific temporal climate.
- Photo coordinates represent sampled occupied environments, not the complete fundamental niche.
- C*/S* remain photo-derived states.
- H3 associations do not identify local adaptation, selection, or causal environmental sorting.
- Failure to replicate an old34 hypothesis under RGFCA C*/S* may reflect a genuine scale/response difference and must not trigger threshold or state redefinition.