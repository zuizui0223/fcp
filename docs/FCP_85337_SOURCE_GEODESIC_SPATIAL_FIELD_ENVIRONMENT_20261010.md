# FCP original 85,337 repeated flower photographs: true-geodesic low-rank spatial field sensitivity (2026-10-10)

## Current evidence and the precise remaining issue

The original global 42,111 one-photo-per-species environmental atlas showed weak environment prediction after simple geography, but full thermal, rainfall, solar, wind, vapor, elevation and 10 SoilGrids covariates did **not** robustly improve over BOTH the geodesic 50/250km space kernel and the real dated LCVP shared-ancestral branch covariance on the limited 342/649 directly tipped species. Full parent-cohort phylogenetic correction remains HOLD.

The independent original **85,337 species×region-cell photos** (39,075 four-colour classifiable) showed a small positive environmental predictive increment after a train-only source species colour baseline and broad global spherical degree-3 geographic coordinates (source ecological CI verified). That nonlinear spatial basis, however, is a GLOBAL trend, not a source-photograph-distance-dependent spatial covariance. This sensitivity is intended to test whether the within-species environmental increment persists after a more FLEXIBLE spatial distance field.

## Source-fixed comparison without adapting photo-colour outcomes

Reopen only checked original taxon+observation+photo IDs from the archived 85,337 taxon×cell colour records and 100,543 expanded WorldClim BIO1–BIO19/solar/wind/vapor/soil original photo-site records. All 42,111 original photo taxa, 39,075 classified source taxon-cell images and 46,262 unclassifiable images remain in the original denominator. No photo is reclassified or substituted, and no iNaturalist/rasters are redownloaded.

Test separate, fixed COMPLETE-CASE populations:

- Climate without soil-selection: original classifiable image sites with elevation, all temperature/precipitation BIO, solar, wind and vapor fields complete; the prior source site-level study had 38,968 such source images, 20,546 photos with an identifiable heldout same-species training mean.
- Climate+soil: the same source fields plus original modelled SoilGrids pH, carbon, nitrogen, clay, available-water, CEC, sand, silt, density and coarse fragments; previous source study had 27,003 classifiable complete original taxon-cell photos and 11,136 originally photo-eligible for heldout same-species training mean.

Use original fivefold **GroupKFold by original equal-area photographed source region cell**. All models score EXACTLY the same held-out original photo IDs in each fixed environment-complete cohort. A species' flower-colour baseline is estimated ONLY using its other photographed regions in the training fold. Species without source-photo replication in other train cells are left outside this within-species estimand (but remain in the original 85,337 ledger). The phylogenetic component that is invariant within one nominal species is fully absorbed in its fixed species intercept; this analysis does NOT estimate a separately identified tree effect.

### Real geographic distance spatial covariance approximation

Within EACH training fold only, use the true latitude/longitude of TRAIN photographed source sites, transformed into a 3D geographic unit sphere. Select **96 fixed unsupervised MiniBatchKMeans source-photo geographic inducing knots** (random state 20261010, KMeans n_init=3). Heldout PHOTO COLOURS AND LOCATIONS ARE NOT USED TO FIT KNOTS; test coordinates are only used to evaluate distance to the fixed training knots.

Construct a finite-rank exponential RBF feature field `exp(-great_circle_photo_to_knot_km / bandwidth)` at original sites. Mathematically the Gram matrix `Phi Phi^T` defines a positive-semidefinite low-rank spatial covariance approximation, rather than solely polynomial latitude and longitude. Two spatial decay assumptions are reported, fixed **250km and 1000km**, never choose the most positive after seeing photo outcomes. Retain the source degree-3 unit-sphere polynomial geographical trend alongside those basis functions.

This is a valid DISTANCE-DEPENDENT SPATIAL BASIS, but is NOT a full GP posterior, exact exponential geodesic covariance process or a proof all residual spatial auto-correlation has been removed.

### Model families

For EACH same cohort/photo ID/fold:
- **M0 species + global nonlinear geographic degree-3**;
- **M1 = M0 + true-original-photo distance-based train-only low-rank spatial field**;
- **M2 = M1 + all temperature, precipitation, solar, vapor, wind, elevation and (if soil-complete) SoilGrids blocks**;
- **M3 leave one block out from M2**, without posthoc feature selection.

Original species-training-colour centred linear probability ridge is used consistently. Compare 4-category Brier, saved heldout error reductions, and 199 SOURCE SPECIES-cluster and original REGIONAL CELL-cluster fixed-fit loss bootstraps, with the same precise test photo IDs in every family. Clusters and genotype data are not independent confirmations; multiple abiotic contrasts are exploratory.

The diagnostic asks only: **does the previously reported original within-species environmental Brier improvement persist after a spatially local distance covariance basis is introduced?** A positive result would be a more credible *conditional photographic prediction* than the polynomial-only study, NOT true adaptation; a negative result would highlight spatial confounding/overfitting, not absence of biological climate effects.

## Nonnegotiable boundaries

- The historical source 42,111 broad species representation and original LCVP subset phylogenetic negative/uncertain results cannot be merged or averaged with this within-species estimand.
- Fitting 96 source site knots without photo colour labels is an exploratory chosen low-rank spatial prior, not complete spatial random-field inference. Next meaningful stricter test would include residual Moran/variogram and truly independently measured source colour and plant breeding population.
- Original 50–100km microgeographic photo-label shuffles of the congeneric source showed no robust independent rainfall-adaptation evidence. Do not reopen that strong claim on any single favourable 250/1000km modelling condition.
- This branch is a draft nested PR; original New Phytologist 1,499 taxa, main and independently held-out future 2,000+730 photo taxa remain unchanged.
