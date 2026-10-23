# FCP: are photographed flower colours multimodal, clinal, or locally complex?

**Date:** 23 October 2026. **Analysis status:** implemented exploratory, not a registered test of inherited morphs. The existing successful nine-colour P2 geographic result remains immutable in the upstream stacked analysis PR #260.

## Core ecological question and separation of estimands

Original photographs from one source (447,500/4,475 taxa; 227,107 four-colour classifiable) reveal local continuous nine-palette hue resemblance, but that does not answer whether the distributions consist of multiple discrete phenotypic modes, continuous geographical clines, or spatial mosaics. No response should assume exactly 2/3/4 genetic colour morphs on the basis of coarse label count. This new analysis treats **two independent dimensions**:

1. **Photographic phenotype geometry:** is the dominant unsupervised Hellinger nine-hue coordinate's *density* convincingly multimodal, against both a single normal and a flexible unimodal skew-normal competitor?
2. **Geographic predictive geometry:** does a broad geographic linear trend predict held-out colour vectors across separate geographic blocks, and/or do locally smoothed residuals predict colours held out by observer?

The second can provide *broad-trend* and *local-nonlinearity* **signatures**. It cannot prove literal monotone clines, discrete patches, mosaics, separate breeding populations or modes' inheritance.

## Original photo provenance and analysis frame

The reader authenticates the immutable 447,500-photo source with original SHA256 8e7f10f0427927208485b2da6fa9c46bb681d6f3611e53b8dfdc606cfc8dabe6; the release source archive is 99f7ee68c3fd009e8ab60fd0cfcdf1fb407e3a8de07edd38341b27889310bc25. Original species with fewer than 40 classifiable photographs are ineligible for the main geometry audit; at least 30 original geovalid pairs within 50 km are needed. The exact strict source baseline of 1,761 eligible taxa and +0.0176149021685411 geographic four-state depletion is an immutable replay gate. Species not passing the gate remain unavailable for this *subset*, not monomorphic.

## Photo mode discrimination (no forced K=2)

For each eligible species, use normalized 9-hue vectors, apply Hellinger sqrt, and calculate the first principal projection without location, observation date or identity. Fit one-dimensional Gaussian mixture densities K=1..4 (BIC), compare the chosen fit with one-Gaussian and unimodal skew-normal models, then examine whether the fitted **density itself** has at least two peaks, each assigned mixture component has at least 5 source photos and 3 distinct observers, the smallest mixture weight >=0.10, and BIC improvement >=10 over both unimodal competitors.

A potential multimodal taxon is then tested on 12 colour-blind 80%-observer subsamples, requiring at least 70% of subsamples to repeat the evidence. Return **recurrent_photographic_mode_candidate**, **unstable_photo_mode_candidate**, or **not_resolved_as_multimodal**, not genetic morph count. Only PC1 is tested; additional modes on PC2/PC3 are not excluded by negative PC1 results. Strong camera, ROI, optical or stage artefacts can themselves create photographic modes.

## Spatial shape (month baseline vs nonlinear model)

The target is the full nine-hue Hellinger photographic vector, not fixed four-category membership. Compare (i) regularized photographed-calendar-month means, (ii) calendar-month + a linear trend in 3D Earth-centred geographic coordinates, and (iii) the linear model plus a Gaussian-kernel local smoother of residual colour. The kernel bandwidth is max(30 km, 0.18 times the training geographic span); it is frozen, with no outcome-tuned parameter search.

Use the same models in two *geographically different* validation conditions: 3 observer-disjoint folds (interpolation into already sampled geographical neighbourhoods) and 3 coordinate-only K-means geographic-block folds (transfer to new geographic blocks). Each test group >=8 photos and train group >=30, and geographic-block transfer requires >=100 km sampled span. Models are compared using holdout summed squared errors relative to the month-only baseline.

**Reporting categories are nonexclusive predictive signatures, and thresholds are exploratory operational rules**: transferable broad linear trend if gain>0.08 in both split regimes; local residual structure if nonlinear increment over the linear model>0.08 in observer holdouts; combinations or unresolved otherwise. Geographic nonlinearity may be a curved cline, real patch mosaic or camera artifact. A per-species category is **not** a detected biological mosaic.

## Specific tests and deliberately restrained conclusions

- A 2-component GMM fit alone is NOT evidence for two modes; the fitted density may be unimodal, skew or unreliable. Likewise, failure to reject a flexible unimodal density cannot establish a true continuous biological phenotype.
- Spatially local prediction is not independent of photo-source concentration; the observer split lessens direct repeat-photo leakage, geographic-block splitting quantifies how much extrapolates to new regions.
- Observed pseudo-categories cannot establish their own prevalence outside the classified sample. Original source classification failures comprise 220,393 photos and remain unfilled.
- The prior strict <=5-km **same-year/month, distinct-observer photographic null** was unsupported (P=0.300/0.275). The new model is not a substitute for that negative test.
- Real genetic FCP, phenotypic plasticity, photographic colour fidelity, pigment pathways, fitness, population structure and gene flow are *not identified* here.
- If few taxa have adequate separated geographic blocks, report that limitation instead of forcing a cline/mosaic split.

## Reproduce

Code: scripts/analysis/fcp_photo_modes_geographic_shape_20261023.py
Workflow: .github/workflows/fcp-photo-modes-geo-shape-20261023.yml
Output: summary.json and per_species_mode_geo_metrics.csv, without exact photo IDs, observers or coordinates. The workflow first passes a synthetic unimodal/multimodal and linear/nonlinear calibration, authenticates both frozen source checksums, reproduces the old strict P2 and only then permits an empirical claim.

**Planned next decision after the empirical run:** determine how many taxa yield *robust observed modes*, *transferable broad trends* and *local nonlinear signatures*, and whether these evidence dimensions overlap. Interpret overlap as a distribution of photographic phenotypes across geography, not a count of natural morphs.
