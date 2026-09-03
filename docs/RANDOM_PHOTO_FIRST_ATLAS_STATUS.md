# Random photo-first flower-colour boundary atlas — current state

## Purpose

This is a **new prospective experiment**, separate from both the frozen six-species Chapter-1 analysis and PR #21's terminal 200-species / 60,000-photo experiment.

The question is:

> When georeferenced flowering photographs are repeatedly resampled without fixing a focal species set, do the same flower-colour transition edges recur more strongly than expected under a species-conditioned null?

The intended sequence is:

`fresh metadata-only candidate pool -> location-blind photo measurement -> repeated photo-first H1 persistence -> H2 environmental concordance`

No step may rescue, subset, relabel, or reopen the terminal PR #21 records.

## Current inferential contract

Primary H1 is frozen before the fresh candidate pool is opened:

- equal-area longitude × sin(latitude) grid: **18 × 9 = 162 cells**;
- **10,000 photographs per replicate**;
- **200 replicate samples**;
- cell-first round-robin sampling;
- maximum **2 photographs per species per cell per replicate**;
- minimum **5 classifiable photographs per cell** for edge evaluation;
- biological coarse-colour states: `white`, `yellow_orange`, `red_pink`, `blue_purple`;
- `mixed_uncertain` is **structural measurement missingness**, never a fifth colour state;
- edge metric: Jensen-Shannon divergence of the four-state cell composition;
- transition edges: exact top 10% of evaluable edge intensities per replicate, with random tie-breaking;
- persistence denominator: only replicates in which both adjacent cells are evaluable;
- global statistic: opportunity-weighted variance of edge persistence around the realized transition rate;
- primary null: 999 permutations of **classifiable morph labels within species**, preserving locations, species frequencies, sampling seeds, and the geographic mask of `mixed_uncertain` measurements.

The active implementation is `fcp_pipeline/photo_first_atlas_v2.py`.

## Fresh candidate-pool freeze

The metadata pool is separately frozen before image pixels:

- source: iNaturalist API v1;
- root taxon: Angiospermae (`taxon_id=47125`);
- Research Grade;
- flowering annotation (`term_id=12`, `term_value_id=13`);
- species-rank observations with photos and coordinates;
- positional accuracy ≤5 km;
- allowed Creative Commons photo licences only;
- one `order_by=random` request per equal-area cell;
- one page per cell, maximum 200 returned candidates;
- maximum raw pool size: 32,400 observations;
- no fixed or targeted species list;
- returned coordinates are re-gridded locally and bbox-boundary duplicates are removed globally;
- exact observation/photo IDs and the candidate-table SHA-256 are frozen before any candidate image pixels are opened.

The freeze is one-shot. If any durable candidate-pool output already exists, the acquisition command refuses to issue a new random query. A rerun can therefore never become a favourable replacement replicate.

## Premeasurement capacity gate

Before candidate image pixels may open, the frozen metadata pool must support the fixed H1 sample size under the species cap.

The gate computes:

`sum over cell × species of min(candidate_count, 2)`

If this capacity is <10,000, the new atlas stops as:

`not_evaluable_candidate_sampling_capacity_before_pixels`

and candidate image pixels remain unopened.

If capacity is ≥10,000, the separately frozen photo-measurement stage may proceed. This capacity gate is about the sampling frame only; it is not a biological result.

## Frozen photo-measurement path

The measurement stage is now frozen before any fresh candidate pixel. It deliberately reuses only two pieces of validated infrastructure, not PR #21's terminal records or species/cohort completeness rules:

1. the independently qualified **ROI-v4** flower detector/mask estimator from PR #21;
2. the generic species-independent **sRGB reference palette** frozen in the completed Chapter-1 Florence pipeline.

The measurement path is:

`candidate image -> ROI-v4 flower mask -> masked RGB pixels -> fixed CIELAB nearest-palette assignment -> four coarse biological colour groups`

The four biological groups are the same H1 states: `white`, `yellow_orange`, `red_pink`, and `blue_purple`. Green, brown and black palette anchors are nuisance anchors and are removed from the biological denominator. The palette is not re-fitted to the fresh images.

The exact mapping code is `fcp_pipeline/photo_first_measurement.py`; its frozen contract is `docs/supporting/random_photo_first_measurement_contract_v1.json`.

The worker is location-blind. Species identity, taxon ID, observation/photo IDs, coordinates, observer, date, climate and pollinator context are unavailable during measurement. All frozen candidate rows must finish in exactly one terminal photo state. Acquisition/ROI/palette failures are retained as `mixed_uncertain`; no image is replaced.

Unlike PR #21, there is **no species-level or cohort-level measurement-completeness gate**. Photo-level non-evaluability remains in the fixed colour-blind sampling frame and is handled later by the edge-specific opportunity denominator.

## Quality-safety correction made before fresh data

The first software prototype allowed `mixed_uncertain` to enter the cell composition as a fifth state. That was rejected before fresh candidate data were opened because geographically structured image/ROI failure could masquerade as a biological colour transition.

The active v2 implementation therefore:

1. samples photographs before looking at colour;
2. retains `mixed_uncertain` rows in the fixed sample;
3. excludes them from biological morph composition and the minimum classifiable-photo count;
4. holds their exact geographic positions fixed under the null;
5. shuffles only classifiable biological morphs within species.

This makes measurement missingness visible without converting it into flower-colour signal.

## Execution state

Implemented and frozen before fresh pixels:

- H1 analysis core;
- fixed-size fail-closed sampling;
- exact top-transition fraction with tie handling;
- opportunity-denominator persistence;
- species-conditioned null;
- quality-safe structural-missing treatment;
- synthetic planted-boundary validation runner;
- metadata-only iNaturalist candidate-pool freezer;
- one-shot query guard;
- premeasurement 10,000-photo capacity gate;
- ROI-v4 + fixed-palette measurement core and tests;
- photo-measurement contract with location blinding and no species/cohort completeness gate;
- CI workflow that validates H1 first and then performs the one-shot metadata freeze only after software validation passes.

Not yet a biological result:

- the fresh metadata pool has not yet been durably frozen in the repository;
- fresh candidate pixels have not been opened;
- the new measurement worker has not run on the fresh pool;
- no new flower-colour persistence statistic or p-value exists;
- H2 climate concordance is not open.

## Legacy evidence retained unchanged

The new experiment does not alter these completed results:

- six-species Stage A: within-species continuous colour is spatially organized (`p = 0.0113`);
- six-species Stage B: one universal shared global transition boundary was not confirmed at the primary scale (`p = 0.0906`);
- PR #21 terminal scale-out: 60,000 location-blind measurements completed, but the experiment was `not_evaluable` under its prospectively frozen species/cohort measurement-completeness gate.

Those results motivate the new photo-first persistence question but are not input data for it.
