# RGFCA Step 8 — expanded photo-first atlas

Date: 2026-09-10 JST
Status: prospectively frozen before any new Step-8 image pixels are opened.

## Aim

Replace the current 1,000-species measured universe with a substantially larger photo-first atlas and estimate atlas-feature saturation against a much larger outcome-blind species universe. The previous 350-species threshold is treated only as an internal rarefaction result within the 363–369 eligible-species tranches, not as a general biological minimum.

## Outcome-blind source frame

Use only the recovered metadata-only capacity scan v3:

- `docs/supporting/global_monte_carlo_capacity_scan_manifest_v3.json`
- `data/frozen/global_monte_carlo_capacity_scan_selected_species_v3.csv`
- required status: `complete_metadata_only_capacity_scan_target_selected_after_transport_recovery_v3`
- selected raw-photo target: exactly 100
- capacity-selected species: exactly 4,730
- flower colour and image pixels were not used to define this frame.

Do not use `data/global_flower_colour_species_ranked.csv` or any literature-derived polymorphism ranking for Step 8 species selection.

## Existing opened cohort

The already measured discovery and reserve cohorts together contain 1,000 species × 100 terminal photo records. Their colour outcomes are opened and they are not an untouched holdout. Step 8 must track overlap explicitly and distinguish:

1. existing opened species;
2. new capacity-frame species whose Step-8 pixels have not previously been opened.

No opened species may be relabelled as an independent validation species.

## Expanded nested atlas sizes

Construct one deterministic, outcome-blind ordering of the 4,730 capacity-selected species using SHA256(`20260910|inat_taxon_id|species`) and define exact nested prefixes:

- N = 500
- N = 1,000
- N = 2,000
- N = 3,000
- N = 4,000
- N = 4,730

These are analysis checkpoints, not separately optimized species sets. A species cannot be swapped after any colour result is opened.

## Photo acquisition and measurement

For every newly measured Step-8 species retain the existing candidate-acquisition rules unless a technical incompatibility is documented before pixels:

- research-grade, flowering observations with photos and coordinates;
- positional accuracy <=5 km;
- allowed photo licences as in the existing global candidate-acquisition contract;
- observer cap 2/species;
- deterministic geographic maximin selection;
- exactly 100 raw photo records/species;
- no replacement after measurement failure;
- fixed existing flower ROI/palette measurement code and fixed background annulus where available.

Raw image bytes should be streamed for measurement and need not be retained in Git. Persist metadata, source identifiers, hashes/receipts where available, terminal status, and derived measurement outputs. Do not create a permanent local 473,000-image archive solely for Step 8.

## Primary atlas features

Estimate saturation separately rather than requiring one conjunctive threshold to define a single minimum species count.

A. Global four-morph composition.
B. Distribution of species-level polymorphism D.
C. Fraction of eligible polymorphic species in photo-derived organization states: C* only, S* only, C*+S*, unresolved.
D. Distribution of within-species geographic colour-organization statistics.
E. Equal-area geographic coverage / occupied atlas cells.

Shared global-boundary recurrence is retained only as a previously tested non-supported estimand; Step 8 does not retune that test.

## Saturation analysis

For each nested N and for within-species raw-photo depths 20, 40, 60, 80, 100, generate 200 deterministic balanced rarefaction realizations. Compare each feature with the largest currently available Step-8 atlas using feature-appropriate distances:

- composition: cosine similarity;
- D distribution: Wasserstein distance and quantile error;
- C*/S* composition: absolute proportion error;
- spatial-organization distribution: Wasserstein distance / rank correlation where defined;
- geography: occupied-cell retention and cell-wise species-count correlation.

Report full saturation curves. Do not collapse them into a single 'minimum' unless a predeclared feature-specific criterion is satisfied.

## Feature-specific operational sufficiency criteria

These are descriptive design thresholds, not biological truths:

- global colour composition: cosine >=0.95 in >=90% of realizations;
- D distribution: Wasserstein <=0.03 and median/upper-quartile absolute error <=0.03 in >=90%;
- C*/S* state composition: each estimable state proportion error <=0.05 in >=90%;
- geographic coverage: occupied-cell retention >=0.90 and cell species-count correlation >=0.90 in >=90%;
- spatial-organization distribution: Wasserstein <=0.05 in >=90% where the statistic is estimable.

The paper must report the saturation point separately for each feature.

## Main interpretation boundary

Step 8 asks how much globally distributed image evidence is needed to recover stable flower-colour atlas summaries. It does not establish a census of all angiosperms, causal selection, focal-species petal truth, phylogenetic signal, or a universal shared flower-colour boundary.

The old 34-species C/S mechanism transfer remains closed after Steps 7C–7F and is not reopened by Step 8.
