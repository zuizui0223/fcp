# RGFCA Step 8C — 42,111-species tiered measurement allocation

Date: 2026-09-11 JST
Status: prospectively frozen before any new Step-8C image pixel is opened.

## Aim

Measure flower-colour information across the complete 42,111-species metadata-discovery universe without forcing an equal per-species photo depth that would discard most species. Preserve breadth and depth simultaneously by assigning an outcome-blind target depth from the recovered metadata-only capacity scan.

## Fixed species universe

All 42,111 species in `data/frozen/global_monte_carlo_species_discovery_combined_v2_species_v1.csv` remain members of the census. No species is removed because of current photo scarcity.

Step 8B established the following capacity counts after the frozen observer cap:

- >=1 photo: 41,322 species
- >=2: 33,944
- >=5: 24,612
- >=10: 18,301
- >=20: 12,985

The remaining 789 species are retained as structural current no-photo capacity.

## Fixed tier allocation

Assign each species exactly one target tier from `after_observer_cap` before opening any new Step-8C pixels:

- capacity >=20 -> target 20 photos;
- capacity 10–19 -> target 10;
- capacity 5–9 -> target 5;
- capacity 2–4 -> target 2;
- capacity 1 -> target 1;
- capacity 0 -> target 0 and retain as structural missingness.

Expected counts from Step 8B are:

- target 20: 12,985 species;
- target 10: 5,316 species;
- target 5: 6,311 species;
- target 2: 9,332 species;
- target 1: 7,378 species;
- target 0: 789 species.

Expected total target records before fresh-acquisition failures: **370,457 photographs**.

## Existing opened cohort

The previous discovery and reserve programs opened colour outcomes for 1,000 species total. Step 8C must identify overlap by species name/ID and mark these species `opened_existing`.

Opened species:

- remain members of the 42,111-species universe;
- may contribute to descriptive pooled atlas summaries only with explicit opened-cohort labeling;
- are not counted as fresh Step-8C validation species;
- are not re-labelled as untouched holdouts;
- need not be remeasured solely to satisfy the Step-8C target if their existing retained measurements are sufficient for the relevant descriptive estimator.

The fresh Step-8C acquisition queue is the target-positive complement after removing species already opened in the previous 1,000-species program.

## Fresh acquisition rule

The capacity scan is planning metadata, not the image-measurement candidate list. For fresh species, actual image selection must be a separate frozen draw using the existing RGFCA query/filter family:

- research-grade iNaturalist observation;
- species rank;
- flowering annotation;
- photo and unobscured coordinates;
- positional accuracy <=5 km;
- allowed photo licences inherited from the global acquisition contract;
- observer cap 2/species;
- deterministic selection independent of colour;
- no target relaxation after outcome;
- no replacement after measurement failure.

A species may return fewer records than its target in the fresh draw. Such shortfall is retained as an acquisition outcome and does not trigger target reduction, species replacement, extra pages chosen from colour outcomes, or post-outcome retries.

## Measurement

Use the existing fixed RGFCA flower-region/palette measurement implementation and retain its known validation limits. New results must not be represented as focal-species petal ground truth. Raw image bytes may be streamed and discarded after hashing/measurement; retain source IDs, provenance/rights fields, terminal status and derived measurement outputs.

## Estimands by depth

- target >=1: species-equal global colour composition / breadth atlas;
- target >=2/5: descriptive within-species heterogeneity, with explicit low-information uncertainty;
- target >=10: species-level polymorphism D screening and rarefaction diagnostics;
- target >=20: primary species-level D distribution and repeated within-species colour sampling;
- >=40 observed classifiable photos is still required for any previously frozen high-depth C*/S* or spatial-organization estimator; Step 8C does not lower that rule merely because its acquisition target is 20.

One or a few photos may never be used to infer biological absence of polymorphism.

## Primary reporting structure

Report separately:

1. the full 42,111-species census;
2. the current-photo-capacity profile;
3. the opened-existing 1,000-species cohort;
4. the fresh Step-8C measurable queue;
5. terminal fresh acquisition success/shortfall;
6. colour atlas summaries stratified by information depth.

Do not collapse these into one denominator.

## Claim boundary

This is an observation-program expansion. The 42,111 frame is not all angiosperms; photo availability is not a plant trait; the existing measurement model has unresolved focal-flower attribution/localization limitations; and low-depth colour observations do not establish presence/absence of intraspecific polymorphism.
