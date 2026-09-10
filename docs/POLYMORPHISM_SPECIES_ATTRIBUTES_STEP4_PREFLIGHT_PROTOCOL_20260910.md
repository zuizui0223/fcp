# Polymorphism species attributes — Step 4 preflight protocol

Date frozen: 2026-09-10 JST
Branch: `analysis/polymorphism-species-attributes-step4b`
Parent: `analysis/polymorphism-spatial-reserve-step5b` at `b08232d7f1e8784a1c94c5213c8f525ecf65fae4`

## Purpose

Freeze the species-attribute covariate panel **before any new association with flower-colour polymorphism D is computed**. This step is source/coverage/taxonomy preflight only.

The biological question is: **which kinds of plant species carry more flower-colour polymorphism?**

This protocol does not alter the already-frozen D estimator and does not reopen the rejected M1/M2 directionality claim.

## Outcome lock

The later Step 4 association analysis must use the 369-species discovery D table already frozen at:

`results/polymorphism_directionality_step1_20260909/species_metrics.csv`

Primary response: frozen four-morph Gini-Simpson diversity `D = 1 - sum(p_k^2)`.

Predeclared finite-sample sensitivity: `D_unbiased = n/(n-1) * D` for `n > 1`.

**This preflight script may read `species` and `n` from the frozen table but must not read or compute any D/covariate association.**

## Block A covariates frozen now

### A1. Genus

Derive as the first token of the frozen binomial species name. No taxonomic backfilling after results are opened.

### A2. Family

Use the frozen iNaturalist taxon ID keyed to each FCP species, then retrieve the family-rank ancestor from the public iNaturalist v1 taxon endpoint. The preflight must save the retrieved family mapping and retrieval receipt before D is opened.

Primary taxon-ID source: `data/frozen/global_monte_carlo_capacity_scan_selected_species_v3.csv`.

If an FCP species has no retrievable family, family is missing; it must not be guessed from genus after D is opened.

### A3. Geographic span / opportunity

Primary proxy: `log1p(maximum_span_km_after_observer_cap)` from:

`data/frozen/global_monte_carlo_candidate_species_audit_v1.csv`

This is explicitly a **sampled geographic span/opportunity proxy**, not a true biological range-size estimate.

Predeclared sensitivity: `log1p(maximum_span_km)` from:

`data/frozen/global_monte_carlo_capacity_scan_selected_species_v3.csv`

### A4. Latitude centroid

For each FCP species, compute the mean latitude from **all colour-blind candidate-photo rows** in:

`data/frozen/global_monte_carlo_candidate_photos_v1.csv`

Primary predictor: `abs(mean_latitude)` (distance of the sampled centroid from the equator).

Descriptive/sensitivity field retained: signed `mean_latitude`.

No filtering by flower classification, morph, D, or image outcome is allowed when calculating this predictor.

### A5. Pollination guild

Reuse the already-frozen `zuizui0223/island` all-master trait campaign; do not run a new LLM trait inference.

Frozen source:

- island workflow run: `29409415292`
- artifact: `all-master-trait-ledger-29409415292`
- artifact ID: `8340477381`
- recorded artifact digest: `sha256:0cbba21f384d2b595272f13ed7be5a4b89d8ab0b77c6a0dd7a355fb3b3be1665`

Required files:

- `all_species_traits.csv`
- `all_species_trait_evidence.csv.gz`

A pollination value is admissible only when the final `pollination_guild` is non-`unknown` **and the same species/value pair has at least one `source_backed=True` row with `field == pollination_guild` in the evidence ledger**. This deliberately excludes final values that cannot be traced to a matching source-backed guild record.

The preflight records raw admissible guild values only. No D association is computed.

### A6. Measurement-quality control

Retain `n_classifiable` (`n` in the frozen D table) as a predeclared control/sensitivity variable. This is not treated as a biological species attribute.

## Pollination ontology for the later association test

This mapping is frozen before FCP coverage is inspected:

- `bees`, `bumblebees` -> `bee`
- `flies`, `birds`, `moths`, `butterflies` -> `other_animal`
- `wind` -> `wind`
- `mixed` -> `mixed`
- `self` -> `self`

The later omnibus pollination test may include only mapped categories with at least 10 FCP species. The deterministic coverage gate is:

1. at least 50 FCP species with admissible source-backed pollination values overall, and
2. at least two mapped categories with `n >= 10`.

If the gate fails, pollination is reported as underpowered/unavailable and **must not be rescued by post-outcome web or LLM backfilling**.

## Block B life form: source preflight before D

Life form remains a separate planned block because no tracked FCP source existed in the prior header audit.

Before any D association is opened, query the **published GIFT v3.2 metadata only** and save candidate trait metadata whose fields contain `life`, `growth`, `habit`, or `wood`. Also save the GIFT version receipt and public-reference availability counts. This preflight does not choose categories by looking at D.

A later life-form analysis is allowed only after an exact GIFT trait ID and categorical mapping are frozen from this metadata. If no adequate source/coverage exists, Block B closes as unavailable rather than being backfilled post hoc.

## Later statistical tests — frozen family

After the preflight covariate panel is committed, the Step 4 analysis may compute exactly these primary tests:

1. sampled geographic span: two-sided Spearman association of D with primary log-span proxy;
2. absolute latitude centroid: two-sided Spearman association;
3. family taxonomic clustering;
4. genus taxonomic clustering;
5. pollination-guild omnibus, only if the frozen coverage gate passes;
6. life-form omnibus, only if Block B is frozen and passes its prospective coverage gate.

For continuous predictors, report raw Spearman rho plus a sensitivity using `D_unbiased` and a partial-rank sensitivity controlling `n_classifiable`.

For family/genus, call the result **taxonomic clustering**, not phylogenetic signal. The exact clustering statistic/permutation implementation must be written into the association runner before it reads D.

Multiplicity: retain raw p-values for all six requested attribute families and apply a single Holm correction after Block B is either tested or formally closed. A closed/failed-gate family receives `p=1` for the final six-slot Holm family; it cannot disappear from multiplicity accounting.

## Hard stop rules

- No new image acquisition.
- No new colour measurement.
- No D association in preflight.
- No trait added, removed, recoded, or web/LLM-backfilled because of a Step 4 outcome.
- `family`, `pollination`, and life-form provenance must be saved before D is opened.
- Existing positive M3 mean-position result remains an auxiliary colour-space descriptor and is **not** part of this requested species-attribute test family.
