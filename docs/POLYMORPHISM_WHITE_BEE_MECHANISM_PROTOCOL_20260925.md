# FCP white–bee mechanism protocol — 2026-09-25

Status: **prospectively frozen after an outcome-blind bee-coverage gate and before opening flower-colour outcomes against bee metrics**.

This is a separate mechanism-analysis line. It does not alter the frozen New Phytologist manuscript, H2 estimand, H2 verdict, or submission package.

## Biological question

Does recurrent white-versus-nonwhite flower-colour variation occur preferentially in plant species with weaker or narrower bee associations?

The specific hypothesis tested here is **bee-mediated chromatic selection**: if diverse bee assemblages favor the maintenance of chromatic floral signals, species with broader bee partner sets should show a lower white fraction.

This is not a test of all pollinators. Moths, butterflies, flies, birds and bats are outside this test.

## Independent bee source

Use Noori et al. 2026 curated GloBI bee–plant interactions:
- Zenodo record 18303036, version 3.1;
- file `GloBI_Curated.csv`;
- harmonized plant field `plant_species`;
- bee field `bee_species`;
- plant family field `plant_family`;
- source database field `database`.

The outcome-blind coverage gate found exact plant-species coverage for 346/499 frozen third-cohort species, exceeding the prespecified minimum of 100.

## Bee metric frozen before colour opening

For every plant species in the curated GloBI table calculate:
- `n_bee_records`: number of interaction rows;
- `n_bee_species`: number of distinct bee species;
- `n_bee_databases`: number of distinct contributing database labels.

Restrict metric calibration to plants with >=5 interaction records.

Fit across the complete curated plant set:

`
log1p(n_bee_species) ~ log1p(n_bee_records) + log1p(n_bee_databases)
`

The standardized residual is the primary predictor:

`
bee_partner_breadth_resid_z
`

Positive values mean more distinct bee partners than expected from interaction-record and source-database effort.

No alternative effort model is selected after colour outcomes are opened.

## Flower-colour outcome

Use the frozen third-cohort measurement artifact and frozen response-blind highlight control.

Eligible photo rows:
- pass `global_classifiable`;
- belong to white, yellow/orange, red/pink, or blue/purple;
- have an available near-clip metric;
- are not in the frozen high-clip set.

Species enter the primary analysis only if they retain:
- >=40 eligible classifiable rows;
- >=5 white rows;
- >=5 non-white rows;
- an exact GloBI plant match;
- >=5 GloBI bee interaction records.

For each species:

`
white_fraction = n_white / (n_white + n_nonwhite)
`

Species are equally weighted in the primary comparative analysis.

## Primary prediction

`
bee_partner_breadth_resid_z ↑  -> white_fraction ↓
`

Primary test:
- Spearman correlation across species;
- one-sided permutation p-value for a negative association;
- 9,999 permutations with seed 20260925.

Minimum inferential support:
- >=100 eligible species.

## Corroborative model

Fit an equal-species-weight fractional binomial GLM:

`
white_fraction ~ bee_partner_breadth_resid_z
               + z(log1p(n_bee_records))
               + z(log1p(n_bee_databases))
`

Use plant-family clustered sandwich uncertainty.

The bee-breadth coefficient must be negative with cluster-robust p < 0.05 to corroborate the primary result.

## Mechanism gate

The bee-breadth mechanism is labeled `gate_pass` only when:
1. >=100 eligible species;
2. Spearman rho < 0;
3. one-sided permutation p < 0.05;
4. fractional-GLM bee-breadth coefficient < 0;
5. family-clustered p < 0.05.

Otherwise the result is `not_supported_under_this_test` or `not_estimable`.

## Robustness outputs

Report but do not use to rescue the gate:
- raw bee partner richness versus white fraction;
- correlation of the primary bee metric with GloBI record effort;
- leave-one-family-out coefficient range;
- number of represented plant families.

## Hard nonclaims

A positive result does not establish:
- pigmented -> white evolutionary transition direction;
- bee preference for individual white versus coloured morphs;
- pollination effectiveness rather than flower visitation/interaction documentation;
- causality;
- a mechanism involving moths, flies, butterflies, birds or bats.

A null result does not reject pollinator-mediated evolution generally; it rejects only this independent cross-species bee-breadth prediction under the frozen design.
