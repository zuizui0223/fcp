# White-state pollinator screen — 2026-09-25

Status: **exploratory post-outcome mechanism screen; not part of the frozen New Phytologist manuscript and not confirmatory evidence.**

## Question

Do recurrent measured white flower states share a broad pollinator context across species?

This screen is deliberately separated from the prospective WorldClim mechanism family. The flower-colour outcome was already known before this pollinator join. The pollinator source therefore cannot upgrade the evidence to a prospective causal result.

## Outcome-blind multi-guild coverage audit

A GloBI-based pollinator coverage extraction was completed without opening flower-colour outcomes.

Across the 499 third-cohort plant species:

- all four query orientations technically valid: 499/499;
- any pollinator record: 423 species;
- >=3 distinct pollinator taxa: 356 species;
- >=2 primary pollinator guilds: 335 species;
- Sphingidae detected: 68 species;
- covered species without Sphingidae detection: 355 species;
- >=10 georeferenced pollinator records and >=2 one-degree cells: 295 species;
- of those, >=2 guilds: 275 species.

The live/API extraction had 180 truncated queries among 1,996 queries. It is therefore suitable for an exploratory screen but not a final stable-source inference.

## Strict exploratory join

To reduce obvious effort artefacts, the analysis retained only species with:

- all four pollinator queries technically valid;
- no truncated pollinator query;
- >=10 distinct pollinator taxa;
- >=20 usable third-cohort flower images after response-blind high-clip exclusion.

This yielded 217 species.

The image outcome was the species-level white proportion after high-clip exclusion. Models adjusted for mean near-clip fraction and standardized pollinator taxon richness as an effort term. A species-equal logit-transformed sensitivity was also inspected.

### P1 — Sphingidae association

Prediction: hawkmoth-associated plants should show more white.

Observed:
- 39/217 strict species had Sphingidae records;
- clustered image-weighted OR = 0.638;
- 95% CI = 0.374–1.088;
- p = 0.099;
- species-equal sensitivity p = 0.105.

The estimated direction is opposite the prediction and is not statistically supported.

### P2 — bee dominance

Prediction: a greater fraction of distinct pollinator taxa that are bees should be associated with less white.

Observed:
- OR per standardized bee-taxon fraction = 0.792;
- 95% CI = 0.624–1.006;
- p = 0.0555;
- species-equal sensitivity p = 0.154.

The direction matches the prediction but the result is weak and does not survive the species-equal sensitivity.

### P3 — Lepidoptera representation

Prediction: greater Lepidoptera representation should be associated with more white.

Observed:
- OR = 0.996;
- 95% CI = 0.788–1.259;
- p = 0.973;
- species-equal sensitivity p = 0.957.

This is effectively null.

## Stable bee-only coverage

A separate outcome-blind audit used the versioned Noori et al. (2026) curated GloBI bee-plant dataset (Zenodo 18303036 / GloBI_Curated.csv).

- FCP third-cohort species: 499;
- exact species-name matches: 346;
- coverage fraction: 69.3%;
- curated source size: 981,982 bee-plant interaction records;
- curated unique plant names: 12,699.

This passes the >=100-species coverage gate for a stable bee-only follow-up. It cannot by itself test hawkmoth, bird or bat alternatives.

## Current inference

Under the current global data:

- **hawkmoth association is not supported** and the point estimate is opposite the whitening prediction;
- **overall Lepidoptera representation is null**;
- **bee dominance shows only a weak tendency toward less white**, not a robust cross-species result.

Therefore there is currently no robust evidence that a single pollinator guild explains the recurrent white/non-white axis.

This does not rule out pollinator-mediated selection within particular systems. Global interaction data are strongly heterogeneous in sampling effort, geography and interaction semantics, and non-detection is not biological absence.

## Comparison with the environmental line

The pre-specified WorldClim test produced one positive signal: white states were modestly shifted toward warmer BIO5 environments across 281 species. That signal weakened under same-observer and local-distance restrictions.

The combined mechanism status is therefore:

> **temperature currently provides the stronger broad-scale sorting clue, while a global pollinator common-cause signal is not supported by the available interaction data. Neither result establishes pigmented -> white evolutionary causation.**

## Next test with the highest mechanistic value

A stronger pollinator test should be spatial and source-stable:

1. use a versioned pollinator interaction/occurrence release;
2. require adequate georeferenced pollinator sampling within each plant species;
3. estimate local pollinator assemblage around white and non-white flower records;
4. compare local bee / moth / bird guild composition within species;
5. test whether the BIO5 white-hot association remains after local pollinator composition is included.

This would discriminate climatic sorting from pollinator turnover rather than comparing coarse species-level labels.
