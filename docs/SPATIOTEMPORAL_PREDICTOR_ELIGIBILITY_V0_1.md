# Spatiotemporal turnover predictor eligibility — v0.1

## Current decision

The generalized programme now separates **response definition** from **predictor acquisition**.

The goal is not merely to ask whether time and space turnover covary. It is to test whether the **same predeclared biological properties** explain both axes without recycling the same biological measurements into predictor and response.

The three predictor families remain unchanged:

1. **abiotic context heterogeneity** — common driver;
2. **biotic partner specialization** — primary decoupler;
3. **direct dispersal / mobility capacity** — secondary decoupler.

No fourth predictor is added.

## 1. Why source independence matters

A predictor can look powerful for purely mathematical reasons if it is derived from the same rows used to define turnover.

Prohibited examples include partner richness from the same Mangal edges used to calculate rewiring, climatic heterogeneity selected using focal-trait outcomes, or relabelling height/seed mass as mobility after seeing results.

Canonical contract: data/spatiotemporal_predictor_independence_contract_v0_1.json

Canonical 4-response × 3-predictor matrix: data/spatiotemporal_predictor_source_matrix_v0_1.csv

## 2. Predictor 1 — abiotic context heterogeneity

Primary source stack: GBIF occurrence geometry plus CHELSA v2.1 climate.

The primary quantity is a frozen climate-dispersion score based on median pairwise Euclidean distance in standardized climate space.

Frozen dimensions: BIO5, BIO10, BIO12, BIO15, and climate moisture index.

No alternative climate PC set or variable subset may replace this after turnover outcomes are opened.

For trait systems, species receive equal weight. For interaction systems, the score is computed from response-independent network-site geography.

Expected direction: TT + / TS + / IT + / IS +.

## 3. Predictor 2 — biotic partner specialization

The biological property is partner specialization rather than network degree derived from the response dataset.

Primary source: stable versioned GloBI interaction archive.

Primary definition: specialization = 1 - normalized Shannon entropy(partner distribution), with higher values meaning narrower partner use.

Minimum support: 10 partner records and at least 2 distinct partners.

For trait responses, GloBI is external to BIEN/AusTraits focal trait values.

For Mangal interaction responses, GloBI records are admissible only after every record traceable to the same Mangal dataset, source DOI, or source namespace has been excluded.

If that exclusion cannot be demonstrated, the test terminates HOLD_SPECIALIZATION_SOURCE_NOT_INDEPENDENT.

## 4. Predictor 3 — dispersal / mobility

This predictor is intentionally not yet qualified for the primary plant panel.

Current terminal status: HOLD_NO_PINNED_SOURCE_INDEPENDENT_DIRECT_DISPERSAL_MEASURE.

Post-outcome rescue is prohibited: do not substitute plant height, seed mass, or an optimized trait PC after seeing turnover outcomes.

Mobility is secondary, so its HOLD does not block context-heterogeneity or specialization tests.

## 5. Trait source calibration is separate from primary evidence

A closed development branch performed a physical-support audit on the Tundra Trait Team archive.

Calibration support: 91,970 rows; 18 traits; 14 passed the paired temporal/spatial support gate. SLA had 899 total species / 88 spatially eligible; LDMC 754 / 64; vegetative height 643 / 79; leaf N 399 / 34.

Canonical receipt: data/spatiotemporal_trait_geometry_calibration_tundra_v0_1.json

This remains calibration only and does not expand the frozen primary trait-source family.

## 6. Interaction response and predictor separation

Current Mangal source qualification found no fresh full time-space interaction system suitable for generalized biological inference.

Havens has structurally unidentified rewiring; Hadfield is retrospective calibration; Ricciardi is space-usable but has temporal focal-guild support below 20.

GloBI is being audited separately as a specialization-predictor source. The first audit opens only dataset namespaces, not interaction rows or taxon identities.

## 7. Four-response × three-predictor eligibility

| Response | Context heterogeneity | Partner specialization | Mobility |
|---|---|---|---|
| trait time | eligible source route | pending GloBI taxon/support audit | HOLD |
| trait space | eligible source route | pending GloBI taxon/support audit | HOLD |
| interaction time | eligible source route | HOLD pending source-overlap exclusion | future transport only |
| interaction space | eligible source route | HOLD pending source-overlap exclusion | future transport only |

This is a source/independence state, not a biological result.

## 8. Prior-art boundary

Geographic distance-decay, functional/phylogenetic turnover, network rewiring, and interaction beta diversity are established component literatures.

The narrower novelty target is: **one common turnover effect convention, one finite biological predictor family, and one response-independence protocol applied jointly to trait and interaction states across evolutionary and geographic separation.**

This remains a working novelty hypothesis until a formal database-level prior-art audit is completed.

## 9. Next executable gate

1. finish the GloBI specialization source-overlap audit;
2. freeze the primary plant trait panel from BIEN or AusTraits without using turnover outcomes;
3. pin a branch-length phylogeny and tip-quality rule;
4. run known-truth detectability on observed support geometry;
5. only then open TT/TS;
6. keep interaction IT/IS closed until a source passes the fresh paired-system gate.

## Governance

- generalized trait/interaction turnover outcomes remain unopened;
- active FCP New Phytologist science is unchanged;
- CHUN Evolution Letters v0.3 is unchanged;
- CHUN/FCP flower-colour results remain motivating anchors only;
- no predictor family may be replaced after its corresponding generalized response is opened.
