# Spatiotemporal turnover coupling — source preflight result v0.1

## Decision

**Both source arms pass their response-blind physical coverage gates.**

This does not test the biological hypothesis yet. It establishes that the programme is not blocked at the source-availability stage.

## Trait arm — Tundra Trait Team

The archived cleaned v1 dataset contains **91,970 rows and 18 trait variables**.

The frozen paired-source gate required a trait to have:

- at least 20 species overall;
- at least 5 species each with >=10 usable observations;
- at least 3 distinct sites per spatially eligible species;
- at least 50 km maximum geographic span per spatially eligible species.

**14/18 traits pass.**

The four primary traits were frozen before any time/space turnover coefficient was computed:

| Trait | Species total | Spatially eligible species |
|---|---:|---:|
| SLA | 899 | 88 |
| LDMC | 754 | 64 |
| Vegetative plant height | 643 | 79 |
| Leaf N per dry mass | 399 | 34 |

These four represent complementary biological axes rather than four post-outcome alternatives.

## Interaction arm — Mangal

The current API v2 preflight recovered:

- **487 networks**;
- **487/487 with geometry**;
- **472/487 with dates**;
- **30 dataset IDs**;
- **38,939 nodes**;
- **144,989 interaction rows**.

The fixed primary interaction families also pass the pre-outcome repetition gate. The gate required at least 20 reference taxa each appearing in at least three georeferenced networks.

| Interaction family | Focal role | Repeated taxa (>=3 networks) |
|---|---|---:|
| Mutualism | either endpoint, undirected pool | **305** |
| Predation | Mangal directed source role (`node_from`) | **132** |
| Parasitism | Mangal directed source role (`node_from`) | **40** |

Herbivory is not promoted into the primary family after the coverage audit.

## What has not been opened

No analysis in this preflight calculated:

- a trait temporal-memory slope;
- a trait spatial-turnover slope;
- partner-profile phylogenetic memory;
- spatial rewiring;
- a time-space coupling statistic;
- an association with context heterogeneity, specialization or mobility.

The programme therefore remains outcome-unopened.

## Next gate

The next gate is **phylogeny + external predictor eligibility**.

For each arm, before any turnover result is opened:

1. construct lineage panels without response values;
2. pin an exact branch-length phylogeny source;
3. establish availability of P1 context heterogeneity;
4. establish availability of P2 specialization;
5. establish availability of P3 mobility/dispersal, or explicitly mark it NOT_EVALUABLE for that arm;
6. test estimator self-detectability.

Only after those conditions pass can the four turnover coordinates be estimated.

## Relation to CHUN/FCP

This development lane is independent of the active FCP manuscript and frozen CHUN EL v0.3.

CHUN and FCP motivated the two distance axes. They are not used as positive observations proving the new coupling hypothesis.

