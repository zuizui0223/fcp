# Spatiotemporal turnover coupling — development programme v0.1

## Main question

> **Do the same biological properties predict how rapidly biological information is lost across evolutionary time and geographic space, and does the answer generalize from phenotypic traits to species interactions?**

This is deliberately broader than flower colour.

Flower colour is retained only as an originating proof of concept through CHUN and FCP. It is not the response domain of this programme.

## Four-response architecture

The framework uses the same distance-decay logic twice for traits and twice for interactions.

| State | Evolutionary time | Geographic space |
|---|---|---|
| phenotype / trait | trait temporal memory loss | within-species trait spatial turnover |
| interaction partner profile | partner temporal memory loss | spatial rewiring among local networks |

Every eligible lineage panel therefore has up to four turnover coordinates.

A missing coordinate is `NOT_EVALUABLE`, never zero.

## Why two arms first

Requiring one dataset to contain phylogenies, repeated individual traits, replicated local interaction networks and external biological predictors would make source availability the dominant scientific result.

The programme therefore starts with two independent arms:

1. **Trait arm:** same trait and lineage panels must support both temporal and spatial estimands.
2. **Interaction arm:** same focal lineages must support both partner-memory and spatial-rewiring estimands.

Only panels that independently qualify all four cells enter the later trait–interaction alignment test.

This prevents a weak cell from being fabricated by proxy substitution.

## Primary biological predictors

Only three predictor families are primary.

### P1 — context heterogeneity: predicted coupler

Environmental and partner-context heterogeneity should expose lineages to more heterogeneous selective and ecological worlds.

Frozen prediction:

> higher context heterogeneity → faster temporal loss and faster spatial turnover.

This is the main common-cause hypothesis.

### P2 — specialization: predicted decoupler

Specialization can preserve a lineage-specific phenotype or partner set through evolutionary time while making local realization sensitive to whether the required habitat/resource/partner is available.

Therefore the prediction is intentionally asymmetric rather than a generic positive correlation.

### P3 — mobility / dispersal: predicted decoupler

Higher movement capacity should homogenize geographic states and interactions more directly than it constrains evolutionary lability.

The prediction is therefore strongest and negative on the spatial axes.

## Central novelty

The programme does **not** ask whether phylogenetic signal, spatial autocorrelation or interaction beta diversity exist. All three are established topics.

The new question is whether the **same external biological properties explain the rate at which information decays along two fundamentally different distance axes**, and whether that coupling changes when the state is an interaction rather than a phenotype.

This makes weak time–space correlation informative if it is explained by predeclared decouplers.

## Interaction-specific safeguard

Raw interaction-network turnover is not the target.

For the spatial interaction arm:

`total network turnover = partner/species turnover + rewiring conditional on shared opportunity`

The primary response is the rewiring component.

If two sites differ only because different partner species are present, that does not count as evidence that the focal species changed partner choice.

Likewise, GloBI or other aggregate record counts are never interpreted directly as preference or interaction strength without an effort model.

## Initial source programme

### Trait feasibility — Tundra Trait Team

The Tundra Trait Team is the strongest first feasibility candidate because its public v1 database contains roughly 92k measurements across 978 vascular plant species; more than 99% of records are georeferenced and most measurements are individual-level.

The first audit is schema/coverage only. It must count how many lineage panels can support both:

- a branch-length phylogenetic temporal test; and
- at least five species with >=10 observations across >=3 sites and >=50 km geographic span for the same trait.

No trait is selected because it gives a strong time or space result.

### Interaction feasibility — Mangal, then GloBI

Mangal is the first candidate for repeated observed local networks because networks are explicit sampling units.

GloBI is a broader transport source because it exposes spatial-temporal interaction records across many interaction types, but raw aggregate occurrence density creates a larger effort problem. GloBI therefore cannot become the first primary interaction arm unless local sampling units can be reconstructed prospectively.

## Eligibility before outcomes

The order is fixed:

1. source identity/version;
2. schema and taxonomic crosswalk;
3. response-blind lineage panels;
4. coverage for both time and space within an arm;
5. external predictor availability;
6. null-model self-detectability;
7. outcome opening.

A source that fails a gate ends in HOLD/STOP. Thresholds are not relaxed after turnover values are seen.

## Relation to CHUN and FCP

CHUN and FCP motivated the two distance axes:

- CHUN: evolutionary proximity can retain phenotype information without one universal privileged resolution.
- FCP: geographic proximity can retain phenotype information without one universal geographic transition map.

The direct CHUN–FCP bridge itself is currently on coverage HOLD, so those two papers are not used as evidence that time and space are already coupled.

This development programme tests that coupling independently, outside the active papers.

## Current next actions

1. run Tundra Trait Team source/schema/coverage preflight;
2. run Mangal network metadata preflight;
3. freeze dated phylogeny sources only after the response-blind panels are known;
4. admit no biological turnover result until an arm passes its paired time+space coverage gate.

## Governance

- Active FCP New Phytologist science is unchanged.
- Frozen CHUN Evolution Letters v0.3 science is unchanged.
- This is not a TTF methods result.
- This is not an IWE meta-analysis extension.
- No claim of a universal coupling is made before independent trait and interaction arms both qualify.
