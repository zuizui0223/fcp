# FCP white–pollinator independent coverage protocol — 2026-09-25

Status: **coverage-only audit frozen before any pollinator variable is joined to white/non-white outcomes.**

This is a separate mechanism line. It does not modify the frozen New Phytologist manuscript or the prospective white-environment result.

## Question

Can an independent plant–pollinator interaction source support a comparative test of whether recurrent white flower states are associated with pollinator guilds, especially nocturnal Lepidoptera/hawkmoths, rather than—or in addition to—broad-scale temperature sorting?

This stage tests **data identifiability only**. It must not read flower-colour outcomes, D, H2 vectors, or environmental-result values.

## Plant frame

Use only the species identities in the frozen third-cohort prospective artifact (499 authorized species). Read the species column only and deduplicate it before querying interactions.

## Interaction source

Coverage discovery uses the current GloBI Web API only as an exploratory index.

Query both semantic orientations:

- pollinator -> plant: `pollinates`, `visitsFlowersOf`;
- plant -> pollinator: `pollinatedBy`, `flowersVisitedBy`.

The API is not a final inferential data source because GloBI explicitly describes live/API results as unstable. If coverage passes, the biological analysis must be rebuilt from a stable versioned GloBI integrated data release (or original source datasets) and archive the exact version/digest.

## Pollinator guild classification

Classify from returned taxonomic paths, not common names.

Primary guilds:
- bee: Apoidea/Anthophila or recognized bee families;
- other Hymenoptera;
- Lepidoptera;
- hawkmoth: Sphingidae subset of Lepidoptera;
- Diptera;
- Coleoptera;
- Aves;
- Chiroptera;
- other/unknown.

No record is converted to pollinator absence. Missing coverage remains missing.

## Frozen coverage gates

### P0 — query integrity
Pass only if:
- >=95% of the 499 plant names return a technically valid API response for all four query orientations, or a valid zero-row response;
- raw query receipt records query URL parameters, HTTP status, row count, and whether the per-query limit was reached.

If P0 fails, stop: `coverage_not_evaluable_api_integrity`.

### P1 — broad pollinator identifiability
Pass only if:
- >=100 plant species have >=1 independent pollination/flower-visitation interaction;
- >=60 plant species have >=3 distinct pollinator taxa.

If P1 fails, no cross-species pollinator mechanism test is authorized.

### P2 — multi-guild identifiability
Pass only if:
- >=40 plant species have records from >=2 primary pollinator guilds among bee, other Hymenoptera, Lepidoptera, Diptera, Coleoptera, Aves, Chiroptera.

If P2 fails, do not fit a generic pollinator-guild diversity model.

### P3 — hawkmoth/nocturnal-lepidopteran identifiability
For the specifically predicted white–hawkmoth mechanism, require:
- >=20 plant species with at least one Sphingidae interaction;
- AND >=20 comparison species with pollinator coverage but no Sphingidae interaction.

The second group is **documented non-detection**, not biological absence. Even if P3 passes, subsequent inference must model/weight interaction-documentation effort.

### P4 — spatial identifiability
For a within-species geographic mechanism test, require:
- >=50 plant species with >=10 georeferenced interaction records;
- >=2 distinct 1-degree occupied interaction cells per retained plant;
- at least 30 retained plants with >=2 pollinator guilds across their georeferenced records.

If P4 fails, species-level guild association may be explored only as a mechanism clue; it cannot explain the within-species geographic sorting of white states.

## Truncation rule

Each live API query is capped only for the coverage audit. Any query reaching the cap is marked `truncated=true`. Truncated records can establish presence and minimum richness but cannot be used for abundance, dominance, guild proportions, or final inference.

## Hard nonclaims

A coverage pass does not establish:
- pollinator preference for white flowers;
- pollinator-driven evolution;
- pigmented -> white transition direction;
- hawkmoth causation;
- pollinator abundance or dominance;
- independence from climate.

The next biological protocol is allowed only after this coverage result is opened and only for dimensions that meet their frozen coverage gate.
