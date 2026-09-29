# General spatiotemporal turnover programme — v0.1

## Decision

This development lane generalizes the CHUN–FCP sister-paper insight beyond flower colour.

The target is not a universal flower-colour rule. It is a general ecological-evolutionary question:

> **Do the same biological properties predict how rapidly biological information is lost across evolutionary time and geographic space, for both organismal traits and biotic interaction identities?**

The programme is deliberately separated from the active FCP New Phytologist paper and from frozen CHUN Evolution Letters v0.3.

## 1. Core object: a four-dimensional turnover signature

For each biologically coherent radiation/system, define four loss rates:

- lambda_trait_time
- lambda_trait_space
- lambda_interaction_time
- lambda_interaction_space

Larger lambda always means faster loss of similarity / faster turnover.

The point is not to force the four rates to be equal. The scientific object is the **shape of the four-dimensional turnover signature**.

### Trait × evolutionary time

How rapidly does trait similarity above a frozen null decline with relative phylogenetic divergence?

CHUN flower-colour memory is one proof-of-concept instance of this coordinate, not its definition.

### Trait × geographic space

How rapidly does the same trait representation turn over across geography using repeated population- or individual-level measurements?

FCP continuous colour organization is one proof-of-concept instance of this coordinate.

### Interaction × evolutionary time

For each focal species, represent interaction identity as a normalized partner profile. Pairwise similarity is based on 1 - Jensen-Shannon divergence.

The primary question is how quickly partner-profile similarity declines with focal-species phylogenetic distance.

This comparison must preserve or explicitly model partner opportunity. A phylogenetic signal caused only by related species occupying the same partner pool is not accepted as interaction-memory evidence.

### Interaction × geographic space

Measure interaction beta diversity across spatially replicated networks.

The programme **requires decomposition** into:

1. species-turnover contribution; and
2. rewiring among species that co-occur in both networks.

The rewiring component is the primary interaction-space response.

Species replacement is never relabelled as partner switching.

## 2. Three predictors only

The first programme version permits exactly three biological predictor axes.

### P1 — context heterogeneity

Role: **shared accelerator**.

This is a pre-outcome equal-weight composite of:

- abiotic environmental heterogeneity; and
- partner-opportunity heterogeneity.

Both components must be present. Missing one component does not authorize renormalization to the other.

Primary directional prediction:

| Response | prediction |
|---|---:|
| trait × time | faster loss (+) |
| trait × space | faster turnover (+) |
| interaction × time | faster loss (+) |
| interaction × space | faster rewiring (+) |

This is the cleanest test of a common time-space driver.

### P2 — specialization

Role: **interaction conservatism / vulnerability decoupler**.

Specialization is frozen from an independent resource/partner opportunity frame, not estimated from the same realized link outcome later used as the response.

Predictions:

- slower evolutionary loss of partner identity;
- less rewiring among co-occurring partners;
- potentially more interaction loss through partner/species replacement when the required partner pool turns over.

Thus specialization predicts a decomposition, not simply more or less total interaction turnover.

### P3 — dispersal / mobility

Role: **spatial buffer / time-space decoupler**.

The mobility axis must come from an independently chosen biological measure before response inspection.

Predictions:

- weaker trait spatial turnover;
- weaker interaction spatial turnover;
- no required directional effect on either temporal-memory response.

A successful P3 result would therefore explain why temporal and spatial stability can decouple.

## 3. Interaction type is a fixed stratum, not another predictor

The programme freezes three interaction strata:

- mutualism;
- antagonism;
- mixed / other.

This is necessary because previous work already establishes that phylogenetic signal in network structure can differ between antagonistic and mutualistic networks.

The programme does **not** claim priority for phylogenetic signal in ecological networks.

The new question is whether the **coupling between temporal memory and spatial turnover**, and its three predictors, differ among interaction types.

## 4. Novelty boundary

The following are established areas and are not claimed as new:

- phylogenetic signal in organismal traits;
- phylogenetic signal in ecological network structure;
- spatial beta diversity;
- interaction beta diversity;
- network rewiring;
- specialization as an ecological axis;
- dispersal as a spatial homogenizing process.

The proposed contribution is their combination into one auditable response architecture:

> **the same biological system is represented by paired time- and space-loss rates, for both traits and interaction identities, and the same predeclared biological predictors are asked to explain aligned versus decoupled turnover regimes.**

A second distinctive feature is that interaction-space turnover is never interpreted without separating species replacement from rewiring.

## 5. Why a linked-system design is necessary

No single standard database currently supplies all required objects with adequate sampling semantics.

The programme therefore permits **linked systems**:

trait source + interaction source + phylogeny + spatial opportunity + predictor source

but only when the sources refer to the same focal radiation under a frozen taxonomic crosswalk.

A linked system is one inferential unit. It does not become five independent replications because five databases were joined.

## 6. Current source-first eligibility result

Canonical matrix:

data/general_spatiotemporal_turnover_eligibility_matrix_v0_1.csv

### Plants: TRY × GloBI / METRIN-KG

This stack is promising because plant traits and pairwise interactions can be linked taxonomically, and TRY includes large numbers of georeferenced trait records.

However, pairwise GloBI records do not automatically provide standardized spatial network effort or partner-opportunity semantics.

Current status: **HOLD_NOT_FULL_LINKED**.

### Mangal-style multi-network systems

Multi-network repositories are promising for spatial interaction turnover because local networks can retain network identity and location.

They still require a frozen inventory of genuinely independent programmes, external phylogeny, external repeated trait data, and an independent mobility axis.

Current status: **HOLD_SYSTEM_INVENTORY_REQUIRED**.

### Rohr & Bascompte 60-network data

The 60 antagonistic/mutualistic networks are important prior art for interaction × evolutionary history.

They are not a time-space candidate because the study was not designed as matched spatial replication.

Current status: **PRIOR_ART_NOT_CANDIDATE**.

### Birds: AVONET-linked route

AVONET is unusually strong for broad phylogenetic morphology, individual specimen measurements, and hand-wing index / migration as a pre-existing mobility axis.

But a fine-scale repeated spatial trait design is not established merely by country/range metadata, and standardized spatial interaction effort is not supplied by AVONET.

Current status: **HOLD_NOT_FULL_LINKED**.

### NEON small mammals

NEON provides repeated individual measurements and a strong spatial sampling structure, making it a useful trait-domain test bed.

It lacks a matched interaction network for the same units.

Current status: **TIER_B_TRAIT_ONLY**.

### CHUN × FCP flower colour

This remains the motivating proof of concept.

It is not used as the confirmatory general test. The direct overlap gate already terminated as insufficiently covered.

Current status: **HOLD_DIRECT_LINK_COVERAGE**.

### IWE phenological interactions

IWE addresses synchrony and reproductive consequences.

Phenological synchrony is neither phylogenetic interaction-memory loss nor geographic network rewiring, so it is not relabelled to fill an empty response cell.

Current status: **NOT_ELIGIBLE_AS_CORE**.

## 7. Programme admission gate

A general-law analysis cannot open until both domains independently qualify.

### Trait domain

Require at least:

- 12 independent systems;
- 3 major taxonomic groups;
- both temporal and spatial trait responses in every admitted system;
- predictor values frozen before response fitting.

### Interaction domain

Require at least:

- 12 independent systems;
- at least 2 interaction types;
- both temporal partner-memory and spatial interaction turnover;
- explicit species-turnover / rewiring decomposition;
- predictor values frozen before response fitting.

If either domain fails, the programme terminates as HOLD rather than lowering the denominator.

## 8. Prespecified biological hypotheses

### T1 — shared context acceleration

More heterogeneous abiotic + partner-opportunity context predicts faster information loss on all four axes.

This is the primary same-property-predicts-time-and-space hypothesis.

### T2 — time-space coupling

Within trait and interaction domains separately:

> systems with faster evolutionary memory loss also show faster spatial turnover.

This is not assumed to be universal. It is tested after P1–P3 are frozen.

### T3 — specialization decoupling

Specialization predicts stronger temporal partner memory, lower rewiring among co-occurring species, and greater vulnerability to partner/species replacement.

### T4 — mobility decoupling

Greater dispersal/mobility buffers spatial turnover but need not change temporal memory.

### T5 — interaction-type modification

Mutualistic and antagonistic systems may differ in T2/T3 coupling.

The direction is not selected after outcomes.

## 9. Turnover regimes

The programme explicitly allows four qualitative regimes rather than treating non-coupling as failure.

| temporal loss | spatial turnover | regime |
|---|---|---|
| slow | slow | conserved / buffered |
| fast | fast | labile / context-tracking |
| slow | fast | historically conserved, spatially contingent |
| fast | slow | evolutionarily labile, spatially homogenized |

The third and fourth regimes are biologically important evidence for decoupling.

## 10. Finite next source search

No open-ended candidate fishing is allowed.

The next source-only search is capped at six programme classes:

1. multi-site plant–pollinator networks linked to repeated plant traits;
2. multi-site plant–herbivore networks linked to repeated plant traits;
3. host–parasite systems with repeated spatial networks;
4. bird resource/interaction systems linked to AVONET mobility;
5. mammal parasite/resource systems linked to repeated individual traits;
6. aquatic food-web systems with repeated local networks and independent mobility/body-size traits.

For each class, inspect only source identity, taxon identifiers, network/site counts, location/date/effort schema, trait field schema, phylogeny availability, and predictor availability.

Do **not** compute any turnover response during source qualification.

Stop when 12 qualifying systems exist in each response domain or when the six-class source family is exhausted.

## 11. Relation to existing repositories

- CHUN supplies a temporal-memory proof of concept.
- FCP supplies a spatial-organization proof of concept and reusable disttrait spatial machinery.
- IWE supplies interaction-type biological context but not automatic response rows.
- TTF supplies abstention philosophy and qualification discipline, not the ecological claim.

This development lane does not alter any of those papers.

## Current decision

**No Tier-A full linked system is promoted yet.**

That is a source-architecture result, not a biological null.

The next legitimate action is the six-class source-only eligibility search above.
