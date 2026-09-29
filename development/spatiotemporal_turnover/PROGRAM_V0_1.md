# Spatiotemporal turnover of traits and biotic interactions — development v0.1

## Current decision

This is a **new general development lane**, not an extension of the active FCP manuscript and not a modification of frozen CHUN Evolution Letters v0.3 science.

Flower colour is retained only as one motivating proof of concept.

The general question is:

> **Do the same biological properties predict how rapidly trait and interaction information is lost across evolutionary time and geographic space?**

The programme deliberately includes both **traits** and **biotic interactions**. It is therefore broader than flower colour, pollination, or any one taxonomic group.

## The four-response architecture

All primary responses use one orientation:

> **larger response = faster turnover / faster loss of similarity with separation**

| | Evolutionary time | Geographic space |
|---|---|---|
| **Trait** | **TT** — trait dissimilarity versus phylogenetic distance | **TS** — within-species trait dissimilarity versus geographic distance |
| **Interaction** | **IT** — partner-profile dissimilarity versus phylogenetic distance | **IS** — local-network rewiring versus geographic distance |

The common rank-scale primitive is:

[
\rho_{turnover}=\operatorname{Spearman}(distance, biological\ dissimilarity).
]

This is not claimed as a new statistic. The novelty target is the **joint biological prediction across axes and data types**.

### TT — trait temporal turnover

For one radiation or biological system:

- retain at least 20 focal taxa;
- require a branch-length phylogeny;
- require at least 80% coverage of a predeclared trait representation;
- calculate pairwise trait dissimilarity;
- associate trait dissimilarity with phylogenetic distance.

CHUN flower-colour memory is a motivating special case, not the general data source.

### TS — trait spatial turnover

Within each focal species:

- require at least 40 georeferenced trait observations;
- calculate pairwise geographic distance and trait dissimilarity;
- estimate one species-specific turnover rho.

A system-level response is the median species response across at least five focal species.

This generalizes the FCP/disttrait spatial architecture beyond flower colour.

### IT — interaction temporal turnover

Each focal taxon is represented by a partner profile.

The partner-profile representation and the metaweb / partner-opportunity definition must be frozen before response calculation.

The response asks:

> Do phylogenetically distant focal taxa use increasingly different partners?

Partner profiles may be binary or quantitative, but the representation must be fixed within a system.

### IS — interaction spatial turnover

For replicated local networks, whole-network dissimilarity is decomposed into:

[
\beta_{WN}=\beta_{ST}+\beta_{OS}.
]

- **beta_WN** — whole-network dissimilarity;
- **beta_ST** — component attributable to species turnover;
- **beta_OS** — rewiring among species shared between networks.

The **primary** spatial interaction response is beta_OS versus geographic distance.

beta_ST is mandatory, because a network cannot be described as rewired merely because one partner disappeared from the local species pool.

The repository implementation uses the additive common-denominator partition. The decomposition itself is prior art; it is not a novelty claim.

## Three frozen predictor families

### AB — abiotic niche breadth

Outcome-independent occurrence data define multivariate climatic breadth.

Primary prediction:

> **broader abiotic breadth buffers turnover on all four responses.**

So the frozen direction is negative for TT, TS, IT and IS.

This is the cleanest direct coupling hypothesis.

### BS — biotic specialization

Specialization is defined from an independent metaweb, external interaction source, or frozen training period.

It may **not** be calculated from the same local networks used to define the focal spatial-rewiring response unless a disjoint training split was frozen first.

Primary prediction:

> **specialization is a temporal–spatial decoupler.**

Frozen directions:

- TT: negative;
- IT: negative;
- TS: positive;
- IS: positive.

The biological idea is that narrow ecological dependence can preserve a lineage-specific state through time while increasing sensitivity to local changes in partner opportunity across space.

### MO — mobility / dispersal

Mobility is a guild-specific dispersal proxy fixed before outcomes and standardized within guild.

Frozen prediction:

- TS: negative;
- IS: negative;
- TT and IT: no directional primary prediction.

Mobility is therefore a **spatial-control / decoupling axis**, not a universal memory predictor.

## Interaction type

Interaction type is not allowed to become a fourth predictor after results are opened.

It is a predeclared moderator with levels:

- mutualistic;
- antagonistic;
- trophic;
- host-parasite / parasitoid;
- other prespecified.

The primary use is two-sided heterogeneity of predictor effects among interaction types.

No universal winner is assumed.

## Main hypotheses

### H1 — common abiotic buffering

The same broad abiotic niche should predict slower turnover in both time and space.

This is tested separately for:

- TT versus TS;
- IT versus IS.

### H2 — specialization decoupling

Specialization should produce opposite temporal and spatial effects:

[
\beta_{time,BS}<0, \qquad \beta_{space,BS}>0.
]

This is tested for both traits and interactions.

### H3 — mobility decoupling

Mobility should reduce geographic turnover for both traits and interactions while receiving no directional temporal prediction.

### H4 — biotic moderation

The magnitude of these effects may differ among interaction types, but that heterogeneity is tested rather than assumed.

## Why this is not already answered by existing network ecology

Several pieces are established prior art:

- interaction-network beta diversity across space/time/environment;
- decomposition of whole-network dissimilarity into species-turnover and rewiring components;
- phylogenetic signal in partner use and network architecture;
- trait-based network assembly;
- traits and environment as drivers of spatial network turnover;
- network rewiring under environmental variation.

The new target is narrower and more synthetic:

> **whether the same predeclared biological properties predict turnover across the trait × interaction and time × space matrix, and whether those predictor effects couple or decouple the axes.**

Do not claim invention of phylogenetic signal, network rewiring, interaction beta diversity, or trait-based network ecology.

## Source strategy

The source registry is frozen separately in
`development/spatiotemporal_turnover/source_registry_v0_1.csv`.

Current status matters:

- TRY v7 is a high-value plant-trait candidate, but its official site currently reports a requested-data download outage. It is therefore a source HOLD, not an assumed available input.
- GloBI is useful as a stable, versioned interaction metaweb/provenance layer, but integrated rows are not treated as one homogeneous survey.
- Web of Life is a candidate for geolocated local interaction networks; file-level identities and original-study provenance must be frozen before use.
- archived GlobalWeb food webs are a candidate trophic source, again conditional on source-study comparability.
- CHUN and FCP are proof-of-concept anchors only and do not count as general-law replications.

## Admission gates

No confirmatory biological response may be opened until all relevant gates pass.

1. **Source identity** — exact version/file/source-study identity.
2. **Response geometry** — minimum taxon/site/observation support and nondegenerate distances.
3. **Predictor independence** — no predictor may reuse its focal response surface without a frozen disjoint training design.
4. **Cross-axis pairing** — direct time–space claims require the same biological system identity on both axes.
5. **Known-truth calibration** — false-positive control and detectability must be calibrated for the planned system count before outcome opening.
6. **Stop rule** — a failed gate is a HOLD, not permission to reduce support thresholds.

## Current empirical status

**No biological outcome has been opened for this general programme.**

The current state is:

`DESIGN_READY -> SOURCE_AUDIT -> KNOWN_TRUTH_CALIBRATION -> ONLY THEN BIOLOGICAL OPENING`

This is intentionally separate from:

- the active FCP New Phytologist paper;
- the frozen CHUN Evolution Letters v0.3 submission;
- the failed direct CHUN–FCP flower-colour bridge, which remains a coverage HOLD.

## What would count as a strong result

A strong positive result is not “all four turnover coefficients are the same.”

Three different outcomes are biologically meaningful:

1. **coupling** — the same predictor has the same directional effect in time and space;
2. **decoupling** — a predictor has opposite temporal and spatial effects;
3. **domain specificity** — a predictor matters for traits but not interactions, or vice versa.

The framework is therefore compatible with the CHUN/FCP lesson that biological structure can be recurrent without one universal global template.
