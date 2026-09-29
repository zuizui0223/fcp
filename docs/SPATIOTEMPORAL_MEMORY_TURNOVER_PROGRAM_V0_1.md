# General spatiotemporal memory–turnover programme — v0.1

## Decision

The next programme is **not a flower-colour generalization paper**.

Flower colour is retained only as the motivating proof of concept. The new biological question is:

> **Do the same biological properties predict how rapidly phenotypic and interaction information is lost across evolutionary time and geographic space?**

The programme explicitly includes both organismal traits and biotic interactions.

Its core 2 × 2 design is:

| Biological state | Evolutionary separation | Geographic separation |
|---|---|---|
| organismal trait | trait memory loss | within-species trait spatial turnover |
| biotic interaction | partner-profile evolutionary turnover | interaction rewiring across space |

All four cells use the same primary effect convention: a **null-centered distance–dissimilarity correlation**. Positive values mean that state dissimilarity rises faster with separation.

The machine contract is `data/spatiotemporal_memory_turnover_design_v0_1.json`.

## 1. Why this is more general than CHUN + FCP

CHUN and FCP already motivate two parts of the matrix:

- CHUN: flower-colour similarity usually declines as phylogenetic separation grows;
- FCP: flower-colour measurements can be geographically organized within species.

But the direct flower-colour bridge failed its outcome-blind coverage gate. That is useful: the general programme should not be built by squeezing the same flower-colour systems until a cross-scale correlation appears.

Instead, the new test starts from independent data resources and asks whether the architecture exists for **other biological traits and for interaction identity**.

CHUN/FCP therefore do not count as generalized replication evidence.

## 2. One estimator across four biological questions

For each admissible system define pairwise separation (x_{ij}) and state dissimilarity (d_{ij}).

The primary effect is

[
Delta ho = ho(x,d) - E_{null}[ho(x,d)].
]

Larger positive (Deltaho) means faster turnover with increasing separation.

This solves an important comparability problem. The programme does **not** compare:

- Blomberg K against Mantel r;
- phylogenetic signal p-values against spatial beta diversity;
- AUC against network modularity;
- null z-scores with different sample sizes.

It compares the same effect concept: **how quickly biological-state dissimilarity increases with separation beyond a frozen null**.

Null SD remains a design/informativeness quantity, not the biological effect size.

## 3. Trait arm

### 3.1 Evolutionary-time response

For one trait at a time:

1. aggregate repeated records to a predeclared species state for the temporal analysis;
2. compute pairwise trait dissimilarity;
3. compute relative patristic separation;
4. compare observed distance–dissimilarity association with tip-state permutations.

Unit: **clade × trait**.

The temporal effect is not called an absolute evolutionary rate unless an independently dated tree supports that interpretation.

### 3.2 Geographic-space response

Use the **same trait definition**, but retain repeated georeferenced observations.

Within each species:

1. pair observations by geographic distance;
2. calculate trait dissimilarity;
3. compare with complete-row permutations among observed coordinates.

Then aggregate species with equal weight to the clade × trait level.

This is intentionally a within-species spatial effect. It prevents ordinary species replacement from masquerading as trait spatial turnover.

## 4. Interaction arm

Interactions are treated as biological states, not merely as metadata attached to species.

### 4.1 Evolutionary-time interaction memory

For each focal guild, represent a taxon by its normalized partner profile.

The primary time response asks:

> **Do closely related focal taxa use more similar partners than expected after declared co-opportunity is respected?**

Partner-profile dissimilarity is Jensen–Shannon distance.

The null permutes complete partner profiles among focal taxa within the frozen opportunity stratum so profile degree/weight structure is retained.

### 4.2 Geographic interaction turnover

Total network turnover is **not** the primary response because it conflates two processes:

1. species disappear or appear between sites;
2. the same co-occurring species change partners.

The primary response is therefore the rewiring component among shared/available species.

The programme retains the classic decomposition:

[
eta_{WN} = eta_{OS} + eta_{ST},
]

where (eta_{OS}) captures changed interactions among shared species and (eta_{ST}) captures network differences forced by species turnover.

Only the rewiring component enters the primary time–space coupling analysis.

## 5. Same biological properties across time and space

Three predictor axes are frozen before generalized outcomes.

### Primary common driver — context heterogeneity

Greater abiotic or partner-opportunity heterogeneity should increase both temporal and spatial turnover.

The important claim is joint:

> **systems exposed to more heterogeneous contexts should lose state similarity faster with both evolutionary and geographic separation.**

The predictor must be frozen from source/exposure data before turnover outcomes are opened.

### Primary decoupler — ecological specialization

Specialization is not forced to have one sign.

The biologically interesting alternative is:

- strong specialization can preserve lineage-specific trait/partner structure through evolutionary time;
- the same specialization can increase spatial turnover when local opportunity changes.

Therefore an opposite-sign time/space coefficient is a **predeclared biological result**, not a failed common-driver test.

### Secondary decoupler — dispersal / mobility

Greater spatial mixing should weaken geographic turnover more strongly than evolutionary-time turnover.

Raw dispersal proxies are never pooled across incomparable taxa. They are standardized only within declared taxonomic groups.

## 6. Interaction type enters explicitly

Source-native interaction labels are frozen before outcomes:

- mutualism;
- herbivory;
- parasitism;
- predation / consumer-resource;
- other admissible categories only if defined by the source schema.

Interaction type is a moderator, not a source-selection criterion.

A prior comparative network study found stronger phylogenetic signal in antagonistic than mutualistic networks (Rohr & Bascompte 2014, DOI 10.1086/678234), so interaction type is biologically justified a priori.

However, the new programme does not assume that this historical result fixes the direction of geographic rewiring.

## 7. Why this is not already answered by existing literatures

Several pieces already exist separately.

- Phylogenetic niche/trait conservatism asks why some traits retain evolutionary history and others do not (e.g. Crisp & Cook 2012, DOI 10.1111/j.1469-8137.2012.04298.x).
- Network phylogenetic-signal studies show that evolutionary history can organize interaction architecture (Rohr & Bascompte 2014).
- Network beta-diversity frameworks separate species turnover from interaction rewiring (Poisot et al. 2012, DOI 10.1111/ele.12002).
- Spatial plant–pollinator studies document geographic interaction turnover and rewiring (Carstensen et al. 2014, DOI 10.1371/journal.pone.0112903).
- Short-term temporal studies show rapid interaction rewiring and effects of phenology/abundance (CaraDonna et al. 2017, DOI 10.1111/ele.12740).

The working gap is narrower:

> **the same response definition, same finite predictor family and same qualification logic have not yet been shown to explain turnover jointly across evolutionary time and geographic space for both traits and interactions.**

This is a working novelty hypothesis, not a priority claim. A formal database-level prior-art audit is required before manuscript priority language.

## 8. Source-first empirical programme

Current source order is frozen in `data/spatiotemporal_memory_turnover_source_matrix_v0_1.csv`.

### Trait discovery arm — BIEN 4.2

BIEN 4.2 is first because it integrates a large plant occurrence/plot/trait system under one taxonomic/geographic validation infrastructure.

Current public metadata report approximately:

- 25.9 million standardized trait observations;
- 54 standardized traits;
- 284 million occurrence records;
- 363,258 vegetation plots.

No trait values are opened for this programme until a metadata-only audit determines which traits actually satisfy both the temporal and repeated-geographic gates.

### Trait independent transport — AusTraits

AusTraits is retained as a separate open plant-trait transport source.

Its public portal exposes raw observations and georeferenced-record filtering. One exact versioned release must be pinned before row values are used.

### Trait high-power route — TRY v7

TRY v7 is attractive because it reports 22,860,407 records for 306,701 plant taxa and roughly half of its data are georeferenced.

It is currently **HOLD_SOURCE_ACCESS** because TRY reported a download/access-rights problem on 2026-09-23. The programme does not substitute scraped outputs while that source route is unstable.

### Interaction discovery arm — Mangal

Mangal is first because its schema directly exposes:

- datasets;
- replicated networks;
- dates;
- network geometry;
- taxonomic nodes;
- interaction types;
- optional traits and environment.

Its documentation reports 172 datasets and more than 1,300 networks.

The first execution is metadata-only: determine how many datasets have enough replicated georeferenced networks before any edge identities are analyzed.

### Interaction external transport — GloBI

GloBI provides stable versioned interaction products with source citations and taxonomic mappings through Zenodo.

It is intentionally second because global integration also creates stronger observation-process heterogeneity. Dataset/study identity and co-opportunity must be frozen before interaction records are analyzed.

## 9. Eligibility, not enthusiasm, controls source entry

Trait system gate:

- at least 20 temporal species;
- at least 8 species with repeated spatial measurements;
- at least 20 georeferenced records per spatial species;
- at least five occupied 25-km cells per spatial species.

Interaction system gate:

- at least 20 focal taxa for the temporal arm;
- at least three partners per focal taxon;
- at least five networks, all five georeferenced for the spatial arm;
- at least ten repeated focal taxa across the network series.

After structural eligibility, each candidate must pass a known-truth informativeness simulation before biological outcomes are interpreted.

The benchmark is frozen:

- median absolute recovery error <= 0.10 in (Deltaho);
- >=80% correct direction recovery when true (Deltaho=0.15).

Failure is a HOLD. Thresholds are not relaxed after seeing biological effects.

## 10. Possible biological outcomes

The programme is explicitly not built so that only positive coupling is interesting.

For each system, time and space turnover define four regimes:

| Time turnover | Space turnover | Interpretation |
|---|---|---|
| low | low | conserved / buffered |
| high | high | labile / context-tracking |
| low | high | historically retained but geographically contingent |
| high | low | evolutionarily labile but spatially homogenized |

The central question therefore becomes:

> **Which biological properties place systems into these turnover regimes?**

That is more general than asking whether one trait has phylogenetic signal or whether one network rewires across space.

## 11. Evidence sequence

The programme must proceed in this order:

```text
source identity
  -> schema/support eligibility
  -> observed-geometry informativeness
  -> freeze trait / interaction representation
  -> freeze biological predictors
  -> open four turnover responses
  -> test time-space coupling
  -> test common drivers / decouplers
  -> independent source transport
```

No step may be reordered because a downstream result looks promising.

## 12. Current boundary

This development lane changes none of the current FCP New Phytologist paper, CHUN Evolution Letters v0.3, CHUN hierarchical-memory extension or IWE meta-analysis.

It is a new cross-program hypothesis family.

If source qualification succeeds, it should be split into an independent repository before response-bearing generalized analyses become the main scientific object.
