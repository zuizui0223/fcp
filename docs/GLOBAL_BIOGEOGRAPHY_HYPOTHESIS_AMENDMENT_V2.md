# Global flower-colour biogeography — hypothesis amendment v2

Status: **prospective amendment frozen on 2026-09-04 while the global species-discovery workflow was still metadata-only and before any global Monte Carlo colour field was available.**

This amendment does not change the frozen H1–H6 results. It reorganizes the next analysis around recurrent colour geography, species-specific heterogeneity, ecological/geographic concordance and between-species community assembly.

Where this document conflicts with the first wording of §10.1 in `GLOBAL_MONTE_CARLO_BARRIER_ATLAS_PROTOCOL.md`, this amendment and `global_edge_mechanism_distance_control_amendment_v1.json` are authoritative: **the primary G3 edge mechanism statistic is distance-controlled partial Spearman, not raw edge Spearman.** Raw edge Spearman is descriptive sensitivity only.

## 1. Current inferential ladder

### G1 — recurrent flower-colour geography

Ask whether adequately sampled species place strong within-species flower-colour discontinuities in the same broad regions more often than expected under species-conditioned colour permutation.

### G2 — stable neutral colour zones

If G1 is supported, require resampling stability and extract persistent zones mechanically as `Z01`, `Z02`, ... before consulting named geography or external ecology.

### G3 — species-specific mechanism concordance

G3 is evaluable even if G1 is null. On fixed within-species graph edges, test whether colour discontinuity tracks external turnover **beyond generic geographic separation**:

`partial Spearman(colour-edge contrast, external-edge contrast | great-circle edge distance)`.

The fixed primary external family is now:

1. direct GloBI pollinator-interaction turnover;
2. independent GBIF pollinator-guild turnover;
3. CHELSA multivariate climate turnover;
4. **SoilGrids multivariate edaphic turnover**;
5. EarthEnv terrain turnover;
6. marine-gap fraction;
7. major-river crossing;
8. RESOLVE biogeographic-boundary crossing;
9. GMBA mountain-system boundary crossing.

The same 999 within-species colour permutations are reused and Holm correction is applied over the fixed primary family.

### G4 — interpretation of supported shared zones

Only after G1+G2 pass may the neutral zones receive non-causal descriptors such as `pollinator-associated`, `climate-associated`, `edaphic-associated`, `terrain-associated` or geographic-barrier-associated.

### G5 — sympatry versus allopatry colour assembly

G5 is a new orthogonal question: **do species that occur sympatrically have more similar or more divergent flower colours than comparable species that do not occur sympatrically?**

G5 uses independent GBIF plant occupancy to define geography and the same location-blind species colour profiles already measured for the atlas. It does not require or rescue G1.

For each focal × sympatric partner pair, construct an outcome-blind matched set of allopatric controls matched on taxonomic distance class, dominant realm, climate centroid, range occupancy and GBIF effort. The primary matched-set contrast is:

`delta = D_colour(sympatric) - mean D_colour(matched allopatric controls)`.

- `delta < 0`: colour convergence in sympatry;
- `delta > 0`: colour divergence in sympatry.

The primary global test is two-sided because both are biologically plausible. The null randomizes which member of each frozen matched partner set is labelled pseudo-sympatric. Focal species are then equal-weighted.

Species can additionally receive FDR-controlled descriptive community-embedding labels:

- `convergent_in_sympatry`;
- `divergent_in_sympatry`;
- `undetected`.

These are phenotypic community descriptors, not evolutionary diagnoses.

## 2. Reassessment of the previous mechanism hypotheses

The older C/S framework remains useful as biological theory, but it should no longer be the sole response variable. Its hypotheses are redistributed as follows.

### Promoted into the new atlas

**Pollinator geographic mosaic.** This is now a direct G3/G4 candidate. The important prediction is not that pollinators explain one C/S category, but that within-species colour discontinuities repeatedly align with independently measured pollinator turnover.

**Abiotic pigment benefit / spatial local adaptation.** Climate alone was too coarse. The new primary family therefore adds SoilGrids edaphic turnover alongside CHELSA. A positive result establishes environmental concordance, not adaptation; reciprocal fitness or genomic evidence is still required for that stronger claim.

**Gene-flow / geographic isolation.** Marine gaps, large rivers and mountain/biogeographic boundaries can now be tested against colour discontinuity. These are barrier proxies, not gene-flow measurements. A pattern dominated by barriers with weak present-day environmental alignment is compatible with historical isolation/drift but is not proof of it.

**Reproductive interference / character displacement (legacy H14).** This becomes a predeclared G5 context hypothesis. Sympatric colour divergence is predicted to be strongest among congeners, especially where flowering overlap is high. This is the cleanest new way to use the user's proposed coexistence/non-sympatry comparison.

### Retained as temporal or trait moderators

**Temporal environmental heterogeneity.** Static BIO4 was a weak/confounded candidate. The decisive test remains occurrence-year / flowering-season dynamic climate. It should be used later as a moderator of species-level spatial colour structure rather than as another static niche-width metric.

**Phenological partitioning.** It enters G5 as a secondary co-flowering qualifier and later into dynamic pollinator/climate analyses. Geographic sympatry alone is not called co-flowering coexistence.

**Autonomous reproductive assurance, clonality and seed-bank storage.** These remain plausible modifiers of whether colour variation persists locally versus spatially sorts, but require balanced trait coverage. SC/SI is still not a substitute for autonomous selfing.

### Kept outside atlas-level causal claims

**Multiple-agent balancing selection, pollinator–herbivore antagonism, negative frequency dependence and overdominance** require morph-specific visitation, fitness or genotype data. A global photograph atlas can prioritize candidate systems but cannot confirm these mechanisms.

**Plastic/nonheritable colour** remains an interpretation guard: the atlas establishes phenotypic colour geography unless heritability is independently known.

## 3. Competing G5 community hypotheses

The coexistence comparison should not assume in advance that sympatry produces one sign.

### G5a — shared-filter convergence

Sympatric species may converge in flower colour because they experience similar pollinator guilds or abiotic filters.

Prediction: `delta < 0`, potentially strongest among noncongeneric species with high pollinator overlap or very similar climate/soil environments.

### G5b — character-displacement divergence

Sympatric species may diverge in colour because discrimination, heterospecific pollen transfer or reproductive interference rewards distinct signals.

Prediction: `delta > 0`, potentially strongest in congeners with high flowering overlap and shared pollinator space.

### G5c — context-dependent sign

A near-zero global mean is not automatically a failure. Convergence and divergence may coexist in different contexts. Therefore the predeclared context strata are pollinator overlap, taxonomic distance, flowering overlap where independently estimable, and abiotic similarity. The sign is estimated rather than chosen after seeing colour.

## 4. What would constitute a strong result

The strongest possible sequence is:

1. G1 detects a recurrent colour field;
2. G2 shows stable `Z` zones;
3. G3 shows that within-species colour discontinuities track one or more external factors beyond geographic distance;
4. G4 shows that the same factor is elevated inside one or more neutral colour zones;
5. G5 independently shows non-random colour convergence or divergence among sympatric species relative to matched allopatric controls.

Even if G1 is null, a combination of heterogeneous species-specific G3 effects plus a supported G5 community-assembly pattern would still produce a coherent biological result: **flower-colour geography is not globally synchronized, but it is structured by local environmental/geographic context and by the identity of neighbouring flowering species.**

## 5. Claim boundaries

Do not translate correlation into `pollinator-driven`, `soil-adapted`, `vicariant`, `competitive displacement` or `mimicry` without direct evidence. The atlas can establish recurrent geography, concordance, community convergence/divergence and species-level heterogeneity. Mechanistic causal language remains reserved for experiments, genetics or direct fitness evidence.
