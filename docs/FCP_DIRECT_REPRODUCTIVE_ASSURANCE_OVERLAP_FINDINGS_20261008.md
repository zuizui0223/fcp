# FCP white–chromatic flower colour and reproductive assurance: direct evidence ceiling
## Source-replayed 2026-10-08 findings

### Completed result

The global FCP source has 149,900 original candidate photographs of 1,499 sampled nominal species, three species-disjoint cohorts. The same originally frozen 40+ classifiable-photograph gate yields **1,109 eligible photo species: discovery 369, validation 363, third 377**.

The new source-controlled analysis tested exact scientific-binomial overlap with **three independent original organism-level reproductive datasets**, without replacing experimental measures by self-compatibility or predicted flower syndrome.

| Original reproductive study | Distinct FCP matched species with measured trait | Discovery | Validation | Third | Matched species with >=5 white and >=5 nonwhite photographed flowers |
|---|---:|---:|---:|---:|---:|
| Goodwillie, Kalisz and Eckert (2005) natural-population multilocus outcrossing | **27** | 7 | 13 | 7 | **15** (4 / 7 / 4) |
| Razanajatovo et al. (2016) experimental autofertility index | **38** | 16 | 12 | 10 | **23** (10 / 5 / 8) |
| Rodger et al. (2021) reproductive fruit/seed output with pollinators excluded | **44** | 14 | 16 | 14 | **21** (9 / 5 / 7) |

These are *three overlapping source panels*, not 109 independent plants. One FCP species can occur in more than one study. Both source original zero reproductive outputs and natural-population genetic outcrossing values were retained appropriately; self-compatibility-only studies were NOT counted as autonomously reproductive plants.

Original source/public DOI:
- Goodwillie, Kalisz & Eckert 2005: 10.1146/annurev.ecolsys.36.091704.175539
- Razanajatovo et al. 2016: 10.1038/ncomms13313
- Rodger et al. 2021: 10.1126/sciadv.abd3524; supplemental direct exclusion dataset 10.6084/m9.figshare.14607882.v1.

The selection of a taxon name mapping for Rodger is traceable to the source implementation in the external island repository: raw genus.species (with underscores replaced by spaces), not automatically cherry-picked among several alternate name fields. Missing, invalid negative and nonnumerical raw reproductive outputs were excluded; zero is a genuine measured result.

### Prespecified cross-cohort feasibility gate: HOLD

A separate new observational FCP-scale cross-species correlation would require, **in one directly measured trait source**:
- >=80 distinct FCP species with that same direct reproductive variable;
- >=20 in each of the three independent photo cohorts;
- >=8 photographed white+chromatic species in each photo cohort.

**NONE** of the direct reproductive panels meets these gates. Goodwillie n=27, Razanajatovo n=38, Rodger n=44. Even the broadest original Rodger panel has just 5 white+chromatic species in the validation photo cohort. Therefore no trait–photo white-frequency association or mechanistic trade-off model is fit in this PR. This is an evidence-coverage and sample-comparability failure, **NOT a negative result for reproduction-based flower-colour selection in nature**.

The eligible source frame of 1,109 is itself biased toward observation-rich species and only measures photographic white/nonwhite states. A source-level species' outcrossing or autonomous seed set is not a morph-specific selection coefficient. The actual white–coloured trade-off would require directly comparing white and coloured individuals of the SAME natural population, under contrasting pollinator pressure, with seeds/survival to estimate their relative fitness.

### What is now known biologically

- Existing FCP photos show within-species four-colour variation geographically clustered even after month/year and observer control, in three source-disjoint photo cohorts. It does NOT identify whether the spatial sorting is maintained by drift, migration, local selection or floral plasticity.
- Literature-confirmed ecological mechanisms include pollinator-mediated rare-colour advantage in Dactylorhiza sambucina (yellow/purple, not biochemical white), local flower-number/pollinator/autonomous selfing trade-offs in Silene littorea, whole-plant vs petal-specific pigment loss in Silene, water/herbivore-dependence in Boechera stricta, and drift/gene flow in Iris lutescens.
- A preceding 11-organism primary source-to-photo audit had only three shared taxa; zero exact source-white-genotype-to-individual-photo correspondence. Enlarging to global *direct* reproductive datasets improved trait candidate coverage to 27/38/44 species but still does not bridge photos to genetically proven competing white/pigment fitness payoffs.
- Consequently the only broad supported evolutionary/ecological statement remains repeated **spatial structure of photographed flower-colour diversity**, with *multiple plausible nonexclusive processes*. An attractive universal pigment maintenance benefit–cost mechanism remains unverified.

### Next biological leverage, conditional on available data

The most informative next evidence is paired colour and *true reproductive assurance measured on the SAME species, population and genetic morph*, or independent biochemistry of petal-specific white variants. This discriminates selection on floral pigment from selection on linked flower production / whole-plant defenses / mating traits.

Broader taxonomic trait collections with source-typed or machine-inferred selfing proxies could raise sample counts, but are **not equivalent** to measured autonomous seed/fruit set. Their use would require a distinct, provenance-explicit study, and they cannot repair the missing morph-specific fitness identity simply by having more rows. No current manuscript H1/H2 or prior geography claim is altered.

### Exact execution receipts

FCP draft PR #132: https://github.com/zuizui0223/fcp/pull/132

Successful source-replay workflow: https://github.com/zuizui0223/fcp/actions/runs/37754332413, job 113235033281. Original FCP measured tables passed all three frozen SHA256 checks. Public independent island source pinned at exact commit 92f007797fbbcb3c6743a4ffb2628fa609329499, with Rodger compressed original projected source SHA256 fa745c578f3537933fafedc1d36b4ea266348cd83d7f6cbb231c253b0f348d3f.

6 synthetic failure-mode checks passed. Output artifact ID 11539227632 digest sha256 bf58647d35cdc9ef2019e100cbb527c13b847ff170e6c954cef47a550c382d2d. Result JSON: results/polymorphism_direct_reproductive_assurance_source_overlap_posthoc_20261008/result.json. Per-species source-matched opportunity rows in workflow artifact file matched_species_direct_trait_opportunity.csv.

### Frozen claim boundary

Do not infer a globally shared balancing selection trade-off from these sample counts, sum different selfing endpoints into one quantitative coefficient, or treat presence of an observed white flower as a validated white anthocyanin allele. HOLD pertains only to the current materialized independent source–photo observational linkage and frozen coverage rule, not the existence of costs and benefits in natural populations.
