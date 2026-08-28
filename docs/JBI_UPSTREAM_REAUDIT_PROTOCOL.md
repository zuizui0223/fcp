# JBI upstream re-audit protocol: local coexistence and spatial segregation of flower-colour polymorphism

## Status and purpose

This protocol replaces the historical binary `within_population` / `among_population` inclusion path for the re-audit. It starts upstream of the historical 34-species freeze and preserves that old freeze unchanged as a provenance/sensitivity analysis.

The canonical replacement is a **single-pass, conservative documented-evidence audit**. Independent duplicate screening infrastructure is retained in the repository but was explicitly waived on 2026-08-26 and is not a canonical gate. The replacement analysis therefore must not be described as an independently double-reviewed systematic review.

Canonical path:

**OpenAlex v2.2 title/abstract retrieval → all-record evidence universe → strict positive-evidence screening → GBIF taxon resolution → C/S source aggregation → versioned C/S evidence freeze → fresh GBIF occurrence reconstruction → WorldClim 2.1 climate metrics → primary C/S models and robustness analyses.**

Climate values, model effects, historical labels, and old manuscript conclusions are never used to assign C or S.

## Why the old result required upstream re-audit

The previous evidence path could:

- convert generic `color/colour morph` terminology into within-population coexistence without direct same-population evidence;
- detect geographic structure with a narrower vocabulary than local polymorphism;
- remove mixed within+among cases before modelling;
- let legacy screening priority affect entry into ecological analyses;
- admit false two-word taxon candidates before taxonomic validation;
- miss true papers when automatic binomial extraction failed;
- retain cultivated, induced, ontogenetic, continuous-only or non-display cases too readily;
- propagate a deterministic upstream label through freeze, models, figures and manuscript without an independent validity check.

The re-audit treats legacy priority, old labels, automatic signals and climatic outcomes as diagnostics only.

## Stage 0 — fixed retrieval surface

The completeness-defining replacement search is OpenAlex title/abstract OQL v2.2:

- 15 conceptual query blocks;
- 13,911 query memberships;
- 12,064 deduplicated works;
- zero truncated v2.2 query blocks;
- no date restriction and no language restriction;
- exact recovery of all 34 historical classification sources;
- all seven prespecified review seeds recovered through the direct/citation benchmark.

The archived 79,242-record systematic map and its 19 truncated legacy shards remain provenance material, not the canonical completeness boundary.

All 12,064 v2.2 records remain visible in the audit universe. Failure of automatic taxon extraction is never itself an exclusion criterion.

## Stage 1 — biological representation: two independent positive-evidence axes

The central variables are **not complements**.

### C — local coexistence

`C_local_coexistence_documented = 1` only when a retained source explicitly supports at least two discrete, naturally occurring floral-display colour variants in the same population or site.

Sufficient positive evidence includes:

- explicit co-occurrence/coexistence within one population/site;
- a population explicitly described as flower-colour polymorphic when the colour context is discrete and natural;
- same-population counts or morph-frequency information showing more than one colour morph;
- mixed populations or equivalent direct population-resolved wording.

Insufficient by itself:

- `color/colour morph(s)`;
- species-level `polymorphic` or `polymorphism` without population resolution;
- multiple colour forms pooled across locations;
- common-garden coexistence without natural same-population evidence;
- continuous colour variation alone.

### S — spatial segregation

`S_spatial_segregation_documented = 1` only when a retained source explicitly supports geographic differentiation in colour morph presence or frequency among populations, sites, localities, regions, islands, elevations or other spatial units.

Sufficient positive evidence includes:

- morph-frequency differences among geographic units;
- monomorphic versus polymorphic populations;
- geographic/clinal/elevational colour differentiation;
- spatial replacement or restriction of colour forms;
- explicit population-level contrasts in colour composition.

Insufficient by itself:

- broad range size;
- sampling many locations without a reported colour contrast;
- an environmental association without spatially resolved colour differentiation;
- generic colour polymorphism.

### Zero semantics

`C=0` means **local coexistence was not documented by this strict pass**; it does not mean coexistence is biologically absent.

`S=0` means **spatial segregation was not documented by this strict pass**; it does not mean spatial homogeneity is biologically established.

After species-level aggregation:

- C=1, S=0 → `local_coexistence_only`;
- C=0, S=1 → `spatial_segregation_only`;
- C=1, S=1 → `coexistence_and_segregation`;
- C=0, S=0 or insufficient evidence → `unresolved`.

`coexistence_and_segregation` is a multiscale biological state, not an intermediate category or classification error. C and S evidence may come from different papers.

## Stage 2 — strict eligibility and false-positive controls

A source may contribute positive evidence only for naturally occurring intraspecific **floral-display colour** variation.

The strict pipeline guards against:

- cultivar/horticultural-line/commercial-variety-only studies;
- artificial or experimental flower arrays;
- induced mutation, transgenics, gene editing, irradiation, mutagenesis or tissue culture;
- ontogenetic or post-pollination colour change without stable among-individual variation;
- interspecific/community-level colour comparisons;
- review/perspective text used as if it were a primary natural population observation;
- continuous-only colour measurements converted into discrete polymorphism;
- negated statements such as morphs that explicitly do **not** coexist;
- non-display sexual-organ colour only, including anther, pollen or androecium colour.

Ambiguous cases remain unresolved rather than being forced into C or S.

Detailed implementation is in the versioned refined builders under `scripts/literature/build_v22_coexistence_segregation_refined*.py`; the final workflow runs the v5 chain.

## Stage 3 — taxonomic resolution

Only C/S-positive or historical-rescue sources require taxon assignment for the strict freeze.

- candidate binomials are extracted from source title/abstract text;
- candidates are resolved against the GBIF backbone;
- a positive source contributes only when it resolves unambiguously to one accepted plant species;
- synonym-derived evidence is merged only after accepted-name resolution;
- unresolved positive sources stay in an explicit unresolved table and do not enter the freeze.

The old 34-species manifest is permitted only as an exact **source → taxon rescue map** for those 34 historical source records. Old spatial labels, old climate-cell counts and old model effects are not passed to the C/S evidence classifier.

## Stage 4 — fixed C/S evidence freeze

The final strict v5 evidence workflow is `JBI v2.2 refined coexistence-segregation audit`, successful run `32926071541`.

The versioned replacement evidence freeze contains:

- **34 species**;
- **11** `local_coexistence_only`;
- **8** `spatial_segregation_only`;
- **15** `coexistence_and_segregation`;
- **26** C-positive species;
- **23** S-positive species;
- 50 positive source rows, with one taxonomically unresolved positive source excluded from the freeze.

Durable files:

- `data/frozen/jbi_cs_evidence_freeze_v22.csv`;
- `data/frozen/jbi_cs_evidence_freeze_v22_manifest.json`;
- `docs/supporting/jbi_cs_positive_source_index_v22.csv`.

Evidence-freeze SHA-256: `be03ea91e7e1a1aa9e9dc8850be7d0b24fbf2f488c440086063e854de1db5836`.

The historical binary 34-species freeze is not overwritten.

## Stage 5 — fresh occurrence and climate reconstruction

Because only 10 of the 34 replacement species overlap the historical binary set, old climatic values are not reused.

The replacement climate workflow reconstructs occurrences from GBIF under the historical sampling conditions:

- maximum 3,000 occurrence records per species;
- page size 300;
- coordinates required and `hasGeospatialIssue=false`;
- accepted occurrence basis types restricted to observations/specimens/material samples used by the historical workflow;
- coordinate uncertainty ≤20 km;
- coordinate deduplication at three decimal places.

Climate extraction uses **WorldClim 2.1, 10 arc-minute** BIO1, BIO4, BIO5, BIO6, BIO7, BIO12, BIO14, BIO15 and BIO17. Identical nine-variable climate vectors are deduplicated as occupied climate cells. At least 20 occupied climate cells are required.

Successful rebuild run: `32926741272`.

Results:

- 47,942 retained GBIF coordinate rows;
- 15,559 occupied climate-cell rows after raster extraction/deduplication;
- all 34 species retained;
- zero species excluded by the ≥20-cell gate.

The five symmetric species-level metrics remain:

1. temperature breadth;
2. moisture breadth;
3. climatic heterogeneity;
4. PCA dispersion;
5. PCA hull area.

These describe realised occupied climate, not physiological tolerance or morph-specific niches.

Durable climate-analysis dataset: `data/frozen/jbi_cs_climate_analysis_v22.csv`.

Canonical climate-analysis SHA-256: `161adbe80ee3b38a60b17cd0ad1e048eb9d454ae1cfea5825192995ef39a9a42`.

## Stage 6 — primary ecological inference

C and S are modelled separately on the same species set:

- `C ~ metric_z + effort_z`;
- `S ~ metric_z + effort_z`.

`effort_z = z(log1p(n_climate_cells))`. Each of the five climatic metrics is fitted separately.

Primary uncertainty and robustness:

- family-clustered sandwich covariance;
- 9,999 response-label permutations;
- leave-one-family-out refits;
- Holm multiplicity context across the five metrics within each outcome.

A C+S species is positive in both models. Neither model treats the other axis as its negative class in a biological sense; zero is a documented-evidence state.

The clean canonical runner is `scripts/run_cs_models_v2.py` and explicitly excludes outcome-path-derived source counts as effort covariates.

### Primary result boundary

Across the five climate metrics:

- C point estimates are all OR<1, but permutation p-values range from about 0.37 to 0.70 and Holm-adjusted permutation p-values are all 1.0;
- S point estimates are near null, with permutation p-values about 0.71–0.97 and Holm-adjusted permutation p-values all 1.0;
- no clustered 95% interval excludes OR=1.

The replacement analysis therefore does **not** reproduce the historical moisture-breadth signal. It also does not establish biological absence of a relationship: the sample is small and the outcomes are conservative documented-evidence states.

## Stage 7 — observation-process and multistate sensitivities

`n_resolved_sources` is outcome-path-derived and is **not** treated as an independent documentation-effort covariate.

Observation-process sensitivity instead uses an outcome-independent literature-attention measure: exact accepted canonical binomial mentions in all 12,064 v2.2 title/abstract records. This measure is conservative for synonym-heavy taxa and is used only as sensitivity analysis.

Additional robustness analyses:

- literature-attention-adjusted C/S GLMs;
- a 23-species stratum restricted to exactly one resolved positive-evidence source, used to equalise positive-source count rather than as an effort adjustment;
- a secondary non-ordinal multinomial comparison of `local_coexistence_only`, `spatial_segregation_only` and `coexistence_and_segregation`.

All three sensitivity families remain compatible with a weak/null climate contrast. No multinomial climate comparison has a 95% interval excluding relative OR=1.

## Stage 8 — relation to the historical result

The old and replacement 34-species sets overlap by only **10 species**. Of those overlapping species, several change representation when multiscale evidence is retained, including:

- *Anemone palmata*: historical within → C+S;
- *Hesperis matronalis*: historical within → C+S;
- *Iris lutescens*: historical among → C+S;
- *Rhododendron arboreum*: historical within → S-only.

The remaining 24 replacement species are different from the historical binary set. Therefore the disappearance of the historical moisture contrast reflects a substantive upstream change in evidence validity, sample composition and state representation, not merely a change of statistical estimator.

The historical binary analysis remains reproducible and should be reported as a provenance-preserving sensitivity analysis, not silently rewritten.

## Archived duplicate-review infrastructure

The repository retains the 384-record Wave 0 material, Reviewer 1/2 packages, kappa gate code and B01–B13 packages. These are archived optional sensitivity/provenance tools. They are not required by the current single-pass canonical path, and no independent reviewer-agreement statistic is claimed.

## Reporting rules

- Describe C and S as **documented positive-evidence axes**, not demonstrated presence/absence of biological states.
- Use association language, not causation.
- Do not equate species-level occupied climate breadth with morph-specific tolerance.
- Do not call C/S states mutually exclusive; C+S is explicitly retained.
- Do not present the replacement audit as independently double-reviewed.
- Do not interpret a null C/S climate model as proof that climate never affects polymorphism organization.
- Keep the historical binary freeze and replacement C/S freeze as separate versioned analyses.
