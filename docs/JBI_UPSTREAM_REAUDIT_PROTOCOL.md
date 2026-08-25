# JBI upstream re-audit protocol: natural floral-colour variation and spatial evidence states

## Status and purpose

This protocol starts **upstream of the historical 34-species dataset and upstream of automatic taxon extraction**. The historical freeze is never silently relabelled or overwritten. A replacement dataset is created only after literature retrieval, record relevance, focal taxon, natural eligibility, and spatial evidence have been independently reviewed before climatic responses are inspected.

The historical 34-species result remains a provenance-preserving sensitivity benchmark until a new adjudicated freeze exists.

## Why the re-audit is necessary

Audit identified several distinct selection/error mechanisms in the previous evidence path:

- generic `color/colour morph` terminology could be converted into a within-population label without direct same-population coexistence evidence;
- geographic evidence vocabulary was narrower than within-population vocabulary;
- mixed within/among candidates were explicitly removed before exploratory ecological modelling;
- legacy P1/P2/P3 screening priority affected entry into candidate analyses even though **15/34** exact historical source papers were outside P1;
- permissive two-word candidate extraction admitted non-taxa before taxonomic validation;
- conversely, requiring an automatically detected binomial can silently lose true source papers when a species is absent or abbreviated in title/abstract;
- cultivated, induced, ontogenetic, continuous, or non-display cases could be excluded before direct source adjudication;
- automated abstract-level signals were allowed to stand in for final biological evidence.

The replacement workflow therefore treats search-query membership, legacy priority, automatic taxon hints, historical labels, and automated within/among signals as **coordinator-only diagnostics**, never biological inclusion criteria.

## Core biological representation: two evidence axes

For every adjudicated eligible species, spatial organization is represented by two independent source-supported axes:

- `local_coexistence_documented`: at least one reviewed source explicitly documents two or more discrete natural floral-display colour variants coexisting within the same population;
- `geographic_structure_documented`: at least one reviewed source explicitly documents among-population, regional, clinal, island, elevational-spatial, or other geographic differentiation in floral colour or colour-morph frequency.

After all retained sources are adjudicated:

- local documented, geographic not documented → `within_evidence_only`;
- geographic documented, local not documented → `among_evidence_only`;
- both documented → `mixed_evidence`;
- insufficiently resolved evidence → `unresolved`.

`mixed_evidence` is not intermediate or ordinal. The two positive axes may be documented in the same source or in different sources for the same taxonomically resolved species.

A zero means **not documented in the reviewed evidence**, not biological absence.

## Stage 0. Replacement search-completeness gate — passed

The archived 79,242-record systematic map is retained unchanged for provenance, including its 19 truncated query/database shards. It is no longer the canonical completeness boundary.

The replacement retrieval is OpenAlex title/abstract OQL **v2.2**:

- 15 conceptual query blocks;
- no date restriction and no language restriction;
- Crossref used for metadata/identifier resolution rather than completeness-defining retrieval;
- non-English blocks broadened only after bounded scope probes;
- lexical omissions repaired generically, including singular/plural `color/colour morph(s)` and geography/geographic-structure wording;
- **13,911** query memberships;
- **12,064** deduplicated works;
- **0** truncated v2.2 query blocks;
- exact recovery of **34/34** historical classification sources;
- all seven prespecified review seeds resolved; benchmark recovery is 34/34 by the direct/citation audit.

No species-specific search term, climatic value, effect direction, or model result was used to tune v2.2.

Passing Stage 0 establishes the retrieval surface only. It does not make any of the 12,064 records biologically eligible.

## Stage 1. Canonical all-record blind screening — before taxon validation

**All 12,064 deduplicated v2.2 records enter record-level screening.**

No record can be excluded because an automatic binomial extractor fails. This correction is essential because the diagnostic automatic taxon universe recovered only 31/34 historical species despite v2.2 containing all 34 exact source papers.

Reviewer-facing material contains:

- record ID;
- title and abstract when available;
- year, journal, work type, language, source URL/identifier;
- blank coding fields.

Reviewer-facing material does **not** contain:

- query membership;
- legacy P1–P5 priority;
- historical-34 benchmark status;
- old within/among label;
- automatic binomial/taxon hints;
- automated natural/cultivated/within/among signals;
- climatic variables or manuscript model results.

For each record, two independent reviewers code:

1. `record_relevance`: `include` / `exclude` / `uncertain`;
2. `natural_intraspecific_variation`: `yes` / `no` / `uncertain`;
3. `floral_display_colour`: `yes` / `no` / `uncertain`;
4. `focal_taxon_text`: taxon text exactly as presented by the source;
5. `full_text_required`: `yes` / `no`;
6. a controlled exclusion reason when exclusion is clear;
7. notes where necessary.

If title/abstract evidence is incomplete, the conservative action is `uncertain` and/or `full_text_required = yes`, not exclusion.

The detailed coding rules are fixed in `docs/JBI_V22_RECORD_SCREENING_CODEBOOK.md`.

## Stage 1a. Wave 0 calibration gate

Before screening all 12,064 records, both reviewers independently code a **384-record blinded calibration sample**.

The hidden strata are:

- all **34** historical benchmark sources;
- **100** records without an automatically detected binomial;
- **100** records with an automatically detected binomial;
- **100** non-English records;
- **50** records with missing abstracts.

The strata are mutually exclusive in Wave 0 and their identities are hidden from reviewers. Review order is deterministic but independent of stratum.

Use separate reviewer copies. Reviewer 1 fills only `reviewer_1_*` fields; Reviewer 2 fills only `reviewer_2_*` fields. Neither reviewer sees the other's sheet before both copies are locked. The coordinator key is not distributed to either reviewer.

After duplicate coding, calculate raw agreement and Cohen's kappa for:

- record relevance;
- natural intraspecific variation;
- floral display colour;
- full-text requirement.

Also report normalized exact agreement for `focal_taxon_text`.

If disagreement reveals a material ambiguity in the codebook, revise the rule **before** full screening and run another calibration wave. Do not silently change rules during the 12,064-record screen.

The same person or AI agent cannot count as two independent reviewers.

## Stage 2. Record-level adjudication and full-text escalation

After independent coding is locked:

1. calculate agreement statistics;
2. list all field-level disagreements;
3. adjudicate disagreements by consensus or a named adjudicator;
4. retrieve full text where `full_text_required = yes` or where adjudication cannot be completed from bibliographic material;
5. preserve an adjudication/correction log;
6. retain `uncertain` records until the evidence gap is explicitly resolved.

Only after this step can a record be excluded from the canonical evidence universe.

## Stage 3. Focal taxon extraction and validation

Taxon validation occurs **after** record screening, using reviewer/adjudicator focal-taxonomy text plus source/full-text information.

Requirements:

1. preserve the taxon spelling/text originally used by the source;
2. resolve candidate names against a taxonomic authority, currently GBIF species matching plus accepted-name resolution;
3. require Plantae and species rank for species-level aggregation;
4. preserve rejected/ambiguous names in an audit table rather than dropping them silently;
5. resolve abbreviated or historical nomenclature using source context where needed;
6. merge synonym-derived evidence only after accepted-name resolution;
7. keep original input names and source identifiers alongside accepted names.

Automatic title/abstract binomial extraction may assist coordinators but cannot define inclusion.

## Stage 4. Source-level biological eligibility

A species-source record is eligible as biological evidence only when it supports naturally occurring intraspecific **floral-display colour** variation, or wild-origin material whose natural colour variation is independently documented.

Exclude when support is limited to:

- cultivars, horticultural lines, breeding collections, or commercial varieties without independent natural-population evidence;
- induced mutation, transgenics, gene editing, irradiation, mutagenesis, or tissue culture;
- ontogenetic colour change within a flower with no stable among-individual/population variation;
- interspecific differences only;
- non-floral colour only;
- non-display floral structures without a display-colour question;
- experimental colour treatments without corresponding natural variants.

Eligibility may remain `unclear` when the source is insufficient. Automated flags do not make the final decision.

## Stage 5. Source-level variation form and spatial evidence

### Variation form

Code natural variation as `discrete`, `continuous`, `both`, or `unclear`.

Continuous colour ITV remains scientifically relevant but is not silently converted into discrete polymorphism. The present local-coexistence state requires discrete variants.

### Local coexistence positive criterion

Code `local_coexistence_documented` only when the source directly establishes at least two discrete natural floral-display colour variants in the **same population**.

Sufficient evidence includes explicit co-occurrence statements, same-population counts containing multiple morphs, or morph-frequency tables with >1 morph in at least one population.

Generic polymorphism, pooled multi-site observations, different-locality samples, or common-garden coexistence without natural same-population evidence are insufficient alone.

### Geographic structure positive criterion

Code `geographic_structure_documented` only when the source directly demonstrates spatial differentiation among populations, sites, regions, islands, or geographic/elevational units.

Sufficient evidence includes morph-frequency differences among geographic units, spatial clines, regional replacement/restriction, or statistical population-level spatial differentiation in colour.

Sampling many locations, broad range, environmental association without spatially resolved colour contrast, or non-colour differentiation are insufficient alone.

## Stage 6. Species-level multiscale aggregation

After source adjudication and synonym resolution:

- aggregate positive local evidence across all retained sources;
- aggregate positive geographic evidence across all retained sources;
- assign within-only / among-only / mixed / unresolved from the two evidence axes.

Evidence for the two axes need not come from the same paper or same portion of the range. This is precisely why `mixed` is retained as a biological multiscale state rather than treated as a classification error.

## Stage 7. Observation-process audit

Because the outcomes are literature-derived documented-evidence states, define observation-process variables before climatic modelling, including where feasible:

- number of retained sources per species;
- years/temporal span of evidence;
- whether full-text evidence was available;
- geographic sampling breadth described by sources;
- occurrence/range support used for climatic calculations.

These variables are assessed as potential detection/documentation covariates. Documentation is not assumed to be perfect observation of biology.

## Stage 8. Ecological inference with mixed retained

Primary analysis fits parallel models for probability that local coexistence is documented and probability that geographic structuring is documented. A true mixed species contributes positive evidence to both outcomes.

If all three informative states have adequate sample size and stable estimation, a secondary **non-ordinal** multinomial analysis models within-only / among-only / mixed. Regularisation or weakly informative priors are used when required by sparse cells.

The historical binary within-vs-among analysis becomes a sensitivity analysis restricted to adjudicated non-mixed, non-unresolved cases. Family/phylogenetic dependence and finite-sample uncertainty must be addressed at least as rigorously as in the historical analysis.

## Stage 9. Replacement freeze and rerun gate

A new statistical freeze is created only after:

- v2.2 retrieval is fixed;
- all-record screening and required full-text escalation are adjudicated;
- focal taxa are resolved;
- source-level natural eligibility, variation form, and both spatial axes are independently reviewed and adjudicated;
- mixed cases are retained;
- reviewer agreement is reported;
- observation-process covariates are defined before outcome modelling;
- climatic metrics are recomputed for the adjudicated species set;
- all models, robustness checks, figures, manuscript text, SI, references/data-source appendix, and DOCX files are regenerated.

The historical 34-species checksum and outputs are never overwritten silently.

## Legacy-priority diagnostics retained only for sensitivity

The archived P1/P2/P3 workflow remains useful for diagnosing screening-path effects, but it is not the canonical inclusion route.

Current diagnostics:

- P1: **209** GBIF-valid candidate species and **41** automated mixed-navigation candidates;
- P2/P3: **522** GBIF-valid candidate species, all navigation-`unresolved` because the legacy P2/P3 screen did not carry directional within/among signals;
- P1–P3 union: **656** candidate species, **447** more than P1 alone, with **41** mixed-navigation cases because only legacy directional signals can generate that automated state.

These counts demonstrate why legacy priority and automated state signals cannot define the replacement biological dataset.
