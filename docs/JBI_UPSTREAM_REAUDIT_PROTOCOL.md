# JBI upstream re-audit protocol: natural floral-colour variation and spatial evidence states

## Status and purpose

This protocol starts **upstream of the historical 34-species dataset**. It does not edit, relabel, or overwrite the historical freeze. Its purpose is to construct a new source-reviewed evidence layer before any climatic response is inspected.

The historical 34-species result remains a provenance-preserving sensitivity benchmark until a new adjudicated freeze exists.

## Why the re-audit is necessary

The repository contained two different discovery paths:

1. historical manuscript chain: 1,075 retained works → 664 candidate species → 72 initial review rows → 111 resolved rows → 34 frozen binary species;
2. later systematic-map chain: 105,249 raw records → 79,242 deduplicated records.

Audit showed that neither path was an adequate final biological classification. Important weaknesses included:

- generic morph terminology could trigger a within-population label without direct evidence of coexistence in the same population;
- geographic evidence vocabularies were narrower than within-population vocabularies;
- provisional two-word species extraction admitted non-taxa before taxonomic validation;
- automated screening priority could affect entry into the ecological candidate set;
- cultivated/induced/ontogenetic conflicts could be discarded too early;
- species carrying both within- and among-population signals were explicitly removed as `mixed` before ecological modelling;
- automated abstract-level signals were allowed to stand in for final source-level adjudication.

The re-audit therefore treats all automated labels and priorities as navigation/diagnostic variables only.

## Core conceptual correction: two evidence axes, not a forced binary state

For every eligible species, spatial organization is recorded on two independent source-supported axes:

- `local_coexistence_documented`: at least one reviewed source explicitly documents two or more discrete natural floral-display colour variants coexisting within the same population;
- `geographic_structure_documented`: at least one reviewed source explicitly documents among-population, regional, clinal, elevational-spatial, or other geographic differentiation in floral colour or colour-morph frequency.

After all retained sources for a species are reviewed, these axes generate three informative documented states:

- `within_evidence_only`: local coexistence documented, geographic structuring not documented in the reviewed evidence;
- `among_evidence_only`: geographic structuring documented, local coexistence not documented in the reviewed evidence;
- `mixed_evidence`: both are documented.

`unresolved` is retained when eligibility or either spatial axis cannot be adjudicated adequately.

The words **not documented** are intentional. A zero on an evidence axis means that the reviewed source set contains no adjudicated positive documentation for that axis; it is not interpreted as proof that the biological process is absent.

`mixed_evidence` is not an intermediate or ordinal category. Local coexistence and geographic structure may be documented in the same paper or in different independent sources for the same taxonomically resolved species.

## Stage 0. Replacement search-completeness gate — passed

The legacy 79,242-record search remains archived unchanged for provenance, including its 19 truncated query/database shards (15 Crossref, 4 OpenAlex). Those truncations are no longer repaired by increasing legacy caps.

A replacement search surface, OpenAlex title/abstract OQL v2.2, is now the completeness-defining retrieval layer:

- 15 prespecified conceptual query blocks;
- title/abstract search rather than ordinary OpenAlex search that may include full text;
- Crossref retained for bibliographic metadata and DOI resolution, not as the completeness-defining retrieval source;
- non-English blocks expanded only after bounded scope probes;
- generic lexical omissions repaired without species-specific query terms;
- 13,911 query memberships;
- 12,064 deduplicated works;
- zero truncated v2.2 query blocks;
- 34/34 historical classification sources recovered by exact identifier in direct v2.2 queries;
- all seven prespecified review seeds resolved and 34/34 historical benchmark sources recovered by direct or prespecified citation pathways.

Search terms were repaired from search-surface diagnostics and benchmark lexical failures only. Climatic values, effect directions, and model outcomes were not used to tune the search.

Passing the search gate does **not** mean that all 12,064 works are biological inclusions.

## Stage 1. Priority-independent source-to-taxon universe

The canonical candidate universe is built from **all 12,064 deduplicated v2.2 works**, not from legacy P1/P2/P3/P4/P5 priority strata.

For each retrieved work:

1. extract candidate binomials from the title; if none occur there, inspect the first 2,000 abstract characters for focal binomials;
2. limit source attribution to a bounded number of focal names so review papers do not become artificial species censuses;
3. validate attributed names against a taxonomic authority before species-level aggregation;
4. keep query membership and all legacy priorities in coordinator-only material;
5. provide reviewers only bibliographic/source material and blank coding fields.

Legacy P1/P2/P3 analyses are retained as **screening-path sensitivity diagnostics**, not as the canonical inclusion path.

## Stage 2. Taxon validation before spatial classification

Candidate names are validated before evidence is aggregated into biological states.

Current requirements:

1. strict GBIF species matching in Plantae;
2. accepted-name resolution for synonyms;
3. species rank and resolved family required;
4. rejected candidate phrases retained in an audit table;
5. all original input names and stable source identifiers preserved;
6. evidence from synonyms merged only after accepted-name resolution.

No species is assigned `within`, `among`, or `mixed` from the taxonomic step.

## Stage 3. Source-level natural-variation eligibility

A species-source record is eligible as biological evidence only when the source supports naturally occurring intraspecific floral-display colour variation, or wild-origin experimental material whose natural variation is independently documented.

Exclude as biological evidence when support is limited to:

- cultivars, horticultural lines, breeding collections, or commercial varieties without independent natural-population evidence;
- induced mutations, transgenic manipulation, gene editing, irradiation, mutagenesis, or tissue culture;
- ontogenetic colour change within an individual flower rather than stable variation among individuals/populations;
- non-display floral organs only;
- interspecific variation only;
- experimental colour treatments without naturally occurring corresponding variants.

A source can be coded `eligible`, `ineligible`, or `unclear`. Automated natural/cultivated flags do not determine this decision.

## Stage 4. Source-level coding of variation form and spatial evidence

### 4.1 Variation form

Record whether natural variation is:

- `discrete`;
- `continuous`;
- `both`;
- `unclear`.

The current spatial-state analysis requires discrete floral-colour variants for `local_coexistence_documented`. Continuous variation remains part of the broader ITV evidence map but is not silently converted into discrete polymorphism.

### 4.2 Local coexistence: positive evidence criterion

Code positive local coexistence only when the source directly supports at least two discrete natural floral-display colour variants **within the same population**.

Sufficient examples include:

- explicit statements such as “co-occur/coexist within populations”;
- population-level counts showing multiple colour morphs in the same named population/site;
- morph-frequency tables where more than one morph occurs within at least one population.

Insufficient by itself:

- “multiple colour morphs exist” with no population resolution;
- pooled observations from several sites;
- two morphs collected from different localities;
- common-garden coexistence without independent natural same-population evidence;
- a species described generally as polymorphic.

### 4.3 Geographic structure: positive evidence criterion

Code positive geographic structure only when the source directly supports spatial differentiation among populations, regions, sites, islands, elevations represented by distinct populations, or other geographic units.

Sufficient examples include:

- morph frequencies differing among named populations/sites;
- geographic or latitudinal/elevational clines in colour or morph frequency;
- regional replacement/restriction of colour forms;
- statistical spatial differentiation in colour phenotype among populations.

Insufficient by itself:

- sampling many localities without reporting colour differentiation among them;
- an environmental association with no spatially resolved colour contrast;
- a statement that a species has a broad geographic range;
- population differentiation in a non-colour trait only.

### 4.4 Mixed/multiscale aggregation

After source-level adjudication, a species is `mixed_evidence` whenever both axes have at least one positive source-supported observation after synonym resolution.

The two positive observations:

- need not come from the same paper;
- need not occur in the same part of the range;
- may describe local coexistence in some populations and geographic restriction or frequency gradients across the wider range.

This rule is essential because multiscale spatial organization is the biological object of interest rather than a classification error to be discarded.

## Stage 5. Duplicate blind review and adjudication

For every species-source unit considered for the replacement freeze:

1. reviewer 1 codes display-colour relevance, natural eligibility, variation form, local coexistence, geographic structure, and evidence passage;
2. reviewer 2 independently codes the same fields;
3. reviewers do not receive climatic values, historical binary labels, automated within/among signals, query membership, or coordinator priority reasons;
4. disagreements are resolved by consensus or a named adjudicator after independent coding is locked;
5. raw agreement and Cohen's kappa are reported separately for eligibility and both spatial axes;
6. all adjudications are written to a correction/adjudication log;
7. only adjudicated species enter the new statistical freeze.

Reviewer-facing values for a spatial axis should distinguish `documented`, `not_documented_in_source`, and `unclear`; `not_documented_in_source` must never be described as evidence of biological absence.

The same agent/person should not be counted as two independent reviewers.

## Stage 6. Review ordering without inclusion bias

Coordinator-only triage may determine review order, never inclusion.

Priority A:
- historical 34 sources;
- species flagged by automated diagnostics as potentially multiscale/mixed.

Priority B:
- candidates supported by multiple independent sources.

Priority C:
- eligibility conflicts, unresolved spatial evidence, cultivated/induced/ontogenetic conflicts.

All remaining v2.2 source-taxonomy units stay in the review universe. Legacy P1/P2/P3 priority membership is a sensitivity descriptor only.

## Stage 7. Ecological analysis with mixed retained

### Primary analysis: two parallel documented-evidence outcomes

For each climatic metric, fit parallel models for:

- probability that local coexistence is documented;
- probability that geographic structuring is documented.

A mixed species contributes positive evidence to both outcomes.

Because these are literature-derived documented-evidence outcomes, models must evaluate/report literature effort and source availability as potential observation-process covariates rather than silently treating documentation as perfect biological detection.

Family/phylogenetic dependence and finite-sample uncertainty should be addressed at least as rigorously as in the historical analysis.

### Secondary analysis: three non-ordinal documented states

If the adjudicated sample contains sufficient information in all three states and estimation is stable, fit a multinomial model for:

`within_evidence_only` / `among_evidence_only` / `mixed_evidence`.

Use regularisation or weakly informative priors if sparse cells make ordinary maximum-likelihood estimates unstable. Do not treat the states as ordinal.

### Legacy binary analysis

The historical `within` versus `among` model becomes a sensitivity analysis restricted to adjudicated evidence-only cases after excluding mixed and unresolved species.

## Stage 8. Replacement freeze and rerun gate

A replacement statistical freeze is created only after:

- v2.2 retrieval provenance is fixed;
- candidate taxa are validated;
- source-level natural eligibility and both spatial axes are independently reviewed and adjudicated;
- mixed cases are retained;
- reviewer agreement is documented;
- observation-process/literature-effort variables needed for inference are defined before outcome modelling;
- climatic metrics are recomputed for the adjudicated species set;
- the complete model, robustness, figure, manuscript, SI, reference/data-source appendix, and DOCX pipeline is rerun.

The historical 34-species checksum and outputs are never silently overwritten.

## Historical and automated diagnostics retained for sensitivity only

The archived 2026-08-02 systematic search and exploratory ecological artifacts remain useful for diagnosing screening-path effects. They showed, among other things, that mixed candidates had been dropped and that legacy P1 priority recovered only 22/34 historical species even though all 34 historical source papers existed in the archived corpus.

These diagnostics motivate the new workflow but do not adjudicate any species.
