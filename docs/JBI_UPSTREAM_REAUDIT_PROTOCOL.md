# JBI upstream re-audit protocol: natural floral-colour variation and spatial evidence states

## Status and purpose

This protocol starts **upstream of the frozen 34-species dataset**. It does not edit, relabel, or overwrite the historical 34-species freeze. Its purpose is to construct a new source-reviewed evidence layer from the broad systematic-search corpus before any climatic response is inspected.

The historical 34-species result remains a provenance-preserving sensitivity benchmark until a new adjudicated freeze exists.

## Why the re-audit is necessary

The repository currently contains two distinct discovery chains:

1. the historical chain used by the manuscript: 1,075 retained works -> 664 candidate species -> 72 initial review rows -> 111 resolved rows -> 34 frozen binary species;
2. the later systematic-map search: 105,249 raw records -> 79,242 deduplicated records.

The later search is the more appropriate upstream starting point, but its exploratory ecological trial was not a human-reviewed replacement for the 34-species freeze. Several implementation details make that trial unsuitable as a final biological classification:

- 19 search shards were reported as truncated; the systematic-search protocol requires partitioning or rerunning truncated searches before the search is declared final;
- provisional species extraction used a permissive two-word binomial regular expression, so non-taxa such as ordinary capitalized two-word phrases could enter the candidate layer;
- the provisional classification script accepted a `--census` argument but did not use it;
- only high-priority P1 records were used to create the exploratory ecological candidate set;
- candidate records carrying cultivated/induced/ontogenetic signals could be automatically discarded rather than retained for source review;
- species supported by both within- and among-population signals were explicitly excluded as `mixed` before ecological modelling;
- automated abstract-level signals were allowed to define the exploratory state, despite the protocol specifying that final natural status and spatial scale require source-level human review.

The re-audit therefore treats all automated labels as navigation signals only.

## Core conceptual correction: two evidence axes, not a forced binary state

For every eligible species, spatial organization is recorded on two independent source-supported axes:

- `local_coexistence_documented`: at least one reviewed source explicitly documents two or more discrete natural floral-display colour variants coexisting within the same population;
- `geographic_structure_documented`: at least one reviewed source explicitly documents among-population, regional, clinal, or other geographic differentiation in floral colour or colour-morph frequency.

These axes generate three informative documented states:

- `within_evidence_only`: local coexistence documented, geographic structuring not documented in the retained evidence;
- `among_evidence_only`: geographic structuring documented, local coexistence not documented in the retained evidence;
- `mixed_evidence`: both are documented.

`unresolved` is retained when neither axis is adequately resolved.

The words **not documented** are intentional. Absence of a published signal is not interpreted as biological absence.

## Stage 0. Search-completeness gate

Before a replacement dataset can be called final:

1. rerun or partition all truncated database/query shards;
2. preserve the original 2026-08-02 corpus and manifest unchanged;
3. merge the repaired retrieval into a new versioned corpus;
4. report database/query contribution, duplicate removal, truncation status, and benchmark recovery;
5. do not use climatic values or the historical model direction to tune search terms.

Until this gate is passed, the existing 79,242-record corpus is a high-value provisional discovery source, not a completed systematic map.

## Stage 1. Record prioritisation without biological exclusion

Automated priorities may determine review order but not final eligibility.

Initial source-level review should begin with P1 and P2 records and include P3 exclusion-review records. Lower-priority records are retained for false-negative auditing and citation chaining. Records with both natural and cultivated/induced signals remain reviewable rather than being silently removed.

For each record, retain separately:

- candidate taxon name(s);
- natural/wild-field signal;
- cultivated/horticultural signal;
- induced/transgenic/mutagenesis signal;
- ontogenetic colour-change signal;
- discrete vs continuous variation signal;
- local-coexistence signal;
- geographic-structure signal;
- display-organ signal;
- DOI/stable identifier and evidence excerpt.

## Stage 2. Taxon validation before species-level spatial classification

Candidate names are validated before evidence is aggregated into biological states.

Requirements:

1. resolve candidate names against a taxonomic authority (currently strict GBIF species matching, with accepted-name resolution for synonyms);
2. require Plantae, species rank, and a resolved family;
3. retain rejected candidate phrases in an audit table;
4. merge synonym-derived evidence only after accepted-name resolution;
5. preserve all original input names and source identifiers.

No species is assigned `within`, `among`, or `mixed` before this step.

## Stage 3. Natural-variation eligibility review

A species-study record is eligible only when the source supports naturally occurring intraspecific floral-display colour variation, or wild-origin experimental material whose natural variation is independently documented.

Exclude as biological evidence when support is limited to:

- cultivars, horticultural lines, breeding collections, or commercial varieties without independent natural-population evidence;
- induced mutations, transgenic manipulation, gene editing, or tissue culture;
- ontogenetic change within a flower rather than among individuals/populations;
- non-display floral organs only;
- interspecific differences only.

The review records the decision and reason; automated flags do not make the final decision.

## Stage 4. Duplicate source-level review and adjudication

For every species proposed for ecological analysis:

1. reviewer 1 codes natural eligibility, variation form, local coexistence, geographic structuring, and evidence passage;
2. reviewer 2 independently codes the same fields without climatic results or model labels;
3. disagreements are resolved by consensus or a named adjudicator;
4. raw agreement and Cohen's kappa are reported for eligibility and both spatial axes;
5. all accepted changes are written to a correction/adjudication log;
6. only adjudicated species enter a new statistical freeze.

## Stage 5. Ecological analysis with `mixed` retained

### Primary analysis: two orthogonal outcomes

The preferred primary analysis does not force mixed systems into one side of a binary contrast.

For each climatic metric, fit parallel models for:

- probability that local coexistence is documented;
- probability that geographic structuring is documented.

A mixed species contributes positive evidence to both outcomes. Models should retain literature-effort and occurrence/range covariates where justified, and family/phylogenetic dependence should be addressed as in the current paper or with a validated replacement.

This formulation directly tests whether climatic breadth is associated differently with local coexistence and geographic differentiation.

### Secondary analysis: three documented states

If the adjudicated sample contains enough information in all three states and the model is stable, fit a multinomial model for:

`within_evidence_only` / `among_evidence_only` / `mixed_evidence`.

Use regularisation or weakly informative Bayesian priors if sparse cells make ordinary maximum-likelihood multinomial estimates unstable. Do not treat the three states as ordinal.

### Legacy binary analysis

The historical `within` vs `among` model should become a sensitivity analysis restricted to evidence-only cases after excluding mixed and unresolved species. It is useful for continuity with the existing manuscript, but it should not define the new primary biological representation.

## Stage 6. Freeze and rerun gate

A replacement freeze is created only after:

- search truncation is resolved or explicitly bounded;
- taxa are validated;
- source-level eligibility and both spatial axes are adjudicated;
- mixed cases are retained;
- reviewer agreement is documented;
- climatic metrics are recomputed for the adjudicated set;
- the complete model, robustness, figure, manuscript, SI, and DOCX pipeline is rerun.

The old 34-species checksum and outputs are never silently overwritten.

## Current audit findings that motivate this protocol

From the archived 2026-08-02 systematic search and 2026-08-03 exploratory ecological artifact:

- 79,242 deduplicated bibliographic records were recovered;
- the search manifest reports 19 truncated query/database shards;
- the exploratory P1 candidate builder produced 588 raw two-word candidate names with directional evidence signals;
- 80 of those raw candidates had both within and geographic signals and were dropped as `mixed` before ecological modelling;
- some of those 80 are not valid taxa, confirming that taxon validation must precede state assignment;
- the exploratory post-GBIF binary dataset retained 120 unambiguous species (30 within, 90 among), of which 107 met the 20-climate-cell model threshold;
- those exploratory classifications remain explicitly `unreviewed` and must not be treated as a replacement for the human-adjudicated dataset.

These counts are diagnostics only. They are not used to adjudicate any species.
