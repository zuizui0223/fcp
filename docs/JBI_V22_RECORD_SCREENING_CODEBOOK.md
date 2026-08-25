# JBI v2.2 blind record-screening codebook

## Purpose

This stage screens the **12,064 deduplicated v2.2 bibliographic records before taxon validation**. It exists to prevent automatic taxon extraction, historical labels, or search-query membership from silently determining inclusion. Reviewers must not see the coordinator key, old within/among labels, climatic values, or manuscript model results while coding.

This stage does **not** code within / among / mixed. Spatial evidence is coded later at source level after a record is retained and the focal taxon is resolved.

## Review unit

One bibliographic record (title + abstract where available). If title/abstract are insufficient, mark `full_text_required = yes`; do not force inclusion or exclusion from missing information.

## Reviewer fields

### `record_relevance`

Allowed values: `include`, `exclude`, `uncertain`.

Use **include** when the available record supports, or clearly studies, naturally occurring intraspecific variation in floral-display colour or colour morphs and should proceed to taxon/full-source review.

Use **exclude** only when the available record clearly falls outside the target evidence domain.

Use **uncertain** when the title/abstract is incomplete, ambiguous, or a full source is needed to decide.

### `natural_intraspecific_variation`

Allowed values: `yes`, `no`, `uncertain`.

`yes` requires evidence that variation occurs among conspecific individuals or populations in nature/wild material. Wild-origin material may remain `uncertain` if the record does not establish that the colour variation itself occurs naturally.

`no` when the record is clearly limited to interspecific comparison, cultivars/horticultural lines without natural evidence, induced mutation/transgenics/gene editing/tissue culture, or ontogenetic colour change within a flower with no among-individual/population variation.

### `floral_display_colour`

Allowed values: `yes`, `no`, `uncertain`.

`yes` for colour/reflectance/pigmentation of floral display surfaces such as petals, corolla, tepals, labellum, or showy floral bracts.

`no` when the evidence is only vegetative colour, fruit colour, or a non-display floral structure with no display-colour question.

Use `uncertain` rather than inferring display function from organ name alone.

### `focal_taxon_text`

Copy the focal taxon name(s) as presented by the source. Do not silently expand abbreviated genera or replace names with accepted taxonomy at this stage. Multiple taxa may be separated by semicolons.

### `full_text_required`

Allowed values: `yes`, `no`.

Use `yes` whenever title/abstract do not support a reliable record-level decision, including missing abstracts where the title alone is insufficient.

### `exclusion_reason`

If `record_relevance = exclude`, use one primary reason where possible:

- `interspecific_only`
- `cultivar_horticultural_only`
- `induced_or_transgenic`
- `ontogenetic_change_only`
- `non_floral_colour`
- `non_display_floral_trait`
- `no_colour_variation_question`
- `not_primary_biological_evidence`
- `other_clear_exclusion`

Add details in `notes` when needed.

## Conservative rule

When evidence is incomplete, choose **uncertain/full-text-required**, not exclusion. A zero or absence of wording in an abstract is never treated as evidence that a biological state is absent.

## Duplicate independent review

Reviewer 1 and reviewer 2 code the same records independently. Neither sees the other's decisions before completion. Use separate copies of the blind sheet: Reviewer 1 fills only `reviewer_1_*` fields and Reviewer 2 fills only `reviewer_2_*` fields. Adjudication fields remain blank until both independent copies are locked.

The coordinator key must not be distributed to either reviewer before independent coding is complete.

Disagreements are adjudicated after raw agreement and Cohen's kappa are calculated for:

1. record relevance;
2. natural intraspecific variation;
3. floral display colour;
4. full-text requirement.

`focal_taxon_text` is additionally checked by normalized exact agreement, but taxonomy is not silently normalized during record screening.

The same person or AI agent cannot count as two independent reviewers.

## Wave 0 calibration gate

Before the full 12,064-record screen, the codebook is calibrated on a blinded **384-record Wave 0** sample containing:

- all 34 historical benchmark source records, with benchmark status hidden from reviewers;
- 100 records with no automatically detected binomial;
- 100 records with an automatically detected binomial;
- 100 non-English records;
- 50 records with missing abstracts.

These strata are mutually exclusive in the calibration sample. Stratum identity, query membership, historical status, and automatic taxon hints are coordinator-only.

The codebook is revised if disagreements reveal a systematic interpretation problem. A second calibration wave should be used after a material rule change rather than silently changing rules during the full screen.

## Transition to source-level spatial review

Only after record-level adjudication should retained/uncertain records proceed to:

1. focal taxon resolution against a taxonomic authority;
2. full-source eligibility review;
3. independent coding of `local_coexistence_documented` and `geographic_structure_documented`;
4. species-level aggregation into within-only, among-only, mixed, or unresolved evidence states.

No climatic values or historical model results are consulted before those evidence decisions are frozen.
