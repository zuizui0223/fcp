# RGFCA Step 8L — literature-derived flower-colour polymorphism registry

Date: 2026-09-10 JST
Status: prospectively frozen before using literature labels to evaluate the expanded photo-first atlas.

## Purpose

Build an independent registry of plant species for which primary literature supports naturally occurring intraspecific floral-display colour polymorphism. This registry is not a species-selection frame for RGFCA. It is an external evidence layer for validating photo-derived polymorphism, auditing literature coverage, and identifying mechanism-rich species.

## Separation from the atlas

The Step-8 photo-first atlas uses only the outcome-blind 4,730-species capacity frame. Literature-positive status, C/S labels, article counts, mechanism terms, and historical 34-species membership must never alter atlas species inclusion, photo selection, colour thresholds, or C*/S* state construction.

After both sides are frozen independently, the registry may be joined to the atlas for validation and interpretation.

## Literature surfaces

Use two retained surfaces with different roles.

1. Broad archived discovery surface: the 79,242-record systematic map and its retained legacy shards. This is a candidate-discovery surface only because some legacy shards were truncated and it is not the canonical completeness boundary.
2. Strict completeness surface: the OpenAlex v2.2 title/abstract universe of 12,064 deduplicated works from 15 conceptual query blocks, with zero truncated v2.2 blocks and recovery of all historical classification sources/review seeds. This is the primary evidence-screening surface.

Do not describe the 79,242-record map as a complete systematic-review denominator.

## Primary positive state P

`P_natural_flower_colour_polymorphism_documented = 1` only when retained source text supports at least two discrete, naturally occurring, intraspecific floral-display colour variants in wild or naturalised populations.

Sufficient examples include explicit flower-colour polymorphism, multiple naturally occurring flower-colour morphs/forms, population morph-frequency observations, or clearly discrete natural floral-display colour variants.

Insufficient alone:

- cultivar/horticultural/breeding lines only;
- induced mutants, transgenics, QTL/mapping populations without natural polymorphism evidence;
- ontogenetic/post-pollination colour change;
- interspecific colour differences;
- continuous colour variation without evidence for discrete morphs;
- anther/pollen/androecium colour without floral-display variation;
- experimental arrays where natural occurrence is not established.

Ambiguous records remain unresolved rather than negative.

## Secondary evidence axes

For each P-positive species retain, when documented:

- `C_local_coexistence_documented` from the existing strict v2.2 C definition;
- `S_spatial_segregation_documented` from the existing strict v2.2 S definition;
- `mechanism_pollinator`;
- `mechanism_abiotic_climate`;
- `mechanism_herbivory`;
- `mechanism_mating_system_selfing`;
- `mechanism_genetic_neutral_gene_flow`;
- `mechanism_pigment_molecular`;
- `mechanism_other_or_unclear`.

Mechanism fields are evidence tags, not causal conclusions.

## Evidence tiers

- Tier A: direct natural-population statement plus discrete floral colour morph evidence in title/abstract or retained primary-source evidence.
- Tier B: direct floral polymorphism evidence with natural occurrence supported elsewhere in the retained source aggregation.
- Tier C: plausible candidate requiring full-text/manual confirmation.
- Excluded: artificial/horticultural, ontogenetic-only, interspecific, continuous-only, or otherwise incompatible.

Only Tier A/B count as literature-positive in atlas validation. Tier C stays visible but outside positive denominators.

## Taxonomy

Resolve binomials to accepted plant species using the retained GBIF taxon-resolution path. Merge synonyms only after accepted-name resolution. Preserve source names and accepted names. Ambiguous taxon assignments remain unresolved.

## Outputs

Create a species registry with at least:

`accepted_name, source_name, gbif_taxon_key, family, P_documented, evidence_tier, C_documented, S_documented, n_supporting_works, supporting_work_ids, mechanism_tags, evidence_excerpt_or_summary, taxon_resolution_status`.

Also preserve a source-level evidence table so every positive species is auditable back to individual works.

## Atlas validation after independent freeze

Against the expanded photo-first atlas report separately:

1. recall of Tier-A/B literature-positive species by photo-derived polymorphism thresholds (including second morph >=10% and >=20% and continuous D);
2. distribution of D among literature-positive versus literature-unresolved atlas species without treating unresolved as biological negatives;
3. correspondence of literature C/S with photo-derived C*/S* where both are estimable;
4. number and identities of high-D photo-derived species with no retained literature-positive record — explicitly labelled candidate novel/under-documented polymorphism, not new biological discovery until checked;
5. fraction of literature-positive species absent from or ineligible for the 4,730-species photo-capacity frame, to quantify citizen-science coverage bias.

## Existing precursor

The retained review-queue builder already classified 664 candidate species and recorded 86 as `natural_polymorphism`, with 72 selected for prioritized review and four known-positive controls recovered 4/4. These are precursor screening outputs, not yet the final registry.

## Claim boundary

The registry estimates documented literature evidence, not the true prevalence of flower-colour polymorphism across angiosperms. Literature-unresolved does not mean monomorphic. The photo-first atlas and literature registry remain independently constructed until the validation join.