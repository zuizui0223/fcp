# JBI v2.2 strict single-pass documented-evidence protocol

## Status

This protocol supersedes the duplicate-independent-review path for the replacement analysis. Independent Reviewer 1/Reviewer 2 calibration was intentionally waived on 2026-08-26. The old review infrastructure remains archived as optional sensitivity/provenance material.

The replacement analysis must therefore be described as a **single-pass conservative documented-evidence audit**, not as a double-reviewed systematic review.

The historical 34-species freeze is never silently overwritten.

## Retrieval boundary

The fixed retrieval surface is OpenAlex title/abstract OQL v2.2:

- 15 query blocks;
- 13,911 query memberships;
- 12,064 deduplicated works;
- zero truncated v2.2 blocks;
- 34/34 historical classification source records recovered.

All 12,064 records remain in the audit universe. Automatic taxon extraction cannot exclude a record.

## Evidence philosophy

The replacement path uses **positive evidence only**.

A spatial axis may be set to 1 only when explicit source wording satisfies a strict positive criterion. Failure to match a criterion is interpreted as **not documented by this pass**, never as biological absence.

Ambiguous eligibility, conflicting source context, or ambiguous taxon attribution remains `unresolved`.

## Natural/display eligibility

A source is `eligible_high_confidence` only when title/abstract material contains:

1. floral-display colour/reflectance/pigmentation context;
2. natural/field/population/geographic context;
3. a colour-variation signal relevant to the spatial audit;
4. no explicit hard conflict such as cultivar-only, horticultural-line-only, induced mutation/transgenic/editing/tissue-culture, or ontogenetic colour-change framing.

A hard conflict does not automatically become an exclusion when other evidence may exist; it is stored as `conflict_unresolved`.

Sources with display relevance but insufficient natural-variation evidence remain `display_relevant_unresolved`.

## Variation form

The pass records `discrete`, `continuous`, `both`, or `unclear` signals.

Local coexistence requires discrete variation. Continuous variation is retained in the audit table but is not silently converted to a discrete colour-morph state.

## Local coexistence positive criterion

`local_coexistence_documented = 1` only when floral-display/discrete context is accompanied by explicit wording equivalent to one of the following:

- colour variants/morphs co-occur or coexist within a population/site;
- the same population/site explicitly contains multiple colour morphs/forms;
- an equivalent same-population statement with direct morph evidence.

Insufficient by itself:

- `color/colour morph(s)`;
- `polymorphic` or `polymorphism`;
- multiple colour forms with no population resolution;
- sampling many sites;
- common-garden coexistence without natural same-population evidence.

## Geographic structure positive criterion

`geographic_structure_documented = 1` only when floral-display colour context is accompanied by explicit spatial differentiation such as:

- differences among/between populations, sites, localities, regions, or islands;
- geographic/geographical/spatial/latitudinal/altitudinal/elevational variation, differentiation, structure, cline, pattern, or gradient;
- colour/morph-frequency variation explicitly among geographic units;
- explicit regional/geographic restriction, replacement, or segregation of colour forms.

Insufficient by itself:

- broad species range;
- many sampling localities with no reported colour contrast;
- environmental association without spatially resolved colour differentiation;
- generic colour polymorphism.

## Taxon resolution

Automatic title/abstract binomial hints are used only for potentially informative records.

Names are checked against GBIF and accepted-name resolution. A source is assigned to a species only when the source resolves unambiguously to one accepted plant species.

If multiple accepted taxa are plausible, the source remains taxonomically unresolved.

The historical manifest may be used only as an exact **source -> canonical taxon rescue map** for its 34 known source records. Historical spatial labels, climate-cell counts, model results, and effect directions are not supplied to the strict evidence builder.

## Species aggregation

Across high-confidence eligible, taxonomically resolved sources:

- any local positive source -> `local_coexistence_documented = 1`;
- any geographic positive source -> `geographic_structure_documented = 1`.

Then:

- local=1, geographic=0 -> `within_evidence_only`;
- local=0, geographic=1 -> `among_evidence_only`;
- local=1, geographic=1 -> `mixed_evidence`;
- neither axis documented -> `unresolved`.

`mixed_evidence` is non-ordinal. Its two axes may come from different papers and different parts of the species range.

## Freeze rule

The first replacement freeze is a **strict high-confidence documented-evidence freeze**.

Primary state analyses use only species with at least one positive documented spatial axis and no unresolved taxonomic assignment for the contributing positive source(s).

Species/sources with conflicts or insufficient evidence remain in explicit unresolved tables and are used for sensitivity/bounds analyses rather than force-classified.

The historical 34-species binary freeze remains a separate historical sensitivity dataset.

## Inference rule

Primary ecological models fit the two positive documented-evidence axes separately. A mixed species is positive on both axes.

A secondary non-ordinal multinomial model may compare within-only / among-only / mixed when all states have adequate sample size.

Observation-process variables such as number of sources per species and source availability should be retained because these are documented-evidence outcomes rather than perfect biological detection.

## Reporting limitation

Methods and Discussion must explicitly state that the replacement evidence audit is single-pass and rule-conservative rather than independently duplicate-reviewed. The tradeoff is deliberate: higher reproducibility and explicit positive-evidence thresholds, but potentially lower recall and no inter-reviewer agreement estimate.
