# Polymorphism species attributes — Step 4 life-form source protocol

Date frozen: 2026-09-10 JST
Branch: `analysis/polymorphism-species-attributes-step4b`

## Why this block exists

The Step 4 preflight found no tracked FCP life-form source, but its outcome-blind GIFT v3.2 metadata scan identified an exact categorical life-form trait:

- GIFT level: `2.3.1`
- GIFT name: `Life_form_1`
- categories: `phanerophyte, chamaephyte, hemicryptophyte, cryptophyte, therophyte`
- metadata type: categorical
- metadata record count: 100390

This source block is frozen **before any new D/covariate association is computed**.

## Frozen GIFT source

Use GIFT stable version `3.2`, direct raw records only.

Reuse immutable source metadata from the already-completed `zuizui0223/island` GIFT campaign:

- workflow run: `29338693288`
- artifact name: `gift-direct-traits-29338693288`
- artifact ID: `8313108327`
- artifact digest: `sha256:ae660ca791bab115f21d8b2ce7c9fdf65446c66365b274b6e33a7af1ada9ccca`

The artifact supplies the pinned GIFT v3.2 versions, trait metadata, reference-trait map, reference metadata, and GIFT `work_ID -> work_species` taxonomy table.

## Direct-record acquisition rule

Query only `trait_ID = 2.3.1` from public, unbiased GIFT reference/trait pairs discovered in the pinned `reference_traits` metadata.

For each pair query:

`traits_raw?traitid=2.3.1&deriv=0&biasderiv=0&refid=<REF_ID>`

Admit a raw record only when all of the following hold:

1. returned `ref_ID` and `trait_ID` equal the requested pair;
2. `derived=0` and `bias_deriv=0` (or equivalent false flags);
3. the reference is public (`restricted=0`) and has a citation;
4. `matched=1` and `resolved=1`;
5. no `cf_genus`, `cf_species`, or `aff_species` uncertainty flag;
6. no infraspecific `subtaxon` scope;
7. the resolved GIFT `work_species` is a binomial;
8. every token in `trait_value` belongs to the frozen five-category ontology.

No logical/taxonomic derivation and no LLM inference is allowed.

## Species-level standardization

The FCP denominator is the same frozen 369-species discovery frame.

For each FCP species:

- collect all admitted direct `Life_form_1` tokens across public references;
- if the union contains exactly one category, assign that category;
- if the union contains more than one category (including slash-combination records), assign `mixed`;
- if there is no admitted direct record, leave life form missing.

This rule is fixed before D is opened.

## Prospective analysis gate

Life form enters the later D association analysis only if:

1. at least 50 of the 369 FCP species have an admitted direct life-form assignment, and
2. at least two standardized categories each have `n >= 10`.

The later omnibus analysis includes only categories with `n >= 10`; the threshold rule is deterministic and outcome-blind.

If the gate fails, the life-form family closes with `p=1` in the predeclared six-slot Holm correction. It may not be rescued by switching to `Life_form_2`, growth form, web search, or LLM trait inference after D outcomes are opened.

## Outcome firewall

This source-acquisition runner may read only the frozen species names from the D table. It must not read the `D` column, compute any D association, or inspect D-stratified coverage.
