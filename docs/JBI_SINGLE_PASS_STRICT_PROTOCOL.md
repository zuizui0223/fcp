# JBI v2.2 strict coexistence–segregation protocol

## Status

This protocol supersedes the forced `within_population` / `among_population` classification for the replacement analysis. Independent duplicate review was intentionally waived on 2026-08-26; the old reviewer infrastructure remains archived as optional sensitivity/provenance material.

The replacement analysis is a **single-pass conservative documented-evidence audit**. The historical 34-species binary freeze is never overwritten and is retained as a historical sensitivity dataset.

## Central biological question

The object of interest is not whether intraspecific floral-colour variation is generically “within” or “among” populations. It is whether discrete colour polymorphism is documented to:

- **C — coexist locally**: two or more natural floral-display colour variants occur in the same population/site;
- **S — segregate spatially**: morph occurrence or morph frequency is spatially structured among populations/sites/regions/islands or along an explicit geographic gradient.

C and S are **independent positive-evidence axes**.

Therefore:

| C | S | documented organization state |
|---:|---:|---|
| 1 | 0 | `local_coexistence_only` |
| 0 | 1 | `spatial_segregation_only` |
| 1 | 1 | `coexistence_and_segregation` |
| 0 | 0 | `unresolved` |

`C=0` does **not** mean that morphs biologically fail to coexist. `S=0` does **not** mean spatial homogeneity. A zero means only **not documented by the strict pass**.

The C=1,S=1 state is biologically important: morphs can coexist in local populations while their frequencies or occurrence remain structured across the wider range.

## Retrieval boundary

The fixed search surface is OpenAlex title/abstract OQL v2.2:

- 15 query blocks;
- 13,911 query memberships;
- 12,064 deduplicated works;
- zero truncated v2.2 blocks;
- 34/34 historical classification source records recovered.

All 12,064 records remain in the audit universe. Automatic taxon extraction cannot exclude a record from that universe.

## Positive-evidence philosophy

A C or S axis may be set to 1 only when explicit source wording satisfies its positive criterion. Missing wording, generic polymorphism terminology, or an abstract that does not resolve spatial organization is never converted into the opposite biological state.

Ambiguous eligibility, conflicting source context, continuous-only colour variation, non-primary evidence, community-level comparisons without an intraspecific signal, or ambiguous taxon attribution remain unresolved.

## Source eligibility before C/S coding

A source can contribute a C/S positive only when all of the following hold:

1. primary empirical source type (`article`, `dissertation`, `preprint`, or `report`);
2. explicit floral-display colour context;
3. explicit discrete colour-polymorphism/morph context;
4. no hard artificial/cultivated/induced/ontogenetic conflict in the same record;
5. not a community-level/inter-specific colour summary unless the record also contains an explicit intraspecific polymorphism signal.

Hard conflicts include cultivar/horticultural-line studies, breeding/mapping populations, transgenic or induced variation, tissue culture, artificial-flower arrays, experimental populations, and ontogenetic colour change.

Reviews, meta-analyses, summaries, definitions of FCP, and background statements do not establish C or S for a focal species.

## C — local coexistence positive criterion

`C_local_coexistence_documented = 1` only when the source explicitly establishes discrete floral-colour variants at the same local unit.

Sufficient examples include:

- colour morphs `co-occur` or `coexist` in a population/site;
- named colour forms are found in the `same population`;
- a study population explicitly contains two or more colour phenotypes;
- multiple morphs are reported in each of several natural populations;
- morphs are explicitly stated to vary `within` populations;
- populations are explicitly described as colour-mixed.

Insufficient by itself:

- `colour morph(s)`;
- `polymorphic` / `polymorphism`;
- a definition saying that FCP means morphs occur within populations;
- common-garden or artificial-array coexistence;
- multiple colours without a population/site unit;
- a general statement that selection can maintain polymorphism within populations.

Negated coexistence statements are not positive C evidence.

## S — spatial segregation positive criterion

`S_spatial_segregation_documented = 1` only when the source explicitly establishes spatial structuring of discrete colour morphs or their frequencies.

Sufficient examples include:

- flower-colour or morph-frequency variation explicitly `among`, `between`, or `across` populations/sites/regions/islands;
- explicit geographic/spatial/altitudinal/latitudinal variation or distribution **of flower colour or colour polymorphism**;
- morph frequencies that vary among geographic units;
- colour morphs that are regionally restricted or replaced;
- polytypism among populations;
- a documented mixture of monomorphic and polymorphic populations across the range;
- colour polymorphism explicitly documented both within and among populations.

Insufficient by itself:

- broad range size;
- many sampling sites;
- spatial variation in a different response such as reproductive success, climate, or pollinator abundance;
- `consistent` morph ratios among populations;
- comparisons to fixed differences in other species;
- generic geographic terminology with no spatially resolved colour contrast.

Negated spatial-difference statements are not positive S evidence.

## Taxon resolution

Taxon resolution is attempted only for C/S-positive sources and for the 34 exact historical source records retained for provenance.

Candidate names are ranked from title occurrences first, then the first 2,500 abstract characters. Names are checked against GBIF and accepted-name resolution. A source contributes to species-level aggregation only when it resolves unambiguously to one accepted plant species.

The historical manifest may be used only as an exact **source → canonical taxon rescue map**. Historical spatial labels, climate-cell counts, and model results are not supplied to the C/S evidence builder.

## Species aggregation

Across taxonomically resolved positive sources:

- any valid C-positive source → `C=1`;
- any valid S-positive source → `S=1`.

Evidence for C and S may come from different papers and different parts of the species range.

The replacement strict freeze contains only species with at least one positive axis. Unresolved records and unresolved taxa are retained separately rather than force-classified.

## Ecological inference

The primary analysis fits the two documented-evidence outcomes separately:

- probability that local coexistence is documented, `P(C=1)`;
- probability that spatial segregation is documented, `P(S=1)`.

A `coexistence_and_segregation` species is positive in both models.

If cell sizes are adequate, a secondary non-ordinal model compares:

`local_coexistence_only` / `spatial_segregation_only` / `coexistence_and_segregation`.

The historical within-vs-among binary analysis is retained only as a sensitivity analysis.

## Reporting limitation

The replacement evidence audit is single-pass and rule-conservative rather than independently duplicate-reviewed. It prioritizes explicit reproducible positive evidence over recall. Consequently, unresolved and zero-coded axes must never be described as demonstrated biological absence.
