# H2 third-cohort prospective selection protocol

Date frozen: 2026-09-16 JST
Status: PREOUTCOME — no biological colour outcome may be opened

## Parent candidate universe

Use only the outcome-blind candidate universe frozen by `results/polymorphism_h2_third_cohort_candidate_frame_20260916/result.json`.

Lineage is fixed as:

- U0: 42,111 species;
- U100: 4,730 species with `after_observer_cap >= 100`;
- legacy high-depth species excluded in the parent P100 selection: 1,000 (500 discovery + 500 reserve);
- P100: 3,730 species;
- all prospectively selected P500 species excluded here: 500, including the one species that later reached only 98 fresh metadata rows;
- third-cohort candidate universe: 3,230 species;
- candidate-universe canonical CSV SHA256: `7fc0337074b55a91c0b50773a8d5c3ef82e877074c8b3cacc21f9af58e0ba77e`.

The parent P100/P500 selection artifact is workflow `34708044962`, artifact `10302477571`, digest `sha256:7e43763e3f1d9fa3ff108154dc6ef9443430b70fdb62532e2f6845582d607abb`.

## Frozen third-cohort sample size

Select exactly **500 species**.

This repeats the existing bounded high-depth resource design and is not tuned to a biological H2 outcome. If fewer than 500 eligible species exist, stop; do not lower the U100 gate or reuse excluded species.

## Deterministic outcome-blind ordering

For every species in the 3,230-species candidate universe compute:

`selection_hash = SHA256("FCP_H2_THIRD_COHORT_20260916_V1|" + inat_taxon_id + "|" + species)`

Sort ascending by `selection_hash`, then by integer `inat_taxon_id` as deterministic tie-break. The first 500 species are the frozen third cohort.

Frozen salt:

`FCP_H2_THIRD_COHORT_20260916_V1`

No family, genus, geography, sampled span, climate, morph, palette, D, white-axis score, H2 statistic, phylogeny, or any recovered P500 H2 information may enter selection.

## Independence rules

The selected set must have zero species/taxon overlap with:

1. all 1,000 legacy high-depth discovery/reserve species already removed upstream; and
2. all 500 P500 prospectively selected species, not merely the 499 species whose pixels were eventually opened.

The P500 sealed measurement table is not an input.

## Required selection receipt

Before any third-cohort image pixel is opened, persist:

- exact candidate count (3,230);
- exact selected count (500);
- ordered selected taxon IDs and species identities;
- selection hashes and ranks or a canonical selected CSV hash sufficient to reconstruct them;
- overlap with all excluded cohorts = 0;
- confirmation that only taxon identity and frozen opportunity capacity were available;
- confirmation that biological colour outcomes remain unopened.

## Next gate after selection

Selection alone does not authorize biological opening. The synthetic end-to-end technical qualification required by `POLYMORPHISM_H2_THIRD_COHORT_PREOPENING_FREEZE_20260916.md` must pass and be durably recorded before a one-run biological authorization can be issued.
