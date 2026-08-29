# JBI Chapter 1 spatial execution status

## Active branch and PR

- branch: `analysis/jbi-global-colour-boundaries`
- stacked base: `audit/upstream-spatial-reaudit` (PR #17)
- active PR: #20

## Primary question

The primary Chapter 1 inference is now:

**global flower-colour map → species-conditioned random-labelling test → species-level transition boundaries → shared-boundary strength → post-discovery geographic correspondence**.

The previous C/S-to-climate comparison is retained as secondary/provenance context and is not the primary spatial discovery analysis.

## Completed on this branch

- species-conditioned permutation core: implemented;
- no between-species colour-label exchange: guarded by tests;
- label-independent boundary detectability denominator `A(x)`: implemented;
- `A(x)=0` / insufficient support: treated as not evaluable rather than biological zero;
- detectability invariance under colour-label permutation: guarded;
- Chapter 1 protocol: frozen;
- deterministic 6 x 200 iNaturalist development-photo acquisition: completed through GitHub Actions;
- acquisition QC: all six species passed observer/spatial/season coverage gates;
- 480/720 calibration/evaluation split: hash-frozen before colour measurement;
- split assignment uses species + stable photo ID only;
- frozen acquisition re-run is locked unless explicit replacement is requested;
- 480 calibration images: 480/480 materialized and decoded successfully;
- conservative technical image flags: 0/480;
- evaluation images opened during technical audit: 0/720;
- six-image calibration-only vision pilot: completed;
- six images x three independent passes = 18 valid repeatability responses;
- repeatability: flower visibility, flower condition, flower region, and segmentation feasibility were unanimous in 6/6 images;
- colour pattern and within-photo consistency were unanimous in 5/6 images;
- diagnostic colour scope was unanimous in 4/6 images;
- operational decision: only the 6/6-repeatable fields are accepted for automatic calibration screening; less-repeatable fields remain advisory and cannot establish final colour states;
- literature-constrained candidate colour codebook: frozen before 480-image semantic screening;
- 480-image semantic screening workflow: implemented and running calibration-only; candidate states remain screening suggestions, not final labels.

## Current empirical state

- development photographs: 1,200 (6 species x 200);
- calibration IDs: 480 frozen;
- evaluation IDs: 720 frozen and unopened for rule tuning;
- calibration technical materialization: 480/480 passed;
- calibration semantic repeatability diagnostic: 18/18 valid responses on six pilot images;
- final flower-colour classifications: 0/1,200;
- final species codebooks: not yet frozen;
- random versus non-random spatial placement: `not_evaluated`;
- species-level colour boundaries: `not_evaluated`;
- shared-boundary surface: `not_evaluated`;
- geographic correspondence: `not_evaluated`.

## Repeatability decision

No numerical acceptance threshold was estimated after seeing the six-image diagnostic. Instead, only fields that were completely stable across all three independent passes for all six pilot images are carried forward for automatic **screening**, not final labelling.

Accepted for automatic screening:

- flower visibility;
- flower condition (`fresh`, `senescent`, `damaged`, `mixed_or_ambiguous`);
- flower region (`single_target_clear`, `multiple_flowers_clear`, etc.);
- segmentation feasibility.

Not accepted as automatic final decisions:

- free-text colour terms;
- colour pattern;
- diagnostic colour scope;
- multiple-flower colour consistency;
- candidate biological colour state.

The latter fields may be recorded conservatively, but discordance forces `unresolved` rather than a biological state.

## Literature-constrained candidate states

Candidate states are declared before the 480-image semantic screen in:

`docs/supporting/jbi_ch1_species_colour_candidate_codebook_v1.json`

Current candidate contrasts:

- `Ipomoea purpurea`: white / pink / blue-purple, with explicit acknowledgement of its multi-locus flower-colour genetics;
- `Raphanus sativus`: white / yellow / pink / bronze;
- `Gentiana lutea`: yellow / orange;
- `Dactylorhiza sambucina`: yellow / purple;
- `Antirrhinum majus`: magenta-pseudomajus-like / yellow-striatum-like / intermediate-or-other, based on whole-corolla pigment distribution rather than hue alone;
- `Lysimachia arvensis`: blue / red.

Every species also has `unresolved`. The screening model is not allowed to invent a new state.

## Current active gate

The active task is **480-photo semantic calibration screening**, not holdout evaluation and not spatial inference.

For each frozen calibration image the workflow records:

1. flower visibility;
2. flower condition;
3. flower-region/segmentation feasibility;
4. within-photo consistency as a conservative unresolved gate;
5. one literature-predeclared candidate state or `unresolved`.

All records are explicitly marked `screening_only=true` and `final_label=false`.

The first semantic-screen attempt failed before opening any image because the candidate-codebook JSON contained a syntax error. The codebook was corrected and a JSON-parse test was added before re-running the workflow. Therefore that failure did not generate or contaminate image labels.

## Next gate after screening

1. aggregate all 480 screening records and quantify usable vs unresolved yield by species;
2. inspect unresolved causes and candidate-state support;
3. independently re-score state-bearing calibration records before final codebook freeze;
4. freeze species-specific final measurement rules and hashes;
5. only after that freeze open the 720-image evaluation set;
6. only after held-out colour states exist run the species-conditioned spatial random-labelling analysis.
