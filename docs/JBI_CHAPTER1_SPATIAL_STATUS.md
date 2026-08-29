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
- six-image calibration-only Copilot vision pilot: completed;
- six images x three independent Copilot passes = 18 valid repeatability responses;
- repeatability: flower visibility, flower condition, flower region, and segmentation feasibility were unanimous in 6/6 images;
- colour pattern and within-photo consistency were unanimous in 5/6 images;
- diagnostic colour scope was unanimous in 4/6 images;
- operational decision: only the 6/6-repeatable fields were accepted for automatic calibration screening; less-repeatable fields cannot establish final colour states;
- literature-constrained candidate colour codebook: frozen before 480-image screening;
- independent two-pass calibration-consensus workflow: implemented but not allowed to run until a complete first-pass calibration screen exists.

## Current empirical state

- development photographs: 1,200 (6 species x 200);
- calibration IDs: 480 frozen;
- evaluation IDs: 720 frozen and unopened for rule tuning;
- calibration technical materialization: 480/480 passed;
- six-image semantic repeatability diagnostic: 18/18 valid responses;
- complete 480-image semantic screen: **not completed**;
- final flower-colour classifications: 0/1,200;
- final species codebooks: not yet frozen;
- random versus non-random spatial placement: `not_evaluated`;
- species-level colour boundaries: `not_evaluated`;
- shared-boundary surface: `not_evaluated`;
- geographic correspondence: `not_evaluated`.

## Copilot quota stop

The first 480-image Copilot screening attempt failed before image processing because of a malformed candidate-codebook JSON; the JSON was corrected and a parse test was added.

The corrected screening then entered real image processing and produced valid records for the first several images in each shard. It subsequently stopped because the GitHub Actions identity reached its **monthly Copilot request quota**. A representative Raphanus shard processed eight calibration photographs successfully and then received `You have exceeded your monthly quota` on three consecutive retry attempts.

Consequences:

- this is a service/quota stop, not a biological or schema failure;
- partial records from aborted shards are **diagnostic only and are not promoted into the calibration dataset**;
- the 480-image dataset will not be assembled by mixing pre-quota Copilot records with a different later method;
- repeatedly rerunning the same Copilot workflow is not a valid recovery path while the monthly quota remains exhausted;
- the 720-image evaluation set remains unopened.

GitHub Models cannot provide a separate fallback because GitHub retired that service on 2026-07-30. The active fallback therefore removes paid/request-quota model inference from the main pipeline.

## Quota-independent image route — active

A six-image calibration-only pilot is running with **Florence-2-base-ft** plus deterministic pixel-colour extraction:

1. Florence-2 open-vocabulary detection localizes a flower ROI;
2. a fixed sRGB reference palette, declared before the pilot result, quantifies colours inside that ROI;
3. species-specific candidate-state scores are computed deterministically from the frozen literature codebook;
4. the result is compared with five predeclared diagnostic expectations from the already-completed six-image pilot;
5. the known senescent `Gentiana lutea` image is retained as a negative control: colour extraction may succeed but cannot validate a fresh-flower state.

Files:

- `scripts/data/run_jbi_ch1_florence_colour_pilot.py`;
- `docs/supporting/jbi_ch1_florence_pilot_expected_v1.json`;
- `.github/workflows/jbi-ch1-florence-colour-pilot.yml`.

The Florence model is used to localize the flower; the named colour suggestion is produced from explicit numeric pixel features rather than unconstrained language generation. This pilot remains `final_label=false`.

## Repeatability decision

No numerical acceptance threshold was estimated after seeing the six-image diagnostic. The prior Copilot diagnostic established only which semantic fields were repeatable enough to use as **screening descriptors**, not final labels.

Previously stable semantic descriptors:

- flower visibility;
- flower condition (`fresh`, `senescent`, `damaged`, `mixed_or_ambiguous`);
- flower region;
- segmentation feasibility.

Not accepted as automatic final decisions:

- free-text colour terms;
- colour pattern;
- diagnostic colour scope;
- multiple-flower colour consistency;
- candidate biological colour state.

## Literature-constrained candidate states

Candidate states were declared before the 480-image screening in:

`docs/supporting/jbi_ch1_species_colour_candidate_codebook_v1.json`

Current candidate contrasts:

- `Ipomoea purpurea`: white / pink / blue-purple, with explicit acknowledgement of its multi-locus flower-colour genetics;
- `Raphanus sativus`: white / yellow / pink / bronze;
- `Gentiana lutea`: yellow / orange;
- `Dactylorhiza sambucina`: yellow / purple;
- `Antirrhinum majus`: magenta-pseudomajus-like / yellow-striatum-like / intermediate-or-other, based on whole-corolla pigment distribution rather than hue alone;
- `Lysimachia arvensis`: blue / red.

Every species also has `unresolved`. New states cannot be invented by the image model.

## Next valid gate

1. finish the six-image Florence localization/colour pilot;
2. inspect localization success and the five predeclared diagnostic state checks;
3. if the open-model route is supported, rerun all 480 calibration images under one uniform quota-independent implementation;
4. quantify state separability and unresolved/localization-failure rates by species;
5. independently verify candidate state-bearing calibration records and freeze final species-specific measurement rules;
6. only after that freeze open the 720-image evaluation set;
7. only after held-out colour states exist run the species-conditioned spatial random-labelling analysis.
