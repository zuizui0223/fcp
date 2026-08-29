# JBI Chapter 1 spatial execution status

## Current state

Branch: `analysis/jbi-global-colour-boundaries`  
PR: #20, stacked on `audit/upstream-spatial-reaudit` / PR #17.

Primary inference remains:

**global flower-colour map → species-conditioned random-labelling test → species-level transition boundaries → shared-boundary strength → post-discovery geographic correspondence**.

Species may disappear from the display, but never from the null model: locations stay fixed and flower-colour observations are permuted strictly within species.

## Frozen empirical sample

- GitHub Actions reacquisition: **1,200 photographs = 6 species × 200**;
- calibration split: **480 = 80/species**;
- evaluation split: **720 = 120/species**;
- split is outcome-blind and hash-frozen;
- **720 evaluation photographs remain unopened for rule tuning**.

## Uniform Florence calibration features — complete

Run `33237848644` completed 24/24 deterministic calibration shards plus aggregation.

- 480/480 calibration feature records;
- 80/80 per species;
- evaluation rows opened: 0;
- final biological labels emitted: 0;
- feature JSONL SHA256: `968648b4d9a4516daa6fb938eaaa0000d665ed18b52801a2f4c6e639b68e7bed`.

`feature_ok` means Florence returned a box; it does **not** establish target-flower ROI validity.

## Reviewer-1 blind ROI/condition audit — complete, not final

Reviewer-1 rules were frozen before decisions. All 480 calibration crops were reviewed without geography, observer, date or colour candidate scores.

Overall ROI validity:

- usable: 370;
- rescue-segmentation candidate: 40;
- invalid: 39;
- ambiguous: 31.

Overall condition:

- fresh: 353;
- senescent: 27;
- damaged: 8;
- mixed/ambiguous: 22;
- not evaluable: 70.

The provisional **usable + fresh** set is **326/480**:

- `Antirrhinum majus`: 60;
- `Dactylorhiza sambucina`: 67;
- `Gentiana lutea`: 23;
- `Ipomoea purpurea`: 62;
- `Lysimachia arvensis`: 66;
- `Raphanus sativus`: 48.

Reviewer-1 is explicitly not an independent final adjudication.

## SAM2 rescue pilot — complete; scale-up rejected

A frozen paired pilot tested one reviewer-1 `rescue_segment` record plus one usable-fresh control for each rescue-bearing species using `facebook/sam2.1-hiera-tiny` with the existing Florence box as the only prompt.

Successful run: `33280003197`.

The before/after mask sheet was judged **before** species/arm mapping was opened. The frozen species-level scale-up gate required both rescue success and control preservation.

Result:

- controls preserved for all five tested species;
- rescue failed for `Antirrhinum`, `Dactylorhiza`, `Gentiana`, and `Raphanus`;
- `Ipomoea` rescue passed technically, but that record was reviewer-1 `mixed_or_ambiguous` condition and there were no remaining Ipomoea rescue records;
- remaining rescue records eligible for SAM2 scale-up: **0**.

Therefore the remaining 35 rescue candidates are **not** batch-rescued with SAM2 and are not silently promoted into the direct colour-calibration set.

## Feature geometry after reviewer-1 filtering — complete, provisional

Run `33280280728` repeated the 200-bootstrap GMM/BIC diagnostic on only the reviewer-1 usable-fresh 326 records. These components remain diagnostic feature-space structure, **not biological morph labels**.

Pre-filter → reviewer-1-filtered support:

- `Antirrhinum majus`: 3 components, bootstrap 0.93 → **0.805**;
- `Dactylorhiza sambucina`: 2 components, 1.00 → **1.00**;
- `Gentiana lutea`: 2 components, 1.00 → **0.995**;
- `Ipomoea purpurea`: 3 components, 0.89 → **0.955**;
- `Lysimachia arvensis`: 2 components, 1.00 → **1.00**;
- `Raphanus sativus`: observed BIC changes 2 → **3**, while bootstrap support for 4 components rises 0.785 → **0.815**.

Thus Raphanus complexity is **not** explained away by obvious reviewer-1 ROI/condition failures. A simple universal discrete-morph representation should not be forced at this stage.

## Reviewer-2 independent reblind package — complete and frozen for review

A second package was constructed from **all 480 calibration rows**, globally reshuffled using a new deterministic hash order and new `r2_id` values.

Generation run `33280339579` completed:

- 480/480 crops;
- 24 contact sheets;
- failures: 0;
- evaluation rows: 0;
- reviewer-2 labels created by script: 0;
- reviewer-facing artifact ID: `9722823556`;
- artifact ZIP SHA256: `1f045cbd174a88935c0ca23c0195ef00e5fb8f77606e4b0b98c9cda29321281e`.

The original run lost only its final git push race. Text records were regenerated without re-downloading images and matched the original run exactly:

- blank reviewer queue SHA256: `bb5df202ced26fe8fbb67f751067da49443adc8a7524226171d2ba9afe65140c`;
- hidden mapping SHA256: `e8108abafe1f46fb414e5648a979ece6b9f0684a438b08d153799542a610c809`;
- recovery run `33280626760`: success;
- recovery commit: `07b3d239f70d5715554be3fb0b466d8f98a78e3c`.

Reviewer-facing sheets expose only review order, `r2_id`, and crop. They hide species, original blind ID, reviewer-1 decisions, colour scores, geography, observer and date. The hidden mapping has **not been opened for reviewer-2 scoring**.

## Reconciliation rule — frozen before reviewer-2 decisions

Direct fresh colour calibration requires record-level agreement:

**reviewer-1 usable + fresh AND reviewer-2 usable + fresh**.

Any disagreement in ROI usability or fresh condition is withheld for third adjudication. A global post-hoc agreement threshold or majority rule will not be introduced after seeing reviewer-2 results.

Reconciliation implementation and tests are prepared but cannot be run until reviewer-2 decisions exist.

## Active gate

1. obtain a genuinely independent reviewer-2 scoring of artifact `9722823556` without exposing the hidden mapping or reviewer-1 decisions;
2. freeze the completed reviewer-2 queue before opening the mapping;
3. reconcile reviewer-1/reviewer-2 record by record and third-adjudicate disagreements;
4. rerun feature geometry on the final consensus fresh/evaluable calibration set;
5. freeze the colour representation and boundary rule — allowing a continuous within-species colour-vector representation where discrete morphs are not defensible;
6. only then open the 720 evaluation photographs;
7. only after held-out measurements exist run the species-conditioned spatial random-labelling and shared-boundary analyses.

## Not yet evaluated

- final flower-colour labels/states: 0/1,200;
- held-out evaluation measurements: 0/720;
- random vs non-random spatial placement: `not_evaluated`;
- species-level colour boundaries: `not_evaluated`;
- shared-boundary surface: `not_evaluated`;
- post-discovery geographic correspondence: `not_evaluated`.
