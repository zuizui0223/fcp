# JBI Chapter 1 spatial execution status

## Current state

Branch: `analysis/jbi-global-colour-boundaries`  
PR: #20, stacked on `audit/upstream-spatial-reaudit` / PR #17.

Primary inference remains:

**global flower-colour map → species-conditioned random-labelling test → species-level transition boundaries → shared-boundary strength → post-discovery geographic correspondence**.

Species may disappear from the display, but never from the null model: locations stay fixed and colour labels are permuted strictly within species.

## Frozen empirical sample

- iNaturalist development sample: **1,200 photographs = 6 species × 200**;
- calibration split: **480 = 80/species**;
- evaluation split: **720 = 120/species**;
- split is outcome-blind and hash-frozen;
- **720 evaluation photographs remain unopened for rule tuning**.

## 480-photo Florence calibration features — complete

GitHub Actions run `33237848644` completed 24/24 deterministic 20-photo shards and the aggregate job successfully.

- 480/480 calibration records extracted;
- 80/80 per species;
- 0 records from the evaluation split;
- 0 final biological labels emitted;
- Florence box returned for 480/480 records;
- feature JSONL: `data/calibration/jbi_ch1_florence_calibration_features_v1.jsonl`;
- feature JSONL SHA256: `968648b4d9a4516daa6fb938eaaa0000d665ed18b52801a2f4c6e639b68e7bed`;
- aggregate commit: `f31e53b08f7a0bdd006e43b24e2bdc01e250223c`.

Direct palette argmax values are diagnostic only. In particular, `Raphanus sativus` is represented by continuous anthocyanin-like and carotenoid-like visual axes rather than treating a single palette argmax as a validated morph classifier.

Important correction: **`feature_ok` means that a Florence box was obtained; it does not prove that the crop is a biologically valid target-flower ROI.** Blinded contact-sheet inspection has already exposed some overly broad/non-target crops, so formal ROI target-validity review is required before any state rule is frozen.

## Pre-condition feature geometry — complete, not morph labels

Run `33238769109` completed GMM/BIC diagnostics with 200 bootstraps while enforcing:

- calibration only;
- evaluation rows opened = false;
- final label = false;
- no geography/date/observer/environment input;
- candidate argmax not used as a training label;
- GMM components are **not** biological morph labels.

Pre-condition feature-space support:

- `Antirrhinum majus`: 3 components selected; bootstrap frequency 0.93;
- `Dactylorhiza sambucina`: 2 components; 1.00;
- `Gentiana lutea`: 2 components; 1.00;
- `Ipomoea purpurea`: 3 components; 0.89;
- `Lysimachia arvensis`: 2 components; 1.00;
- `Raphanus sativus`: observed BIC selects 2, but bootstrap selects 4 in 0.785, 3 in 0.18, and 2 in 0.035 — therefore the pre-condition structure is not stable enough to collapse into a simple morph rule.

These diagnostics must be repeated after independent ROI/flower-condition validation.

## Flower condition — automatic gate rejected; blind package ready

The quota-independent SigLIP pilot reproduced only **3/6** predeclared condition checks. It is therefore rejected as the automatic `fresh/senescent/damaged` gate; prompts and thresholds were not tuned post hoc.

A separate blinded review package is complete:

- workflow run `33238779685`: success;
- **480 blind rows, 24 contact sheets, 4 sheets/species**;
- artifact ID `9710738745`;
- failures: 0;
- geography, observer, date and colour candidate scores hidden from reviewer;
- labels created automatically: 0;
- evaluation rows opened: 0.

Allowed review outcomes are `fresh`, `senescent`, `damaged`, `mixed_or_ambiguous`, and `not_evaluable`.

## Active gate

1. review all 480 blinded crops for **target-ROI validity first**;
2. assign flower condition only when the target flower is evaluable;
3. exclude/mark non-target, overly broad and ambiguous ROIs without using colour state to make that decision;
4. repeat within-species feature-geometry diagnosis on the independently validated fresh/evaluable calibration subset;
5. freeze species-specific colour measurement/state/unresolved rules and hashes;
6. only then open the 720 evaluation photographs;
7. only after held-out colour states exist run the species-conditioned spatial random-labelling and shared-boundary analyses.

## Not yet evaluated

- final flower-colour states: 0/1,200;
- held-out evaluation measurements: 0/720;
- random vs non-random spatial placement: `not_evaluated`;
- species-level colour boundaries: `not_evaluated`;
- shared-boundary surface: `not_evaluated`;
- post-discovery geographic correspondence: `not_evaluated`.
