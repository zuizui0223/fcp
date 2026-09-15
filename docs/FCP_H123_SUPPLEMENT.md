# FCP H1–H3 companion supplement

Working evidence inventory — incomplete, not submission-ready.

## S1. H1 source hierarchy and chronology correction

The direct observer-disjoint protocol was committed as
743e13d3597c88ea7a725c596ea733b30d978024 at 2026-09-13 01:53:54 JST.
The repeated-partition protocol was committed as
92488ada5f4f535f6731bba6519101d8ed8a7e1a at 2026-09-13 02:03:49 JST.
Thus the older evidence-ledger description of the repeated-partition protocol
as the earliest protocol is not used in this manuscript. Commit chronology
alone is not authenticated evidence of every historical non-access event.

Direct result: `results/polymorphism_h1_observer_disjoint_d_20260913/result.json`,
copied without modification from c34c5e85504bee76aba7c5f780ef162863aef09e.
CI run 34706708406, artifact 10302770009. This is the direct current-D result,
not the historical approximate 0.971 summary.

Other H1 analyses remain retained:
`results/polymorphism_h1_observer_disjoint_reliability_20260913/result.json`
and `results/polymorphism_h1_observer_disjoint_D_reliability_20260913/result.json`.
Their different partition rules and the latter's failed threshold must remain
visible. No frozen result or threshold is changed by this correction.

## S2. H2 source and retrospective status

`results/polymorphism_white_axis_targeted_test_20260912/result.json` records the
targeted primary/strict results and its retrospective claim boundary.
CI run 34675697583, artifact 10292077657. The coarse-state-preserving null
must not be replaced with the easier isotropic comparison. Full mode/null
construction and residual diagnostics require dedicated reproducibility tables.

## S3. H3 explanatory boundaries

H3a source inventory:
`results/polymorphism_h3a_phylogenetic_signal_20260912/frozen_result_manifest.json`.
H3b result note: `docs/POLYMORPHISM_H3B_RESERVE_SPAN_RESULT_FREEZE_20260912.md`.
The manuscript figures must include reserve results and preserve all placement
scenarios. Missing trait coverage is not a biological negative.

The CI artifact summary tables were retrieved without rerunning either test:

- H3a: run 34677042793 / artifact 10292218669; retained as
  `docs/supporting/h123/h3a_signal_by_scenario.csv`. SHA256
  `8b37b0edcdf0823351aad82d443d2ff6686d3684294855f20900e1770614f986`
  matches the existing frozen manifest and retrieved artifact bytes.
- H3b: run 34677468362 / artifact 10292399238; retained as
  `docs/supporting/h123/h3b_span_summary.csv`. Retrieved SHA256
  `a8f58d2cf7282ec8aa55fd5741eb34e9440112fd319710654c4f03c12551460b`.

These copies use LF newlines; CI checks LF-canonicalized hashes to accommodate
Windows checkout conversion. This does not rewrite any original exact-artifact
contract. Manuscript H3 values are checked against these tables. This verifies
reporting, not independent reconstruction of tree models or permutation draws.

## S4. Measurement failures and unqualified gates

Retain `docs/RGFCA_MONARDA_REGION_AGREEMENT_RESULTS.md` and
`docs/RGFCA_ROI_QUALIFICATION_AUDIT.md` as measurement limitations.
P500 implementation qualification is documented in
`docs/P500_WHITE_CONTROL_IMPLEMENTATION_ERRATUM_20260915.md`.
The static gate does not authorize opening. No P500 biological results exist
in this package. CI success is not a measurement-accuracy result.

The [measurement-reference audit](FCP_H123_MEASUREMENT_REFERENCES.md) records the
primary sources supporting the main text's photography/calibration limitations.
Their validation does not transfer automatically to this estimator or cohort.

### Table S1. Completed Monarda operational measurement diagnostic

| Criterion | Observed | Fixed minimum | Decision |
| --- | ---: | ---: | --- |
| Pooled prediction precision | 0.56824250 | 0.70 | Failed |
| Pooled reference recall | 0.42608787 | 0.35 | Passed |
| Median annotated-image prediction precision | 0.51480059 | 0.70 | Failed |

All three conditions were required; the overall decision is **failed**.
The denominator is 109 annotated images, including 15 empty predictions, from
110 accounted-for images. One annotation-unknown image is not a verified
negative and does not enter overlap metrics. Pooled reference, prediction and
intersection counts are respectively 28,428,697, 21,316,820 and 12,113,123 pixels.
Pooled precision divides intersection by prediction; recall divides intersection
by reference. Per-image median precision retains empty predictions with zero
precision under the frozen convention. Larger regions contribute more weight
to pooled metrics than to the image median. These are descriptive operational
metrics, not independent pixel-level confidence estimates.

The original source files remain unchanged:

- [Saved result](../data/validation/monarda_region_agreement_v1/monarda_region_agreement_result_v1.json),
  exact SHA256 `49b2f017a6accb305867d62e7b4afe4e8ed96fcd03f8aafde86bda2da915fa4c`.
- [All-image count table](../data/validation/monarda_region_agreement_v1/monarda_region_agreement_rows_v1.csv),
  exact SHA256 `67cbeaa2cd0ea1f1d3f98516b4fe9f7bb1abc8648ac80268e6dc4ec3e754dc22`.
- [Original independent verification](../data/validation/monarda_region_agreement_v1/independent_verification.json).

Publication checks rederive Table S1 from retained integer pixel counts and
compare the saved result, exact byte identities, image census and table cells.
They do not verify annotation truth or rerun the estimator. Generic flower
polygons are not exhaustive focal-taxon petal labels. The observed disagreement
cannot separate localization error, incomplete annotation and target mismatch;
no such explanation is selected as a proven cause. Provider train/valid/test
labels are not FCP holdouts, and this opened export cannot serve as a new
untouched validation set. Neither these metrics nor the JRC box qualification
calibrates reflectance, ultraviolet signal or biological morph frequencies.

## S5. Package completion checklist

- [ ] All source-to-number mappings checked automatically.
- [ ] Dedicated main/supplement figures generated and visually inspected.
- [ ] Methods and bibliographic references audited together.
- [ ] Measurement interpretation and P500 outcome/STOP finalized.
- [ ] Main text and supplement validated as one submission package.

The retired 34-species comparison and six-species spatial analysis are excluded
from this supplement, while their original files remain unchanged in history.
