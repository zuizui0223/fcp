# FCP v2 measurement-validity terminal claim freeze — 2026-09-25

Status: **TERMINAL v2 INTERPRETATION FROZEN**

This document freezes the interpretation of FCP v2 after the one-shot fresh cohort, response-blind technical seal, source-byte-verified biological Pass B, and MV1–MV4 summaries completed.

It does **not** redefine H1–H3, add a new null family, refit q_white, or reopen the current New Phytologist manuscript beyond the explicitly limited D-transport paragraph described in section 8.

## 1. Chain of custody

Fresh cohort:

- Panel P: 200 previously studied species with fresh photo IDs;
- Panel N: 200 previously unused species;
- 100 rows/species;
- 400 species / 40,000 fresh rows.

Technical Pass T:

- 256/256 partitions completed;
- image acquisition = 40,000/40,000;
- technical metrics available = 34,616 rows;
- technical metrics unavailable = 5,384 rows;
- heavy counterfactual subset = 8,000 rows;
- biological outcomes remained sealed;
- technical table SHA256 = `f757c90eddcb00f43d504a2ae6ccaeff4639d7ce468f231bc3b75675d66440e9`.

Biological Pass B:

- 256/256 partitions completed;
- terminal rows = 40,000;
- biological measurement complete = 39,999;
- image acquisition failures = 1;
- source-byte drift = 0;
- replacements = 0;
- new null generated = false.

Canonical terminal artifact:

- workflow run = `35941316633`;
- artifact = `fcp-v2-biological-and-validity-result-v1`;
- artifact ID = `10836260454`;
- artifact digest = `sha256:8c83ca3d4b5d548c5ce345bb6108068e463a837258558a4286ddc9879b9e100d`.

## 2. Terminal synthesis

The best-supported effect-size description is:

> **The FCP image-measurement system is technically sensitive at the image level, especially to exposure and upstream ROI/background admission, but the species-level D phenotype and the recurrent q_white/W geometry are substantially more stable than individual image states.**

This is a descriptive synthesis, not a new omnibus pass/fail verdict.

## 3. Fresh-run D transport is strong

Among 136 Panel-P species evaluable in both the earlier frozen analysis and the fresh v2 run:

- Spearman rho = **0.9681064**;
- Lin CCC = **0.9720478**;
- calibration slope = **0.9690959**;
- calibration intercept = **0.0106821**;
- median absolute D change = **0.0147619**;
- mean signed fresh-minus-prior D change = **+0.0031680**.

This is stronger evidence than observer-disjoint splitting alone because the comparison crosses a fresh photo set and a separate measurement execution while keeping the D definition and measurement system fixed.

The defensible claim is:

> D behaves as a reproducible species-level phenotype within the same image-measurement system across independent fresh-image execution.

It is **not** independent-source replication because the source/opportunity universe and measurement system remain the same.

## 4. Exposure causes real image-level measurement movement

Fixed-mask exposure shifts produced:

- EV -0.5: complete-state disagreement 5.73%, morph flip 2.00%, mean palette Hellinger 0.104;
- EV +0.5: complete-state disagreement 6.15%, morph flip 2.38%, mean palette Hellinger 0.095;
- EV -1.0: complete-state disagreement 10.53%, morph flip 8.44%, mean palette Hellinger 0.227;
- EV +1.0: complete-state disagreement 10.65%, morph flip 8.02%, mean palette Hellinger 0.172.

Species-level D remains much more stable:

- EV -0.5: Spearman 0.944, CCC 0.944;
- EV +0.5: Spearman 0.921, CCC 0.902;
- EV -1.0: Spearman 0.821, CCC 0.801;
- EV +1.0: Spearman 0.758, CCC 0.709.

The union of the pre-frozen primary technical exclusions leaves D highly stable:

- Spearman = **0.985**;
- CCC = **0.991**;
- median absolute D change = **0.0103**.

## 5. Technical exposure preferentially moves colour along q_white

For fixed-mask same-image exposure perturbations, the technical colour displacement itself has substantial squared projection on the already frozen q_white direction:

- EV -0.5: mean T_white = **0.298**, mean white-fraction change = -0.0276;
- EV +0.5: mean T_white = **0.355**, mean white-fraction change = +0.0394;
- EV -1.0: mean T_white = **0.293**, mean white-fraction change = -0.0700;
- EV +1.0: mean T_white = **0.397**, mean white-fraction change = +0.0904.

The isotropic reference for a squared projection in the eight-dimensional zero-sum colour subspace is 0.125. No new v2 null distribution was generated.

Therefore:

> q_white is partly aligned with a preferred direction of the image-exposure measurement process.

This materially narrows H2 interpretation. q_white cannot be treated as a purely biological direction independent of image formation.

## 6. H2 geometry is exposure-sensitive but structurally persistent

Fresh-v2 baseline:

- vector species = **112**;
- observed W = **0.5454**.

Under fixed-mask exposure:

- EV -0.5: own-set W = 0.5293; common-species W change = **-0.0370**;
- EV +0.5: own-set W = 0.5651; common-species W change = **+0.0313**;
- EV -1.0: own-set W = 0.5163; common-species W change = **-0.0370**;
- EV +1.0: own-set W = 0.6060; common-species W change = **+0.0672**.

For ±0.5 EV, the median cosine of common species displacement vectors remains high:

- EV -0.5: **0.9813**;
- EV +0.5: **0.9904**.

Yet the lower tails include sign-reversed species, so not every species is stable.

The union of the primary technical-stratum exclusions changes W little:

- own-set W change = **-0.0073**;
- common-species W change = **-0.0119**;
- median displacement-vector cosine = **0.9991**.

Thus:

> The white-axis signal is partly measurement-dependent in magnitude, while much of the species-level displacement geometry persists.

The original structured-null median 0.457 remains historical provenance only. v2 does not create a new structured null or a new confirmatory H2 p-value.

## 7. Spatial organization is comparatively stable

The baseline fresh-v2 D–spatial descriptive association is rho = **0.0980**, close to the earlier reserve result.

For the species-specific spatial statistic itself:

- EV -0.5: Spearman agreement with baseline = 0.889, CCC = 0.909;
- EV +0.5: Spearman = 0.876, CCC = 0.892;
- EV ±1.0: Spearman ≈0.763–0.769;
- union of primary technical exclusions: Spearman = **0.952**, CCC = **0.955**.

No new spatial permutation null or p-value was generated in v2.

## 8. Scope split across papers

### Current New Phytologist manuscript

Only the fresh-run D transport result is imported because it directly strengthens the paper's H1 measurement-reproducibility claim without changing the H2 inferential spine.

Permitted imported statement:

> In a later fresh-image execution using the same frozen D definition and measurement system, 136 overlapping species retained strong D agreement with the earlier estimates (Spearman rho = 0.968, Lin CCC = 0.972, calibration slope = 0.969; median absolute change = 0.0148).

Required qualifier:

> This is fresh-image / fresh-execution transport within the same iNaturalist measurement system, not independent-source replication.

The current paper does **not** import:

- the 17-condition counterfactual battery;
- MV3a technical q_white displacement;
- MV3b W perturbation sensitivity;
- MV4 counterfactual spatial sensitivity;
- heavy-subset ROI/background results.

The existing direct-highlight limitation remains sufficient for the current paper's H2 boundary.

### Separate measurement-validity paper

FCP v2 itself is reserved as a separate measurement-validity study, potentially paired with disttrait as the generic implementation layer.

Its central story is:

1. image states are technically sensitive;
2. D is much more stable than individual image states;
3. exposure perturbations preferentially project onto q_white;
4. W changes quantitatively under exposure but remains structurally persistent;
5. spatial organization is comparatively stable;
6. the next unresolved step is independent measurement-system transport (FCP v3).

## 9. ROI-jitter interpretation boundary

ROI-jitter counterfactuals classify deterministic prompt masks without recomputing the baseline horizontal-flip admission gate. Their large apparent classifiability increase is therefore **not baseline-comparable**.

For ROI-jitter, retain only:

- palette displacement;
- Hellinger distance;
- q_white displacement direction.

Do not interpret jitter complete-state disagreement or classifiability change as direct biological admission improvement.

## 10. Final boundary

FCP v2 resolves the earlier question from:

> Is there measurement confounding?

to:

> How large is the measurement effect, which estimands are sensitive to it, and which species-level structures persist despite it?

The answer is now quantitative: image-level colour states are measurably exposure-sensitive, q_white is partly aligned with the exposure direction, but species-level D is strongly reproducible and much of the recurrent and spatial structure persists.
