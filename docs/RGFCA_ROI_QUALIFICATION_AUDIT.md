# RGFCA flower-region qualification: verified scope and transfer limits

Audit: 7 September 2026. **Share with caveats as completed estimator
qualification; not independent ecological validation or submission clearance.**
This note reconstructs the saved JRC qualification, without training, new model
predictions, image downloads or reserve outcomes. It changes no frozen rule or
result. The related [training-source audit](RGFCA_TRAINING_SOURCE_AUDIT.md)
addresses data provenance and published sampling; the
[provider audit](RGFCA_MEASUREMENT_PROVIDER_AUDIT.md) identifies the models.

## 1. What was actually tested

The fixed YOLO11n detector was fine-tuned on all 400 public JRC development
images, using 6,991 positive-area clipped boxes. One source box with zero image
intersection remains recorded, not silently repaired. Training used 50 epochs
and the last epoch, without validation-based epoch selection. EfficientSAM-Ti
converted predicted boxes into candidate masks. The predeclared admission rule
required a retained instance, at least 100 flower-mask and 100 background pixels,
horizontal-flip mask IoU at least 0.5 and colour difference at most 5.0.
These are model stability rules, not reference flower-colour accuracy measures.
[Frozen estimator contract](https://github.com/zuizui0223/fcp/blob/2449c77597cc4c57f18eb0aaa446211143ce3be5/docs/supporting/jbi_atlas_roi_estimator_contract_v4.json).

The locked test contains 100 images and 2,524 evaluable reference flower boxes.
All 100 rows, including 15 admission failures, enter the saved aggregate
qualification. This is not a summary restricted to the successful 85 images.
JRC annotations supply flower bounding boxes, not reference petal masks or
spectral reflectance. The public train directory pools the source study's
training and validation slices; FCP does not count its training-domain check
as independent validation. Source-study sampling and geographical scope are
documented in the [training-source audit](RGFCA_TRAINING_SOURCE_AUDIT.md).

## 2. Independent arithmetic reconstruction of the saved test

The following values were recomputed from the 100 committed rows, not obtained
by a new model run. All eight original `>=` qualification decisions reproduce.

| Quantity | Exact denominator or definition | Reconstructed value | Frozen minimum |
| --- | --- | ---: | ---: |
| Images | All test images | 100 | 100 |
| Admitted images | 85/100 | 0.850000 | 0.80 |
| Detector precision at box IoU 0.5 | 2,008 matches / 2,749 predictions | 0.730447 | 0.70 |
| Detector recall at box IoU 0.5 | 2,008 matches / 2,524 reference boxes | 0.795563 | 0.35 |
| Medium-object recall | 112/135 reference boxes | 0.829630 | 0.35 |
| Large-object recall | 2/4 reference boxes | 0.500000 | 0.50 |
| Pooled mask containment in reference boxes | 5,543,384 / 6,448,251 predicted mask pixels | 0.859672 | 0.70 |
| Median image mask containment in reference boxes | Median of 100 per-image containment fractions | 0.913904 | 0.70 |

There are 741 unmatched detector predictions and 516 missed reference boxes.
The 15 admission failures remain present with their original flip-instability
reasons. The 2,385 small, 135 medium and four large reference boxes sum to the
complete 2,524-box denominator. **The large-object gate passes exactly at its
threshold, but its 2/4 denominator spans only three images.** It cannot support
a precise or broad claim about large-flower performance. We add no nominal
box-independent confidence interval: flowers are nested within images and
source survey points.

The size-bin definition was supplied by an explicit amendment during the first
training epoch, before completed weights, development predictions or locked
test outcomes. It uses clipped reference-box areas on a 512-equivalent canvas:
small <1,024, medium 1,024–<9,216, large >=9,216 pixels. This timing is preserved;
the definition is not misdescribed as present in the initial contract.
[Immutable amendment](https://github.com/zuizui0223/fcp/blob/2449c77597cc4c57f18eb0aaa446211143ce3be5/docs/supporting/jbi_atlas_roi_v4_reference_size_amendment_v1.json).

## 3. What these numbers do not validate

1. **Box containment is not petal segmentation IoU.** Its numerator is predicted
   mask pixels inside the union of annotated boxes; its denominator is all
   predicted mask pixels. Background or non-petal tissue inside a reference box
   can count as contained. Missing true petal pixels are not measured by this
   fraction. No reference-mask union appears in its denominator.
2. **JRC-domain qualification is not global taxon-uniform validation.** European
   grassland survey patches differ from worldwide iNaturalist photographs in
   composition, flower form, scale, camera and observer behaviour. The locked
   test provides a useful localization gate, not an error rate for every atlas
   species, continent or colour. The 85% JRC admission fraction is not a forecast
   of global atlas classifiability. The source's flower-present selection also
   does not supply flower-absent specificity validation.
3. **Flip consistency is repeatability, not colour calibration.** A stable wrong
   mask or camera colour cast can pass it. Ordinary RGB measurements do not
   establish spectral reflectance, ultraviolet traits, pigments or pollinator
   perception. Source-paper Faster R-CNN/Pl@ntNet benchmarks are not the scores
   of this different YOLO11n/EfficientSAM composite.
4. **A reproduced photo bar is not an independent annotation benchmark.** Its
   24 CC0 crops reproduce retained discovery measurement summaries. Selection
   and display do not create biological ground truth. Exact source-byte and
   palette-count replay is distinct from anatomical mask accuracy.
5. **Reserve/background controls are necessary but not universal measurement
   validation.** Their completion can test the prospectively stated transfer
   and image-background contrasts. They do not prove the model uniformly measures
   petals or eliminate every systematic error shared by regions and taxa.

These caveats limit the inference; they do not retrospectively change the
passing qualification, reject selected images, or authorize threshold changes.
The older SegFormer ROI v3 development STOP (17/400 admitted; recall
0.1207266485) is preserved separately. No new model search is opened by this audit.

## 4. Immutable evidence and reproducible checks

The complete qualification was committed at
`2449c77597cc4c57f18eb0aaa446211143ce3be5`. Its recorded executable revision was
`5d1c033f84b00a2bbc4276a33c4f4f4b82ec5b3d`. The later RGFCA runtime source,
`9fae6ccdf684a46026f72ba12e98de2c5c54bf2a`, contains the identical locked-result
bytes. Qualification rows need not be recreated on the current branch: tests
read their literal committed Git objects, not a result-directory search.

| Evidence at the qualification commit | Exact SHA-256 |
| --- | --- |
| [100 measurement rows](https://github.com/zuizui0223/fcp/blob/2449c77597cc4c57f18eb0aaa446211143ce3be5/data/atlas/qualification/roi_v4_locked_test/jrc_roi_v4_locked_test_rows.csv) | `9f834de4bcdd530581d7713376cc918aa587e80726bba32021fa2dfa1034b248` |
| [Locked result](https://github.com/zuizui0223/fcp/blob/2449c77597cc4c57f18eb0aaa446211143ce3be5/data/atlas/qualification/roi_v4_locked_test/jrc_roi_v4_locked_test_result.json) | `420ae571e0950bfe252d5f0de34b4e047785c3efb2968d9f14abfc2161c08aea` |
| [Evidence manifest](https://github.com/zuizui0223/fcp/blob/2449c77597cc4c57f18eb0aaa446211143ce3be5/data/atlas/qualification/roi_v4_locked_test/gate_evidence_manifest.json) | `dc467f9ac4e12325c86c5c4e305216b02adc6cc97b9fbebb67fd28236fbb89e8` |
| [Frozen gate executable](https://github.com/zuizui0223/fcp/blob/2449c77597cc4c57f18eb0aaa446211143ce3be5/data/atlas/qualification/roi_v4_locked_test/gate_executable.py) | `f43048839f5dde8ec1acd56c1163d988a8a7cddb024a289a6a2ba96b76c84074` |

The [eight publication audit tests](../tests/test_rgfca_roi_qualification_audit.py)
check seven exact source identities, runtime/result identity, the complete image
and object census, all saved metric arithmetic, all eight frozen decisions,
image admission, reference-size denominators and the retained historical STOP.
They do not independently reproduce model predictions, annotations, original
training or the source point split. They require the repository's full Git
history, as provided by the [publication workflow](../.github/workflows/rgfca-publication-figures.yml).

This source/arithmetic audit closes a documentation and provenance gap. The
scientific outcome remains a small exploratory image-derived association;
independent replication, flower-specific robustness and the final submission
package remain unfinished.
