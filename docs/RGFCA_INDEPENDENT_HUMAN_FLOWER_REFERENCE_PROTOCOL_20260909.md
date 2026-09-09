# RGFCA independent human-adjudicated flower reference — prospective measurement gate

Date frozen: 2026-09-09 JST

This protocol defines the next measurement-validity gate required before the independent local-signal replication in `RGFCA_INDEPENDENT_LOCAL_SIGNAL_REPLICATION_DESIGN_20260909.md` may open ecological colour outcomes. It does not authorize model tuning, ecological inference, or reuse of outcome-opened validation frames.

## Target

The validation target is the **visible floral display of the focal observation taxon** in an iNaturalist photograph.

Human reference pixels use four mutually exclusive labels:

1. `focal_flower`: visible flower/display tissue belonging to the observation taxon;
2. `nonfocal_flower`: visible flower/display tissue belonging to another taxon;
3. `ambiguous`: tissue whose focal/nonfocal or flower/background status cannot be adjudicated from the image;
4. `background`: all remaining scored pixels.

This target is deliberately stronger than generic flower foreground. Nonfocal flowers are false positives for the focal target. Ambiguous pixels remain unknown and are excluded from pixel numerators and denominators rather than relabelled background.

No claim of petal-only anatomy, calibrated reflectance or pollinator vision follows from this target.

## Independent validation frame

Freeze a metadata-only manifest before any validation image pixels or model predictions are opened.

The untouched validation frame contains exactly **300 photographs from 50 species, 6 photographs per species**.

For every species:

- photographs must come from at least 3 distinct observers;
- no observer may contribute more than 2 photographs for that species;
- all 6 `photo_id` values must be distinct;
- all 6 `observation_id` values must be distinct.

Across the full frame, every `photo_id`, `observation_id` and exact `observer_id` must be disjoint from:

- the 50,000-photo discovery and reserve atlas frames;
- the 110-image Monarda validation frame;
- JRC flower-detector development/qualification frames;
- the six previously schema-opened FlowerMask examples.

Known overlap is disqualifying. Unknown foundation-model pretraining overlap is reported as unknown rather than encoded as independent.

Species are selected from metadata only. No segmentation, palette, M2 score, candidate-cell outcome, island status or environmental response is available to selection. The 50-species frame is not enriched for the Canadian Rockies candidate; this is a measurement-transfer gate, not ecological confirmation.

If a metadata-only candidate pool cannot supply the complete 50 x 6 frame under these rules, record `human_reference_frame_qualified=false`; do not reduce the species count, photos per species or observer requirement after inspecting pixels.

## Annotation procedure

Before validation annotation, annotators may train on a separate set of non-validation example images that can never enter the 300-image frame. Model predictions remain hidden.

Each of the 300 validation images is labelled independently by two human annotators using the fixed four-class ontology above. Annotators are blinded to:

- model boxes/masks;
- colour palette measurements;
- M1-M3 scores;
- candidate-cell status;
- island/mainland classification;
- ecological hypotheses.

After both independent annotations are frozen, disagreements are presented to a third adjudicator who sees the image and the two human labels but still not the model prediction or ecological metadata. The adjudicated four-class mask becomes the sole scoring reference.

Retain the two original annotations, disagreement mask and final adjudication. Do not overwrite disagreement history.

## Frozen model and execution

Score the exact already-frozen ROI-v4 estimator identity used by the atlas. Before execution verify the detector and EfficientSAM model hashes against the existing contract. No threshold, detector, prompt rule, mask postprocessing, image transform or failure handling may change after the human reference is opened.

No development or tuning occurs on the 300 validation images. A failed gate closes this model/target pairing on this frame.

## Scoring semantics

For focal-taxon validation:

- reference positive = `focal_flower`;
- reference negative = `nonfocal_flower` + `background`;
- `ambiguous` pixels are excluded from both numerator and denominator;
- model-predicted pixels on `nonfocal_flower` count as false positives;
- model-predicted pixels on `background` count as false positives;
- unpredicted `focal_flower` pixels count as false negatives.

Report pooled TP/FP/FN and per-image precision, recall, IoU and Dice. Empty model predictions remain in the image-level denominator under the same convention used by the closed Monarda gate; they are not dropped as technical negatives. Technical failures remain in a terminal ledger and are not silently excluded.

Also report, as diagnostics rather than substitute endpoints:

- generic-flower agreement, treating `focal_flower + nonfocal_flower` as positive;
- fraction of predicted scored pixels landing on `nonfocal_flower`;
- pre-adjudication human-human agreement;
- species-equal distributions of focal precision and recall.

## Pass/fail gate

To avoid outcome-tuned standards, retain the same three operational floors used prospectively for the closed Monarda region-agreement gate:

1. pooled focal-target precision >= **0.70**;
2. pooled focal-target recall >= **0.35**;
3. median annotated-image focal-target precision >= **0.70**.

The measurement gate passes only if all three conditions pass on the complete 300-image frame. No alternate successful subset, taxon subset, observer subset, threshold or relabelling can replace the complete-frame result.

These are operational continuity floors, not universal botanical standards. Passing supports use of the stated visible-focal-flower photographic target under this pipeline; it does not validate calibrated colour or pollinator perception.

## Qualification and opening order

The required order is:

1. freeze and hash the 300-row metadata manifest;
2. verify all identity/disjointness/support rules without opening validation image pixels;
3. freeze annotation codebook and annotator roles;
4. acquire/open the exact 300 validation images and complete two independent annotations plus adjudication;
5. hash and freeze all human-reference artifacts;
6. verify the unchanged ROI-v4 model/runtime identities;
7. execute exactly one primary model-scoring run;
8. retain all 300 results and the terminal ledger;
9. independently recompute the three gates from retained pixel counts.

Ecological replication images and colour outcomes remain closed until step 9 reports `measurement_gate_passed=true`.

## Stop rules

If the human-reference frame cannot be fully qualified, stop before pixels.

If the measurement gate fails, do not:

- tune ROI-v4 on these 300 images;
- lower the three floors;
- relabel nonfocal flowers as acceptable focal signal;
- drop low-performing species/images;
- reopen Monarda, FlowerMask or another convenient reference to choose a passing model;
- open the ecological replication outcome.

A later model revision would require a separately named development dataset and a genuinely new untouched validation frame.
