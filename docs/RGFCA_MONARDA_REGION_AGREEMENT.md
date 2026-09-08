# RGFCA Monarda flower-region agreement — qualified execution checkpoint

## Decision first

**Execution clarification:** the earlier qualification below did not resolve
the full runtime's internal colour computation or bind local execution to its
checked sources/environment. The [pre-outcome amendment](RGFCA_MONARDA_EXECUTION_AMENDMENT.md)
preserves this original contract, clarifies that scope and adds single-execution
guards. Actual image scoring remains closed pending a separate authorization
for the amended, qualified implementation. Do not run the historical CLI.

The limited Monarda flower-region agreement diagnostic is **prospectively frozen and pre-pixel qualified**, but the actual 110-image model/reference comparison has **not yet run**. The pinned user-supplied archive is intentionally not committed to this repository, and no Monarda image was decoded by the qualification workflow.

This checkpoint does not change any ecological result. In particular, the independent reserve matched flower-minus-background differential remains **not supported (p = 0.087)** and `flower_specific_robust_replication=false`.

## Frozen design

The prospective contract is:

- `docs/supporting/rgfca_monarda_region_agreement_contract_v1.json`
- specification commit `336b4b8052c54d396f24a5ccd26d0ca2fee3669f`

The denominator is fixed at all **110** images in the received Roboflow COCO-segmentation v1 export:

- 109 images have one or more generic `flowers` polygon annotations and contribute region-agreement metrics;
- one image has zero annotations and remains `reference_unknown_empty`; it is retained in the 110-image census but is **not** treated as a verified flower-negative.

The reference is the deterministic union of all COCO polygon components for the generic `flowers` category. It is not verified focal-taxon petal truth.

The FCP estimator is the unchanged ROI-v4 detector-plus-EfficientSAM runtime already frozen and used by the atlas. No Monarda training, fine-tuning, threshold search, morphology, component filtering, colour prior or provider-split parameter selection is allowed.

The limited operational gate is conjunctive:

- complete reference alignment and all 110 images accounted for;
- pooled prediction precision >= 0.70;
- pooled reference recall >= 0.35;
- median positive-image prediction precision >= 0.70.

These numerical floors are borrowed unchanged from the previously frozen JRC locked-test operational gates. They are used here only as bounded adequacy floors and are **not universal segmentation-accuracy standards**.

## Failure and unknown handling

Before any model execution, every archive member used by COCO must decode with stored JPEG dimensions exactly equal to the COCO declarations, and EXIF orientation must be absent or 1. Any decode, dimension, orientation, polygon-rasterization or archive-identity failure makes the complete diagnostic not evaluable rather than adapting the denominator.

For the 109 positive-reference images, an empty prediction or model runtime failure remains in the denominator with zero precision, recall, IoU and Dice. The zero-annotation image receives only descriptive model status and predicted-pixel count; it never becomes a verified false-positive example.

## Implementation and pre-pixel qualification

Implementation:

- runner commit `132726c891ea1062bf8943c1c531966169215cf6`
- focused test file added at `tests/test_rgfca_monarda_region_agreement.py`
- read-only qualification workflow `.github/workflows/rgfca-monarda-region-agreement.yml`

Two purely technical pre-outcome failures occurred and were recovered without changing the frozen design:

1. workflow expected a nonexistent nested `counts` key in the already committed source-mapping receipt;
2. a synthetic unit-test helper accidentally forced 100 predicted pixels, so its requested high-precision example was mathematically only 0.5/0.2 precision.

Neither failure decoded a Monarda image or loaded the model. The helper was corrected before any target-domain outcome.

Final pre-pixel qualification:

- commit `846b2ad81bad736b8f74675ab2f483488bb605f9`
- GitHub Actions run `34189809578`
- **8/8 focused tests passed**
- frozen ROI-v4 code blobs and detector-weight SHA matched
- Monarda reference pixels decoded: **false**
- model loaded/executed on Monarda: **false**

## Exact input required for the one permitted execution

The diagnostic runner requires the same received archive audited at intake:

`monarda_fistulosa_segmentation.v1i.coco-segmentation.zip`

SHA-256:

`7c9d213514e319b11cc15623747277f88eb37c453722f665b8570b40b9f0dd0d`

Size: 41,055,963 bytes.

The archive is deliberately not committed because the original photo licences vary and the repository does not publish the JPEGs or polygon vertices.

Once the exact archive is available in an execution environment, the runner must execute once under the frozen contract. All 110 row-level statuses and all aggregate metrics must be retained whether the limited gate passes or fails. A failure must not trigger model retuning or a second parameter search on these same images.

## Claim ceiling

A pass would support only bounded agreement of the unchanged FCP ROI-v4 localization with generic `flowers` polygons in this 110-image Monarda export. It would not establish focal-petal accuracy, annotation completeness, calibrated flower colour, taxonomic generality, ecological replication, causal biology or passage of the reserve flower-specificity gate.

A fail would be retained as measurement-validation evidence and would further limit interpretation of atlas colour results. It would not erase the already observed weak photo-level distance-colour associations, but it would make biological localization claims less defensible.
