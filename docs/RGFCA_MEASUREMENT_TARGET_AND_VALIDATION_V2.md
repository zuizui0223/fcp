# RGFCA measurement target and independent validation: next design

Design recorded **8 September 2026, after the Monarda outcome**.
**No new model, threshold, acquisition cohort or benchmark execution is authorized
by this document.** This is the next measurement-development design, not a
retroactive amendment of any completed analysis.

## Decision and target

Continue the image-first, species-conditioned atlas, while keeping its measured
quantity explicit. The next engineering target is **visible flower-region
localization**, evaluated against a reference that actually annotates that
target. It is not renamed petal-pixel truth to gain a stronger interpretation.
Independent focal-taxon attribution and colour validity remain necessary for
interpreting photo-labelled measurements as species' floral phenotypes.

| Layer | Quantity / unit | What reference evidence could support it |
|---|---|---|
| Frozen atlas measurement | A labelled photograph's apparent RGB-derived colour from all retained model flower regions | Existing result and exact implementation; no anatomical or focal-taxon certification implied |
| Next bounded localization diagnostic | Visible generic flower foreground, with source-defined unknown/ambiguous tissue preserved | Paired flower-region annotations with documented organ inclusion, alignment, completeness and split provenance |
| Stronger species-trait interpretation | The focal observation taxon's floral display, including an explicit organ definition | Independently adjudicated taxon-to-region correspondence; a generic flower mask alone cannot supply it |
| Colour validity / ecology | Apparent photographic colour versus calibrated reflectance; species-conditioned spatial or environmental effects | Separate colour-error evidence and prospectively fixed independent ecological validation; overlap scores do not validate either |

No RGB-photo workflow here is claimed to recover reflectance or pollinator
perception. Continuous photographic colour can still be the stated endpoint;
its camera, illumination and tissue-selection limitations must remain explicit.
The species-free map/photo bar is presentation. Conditioning inference on the
photo's taxon label does not itself validate which taxon's tissue was measured.

## Completed evidence remains closed

- Monarda: all 110 images completed, 109 annotated plus one unknown, zero model
  runtime failures. The original three-part agreement gate did not pass.
  See [the complete result](RGFCA_MONARDA_REGION_AGREEMENT_RESULTS.md).
  No re-scoring, successful-subset analysis, changed floor or alternate model
  selection on these images is part of this route.
- The two 50,000-photo atlas cohorts are outcome-opened. The weak reserve
  association replicated, but its matched differential remained unsupported
  (`p = 0.087`). Repartitioning those photos cannot create a new holdout.
- JRC's box test, FlowerMask's six schema inspections and sharedness-v2 geometry/
  synthetic work retain their distinct evidential roles. None establishes
  target-domain focal-petal accuracy. Legacy six/34-species results stay frozen.

## Reference decisions, before choosing a new model

The [three-family primary-source audit](RGFCA_INDEPENDENT_REGION_REFERENCE_AUDIT_20260908.md)
finds no ready joint focal-petal/iNaturalist reference. Its narrower decisions
are not a requirement to prove impossible universal independence:

1. **FlowerMask:** retain as a candidate for a documented generic-flower target.
   Use the existing 300 exposed annotation IDs as a bounded source frame, not
   a claim of complete access to all 3,600 advertised records. A complete public
   paired manifest and annotation-scope evidence are still needed.
2. **Oxford:** exclude 102 algorithmic segmentations as human pixel truth.
   Oxford 17's manual trimaps may support only their labelled region; unknown
   pixels must not become background. Current mask membership and rights remain
   unverified, as do event/observer groups.
3. **USDA:** retain an orchard, human-guided/algorithm-refined reference lead.
   Its public train/validation lists have 100/30 rows but share `IMG_0339.JPG`:
   exact union 129. This is a current filename conflict, not proof of historical
   model leakage. Preserve it; do not silently remove the row or invent a split.

### Reproducible USDA metadata check

The [four-file retained bundle](../data/validation/usda_flower_split_audit_v1/)
contains the two exact provider text files, a source manifest and the complete
audit. Their original mixed line endings are retained. The script verifies
provider MD5, locally computed SHA-256 and byte counts before parsing rows;
it reports exact duplicates and casefold collisions separately, without changing
the source keys. The single overlap affects 1/30 listed validation filenames
(3.33%), not an estimated model-bias magnitude.

```text
python scripts/analysis/audit_rgfca_reference_splits.py data/validation/usda_flower_split_audit_v1 --verify data/validation/usda_flower_split_audit_v1/result.json
```

The CI exit code certifies faithful reconstruction, **not** disjointness:
`filename_integrity_gate_passed=false` and `benchmark_execution_authorized=false`
remain required in this retained result. A filename-only pass elsewhere would
still not prove image-content, plant, observer, event or pretraining independence.
This script does not request network resources, read images, load models or
join coordinates. The audit is retrospective source checking, not preregistered
scientific inference. The cause and historical relevance of the source conflict
are unknown; no author contact was made.

## Next executable gate and stop rules

The next useful action is **public metadata-only paired-manifest feasibility
for the 300 already exposed FlowerMask annotation IDs**. Before requesting
further metadata, freeze the complete existing ID frame and a bounded request
plan, retain all 300 rows, and retain the six already schema-opened IDs as such.
Inspect provider file names/IDs, sizes, checksums and declared matching paths;
do not acquire another annotation body with embedded image data or image/mask
payload in this gate. Do not replace inaccessible IDs or rank records by model
performance, colour or geographic outcomes. No author contact or rights change
is authorized. Public metadata can establish documentary correspondence only,
not actual image/mask alignment or anatomical truth.

If this gate cannot establish the paired frame or annotation target, preserve
that specific unresolved outcome; do not proceed to a convenient subset score.
A later, separately specified inspection may be needed for pairing/alignment
and duplicate-content evidence. Do not restart the closed six-document job.

Before any later measurement-development benchmark:

- Record source version, source IDs, image/annotation hashes, rights scope,
  target organ rules, image transforms, unknown-label semantics, missingness
  and completeness. Unannotated does not mean biologically absent.
- Declare development and untouched validation membership using documented
  independent source/event groups where available. Audit known overlaps with
  JRC, Monarda and all opened FCP frames. Unknown foundation-training or observer
  overlap is reported as unknown and limits the claim; do not mark it verified
  merely because filenames differ. If groups are unavailable, do not claim an
  event-independent holdout or population-level uncertainty from image counts.
- Freeze the model, error budget, denominator, failure handling and descriptive
  strata before scoring. Any operational floor needs a rationale for the new
  target and cannot be selected to outperform the failed Monarda gate. A
  separate development revision must name its data and validation route.
- Keep all admitted rows and technical failures in a durable terminal ledger.
  Keep ecological coordinates closed during measurement validation. A localized
  reference can support a localized methods result, not a global coverage claim.

Reference-access progress does not complete the biological publication goal.
The publishable evidence already in hand includes a bounded photo-association
result and an explicit measurement limitation; stronger floral or mechanistic
claims require their own independent evidence.
