# Monarda reference: complete archive intake and source-ID audit

Completed 8 September 2026 (Japan). The user supplied the finished COCO ZIP;
the earlier incomplete-download/access blocker is resolved. This is a real
image/annotation-pair census and source-ID audit, **not a model benchmark or a
new ecological result**. No image pixels were decoded or viewed, masks
rasterized, model run, colour extracted, or geographic fields joined.

## Received material

`monarda_fistulosa_segmentation.v1i.coco-segmentation.zip` is 41,055,963 bytes,
SHA-256 `7c9d213514e319b11cc15623747277f88eb37c453722f665b8570b40b9f0dd0d`.
All 115 members passed full read/CRC verification, totalling 41,434,856
uncompressed bytes. ZIP decompression and JPEG-byte hashing are not JPEG pixel
decoding. The local archive is not extracted or republished.

| Provider split | JPEG/COCO pairs | Annotations |
| --- | ---: | ---: |
| train | 77 | 537 |
| valid | 22 | 134 |
| test | 11 | 117 |
| Total | 110 | 788 |

All image members have exactly one COCO image row. There are no exact image-byte
duplicates or duplicate original filenames. COCO image IDs restart in each
split; joins must use `(split, image_id)`, not `image_id` alone. Each split
declares category IDs 0 and 1 with the same name `flowers`, but all annotations
use ID 1. All 788 polygon components pass the implemented finite-array and
declared-dimension bounds checks; this does not verify area, topology, image
alignment, botanical identity, annotation completeness or accuracy.

One image has **zero annotations**: train ID 69, original `22132.jpg`, mapped
photo ID `212509641`. Its status is unknown: retain it in the 110-image census,
do not declare it flower-absent, discard it silently, or score a predicted
flower there as a verified biological false positive.

Every COCO `date_captured` equals the export creation timestamp. These are not
verified observation dates and must not be used for season/phenology controls.
Image dimensions remain COCO declarations, not decoded-pixel or JPEG-header
measurements. The provider README reports auto-orientation with EXIF orientation
stripping and no augmentation. Its dummy `object-detection/undefined` information
URL is not a provenance identifier.

## Documentary source correspondence and overlap

The ZIP's dataset README names the
[original Monarda project](https://universe.roboflow.com/patricks-dashboard/monarda_fistulosa_segmentation),
not the receiving user's fork. Source correspondence uses the author's
[immutable repository](https://github.com/pmckenz1/monarda_fistulosa_color/tree/49812a431abce78c871f23aff6bc7978573dd502).
Notebook 1 writes photos as the zero-based pandas row index of `multimedia.txt`;
COCO `extra.name` preserves that numeric filename. Photo IDs come from the
multimedia `identifier` URL. Observation IDs are joined through `gbifID` to
`occurrence.txt` / `occurrenceID`;
the multimedia `references` field points to **photos**, not observations.

- **110/110** exported images map to 110 unique photo IDs and 110 unique
  observation IDs under that documented filename convention.
- Full immutable discovery and reserve frames each contain **50,000 rows**.
  Both photo-ID and observation-ID intersections are **zero in each frame**.
- Comparison reads exact Git blobs, not newline-translated working copies.
  The report pins the two commits, paths, full-input hashes and denominators.
- This is documentary filename/index linkage, not byte identity against
  reacquired original images or perceptual/event/observer independence. JRC,
  foundation-model training, legacy and other datasets are not thereby cleared.

The author metadata files are pinned by both Git blob and SHA-256. The mapping
retains photo/observation identifiers and historical source rights, not raw
occurrence tables, coordinates, dates, collectors or observer fields. The
report links the author notebooks and records their immutable blob identities.

Selection remains partly unresolved: notebook 3 samples 200 GPT flower-present
indices with `np.random.choice`, without an explicit seed or `replace=False`.
This does not reproduce or explain the released 110-image selection. Provider
train/valid/test labels concern its own model; they do not establish an untouched
FCP holdout. Do not tune on these images and then describe the same scores as
independent validation.

## Source rights and publication boundary

The source metadata archive records 97 photos as CC BY-NC 4.0, one as CC
BY-NC-SA 4.0, nine as CC BY 4.0 and three as CC0. Thus **98/110 have historical
noncommercial conditions**. These are source-record observations, not current
rights verification or a legal clearance. The export's blanket CC BY notice
does not establish that the original photos can be republished under that
licence. No ZIP, JPEG, original occurrence table or polygon vertices are
committed or uploaded by this audit. CI publishes only code, tests and the
bounded metadata/hash receipts.

## What this unlocks, and the next gate

The archive-access, complete pair-census, documentary photo correspondence and
two-frame ID-overlap gates are now complete. Generic `flowers` labels do not
establish focal-species petal truth. The next prospective gate must state the
narrow **flower-region agreement diagnostic**, retain the full 110-image
denominator and unknown empty-label case, and fix image orientation/alignment,
the already frozen FCP runtime, mask aggregation and scoring rules before any
pixel decoding/model execution. It must report every failed or unscorable image
and cannot relabel the reference or retune the model using its scores.

No benchmark pass threshold or measurement-accuracy result is asserted here.
An eventual flower-region score cannot establish calibrated reflectance,
pollinator colour vision, focal-petal accuracy, global taxonomic transfer or
flower-specific ecological replication. Any stronger ecological claim needs
genuinely fresh prospective validation; both existing 50,000-photo cohorts are
already opened. The evaluable reserve differential `p = 0.087` and discovery
background-recovery STOP remain unchanged.

## Reproduction and retained receipts

Both committed JSON files use exact UTF-8/LF bytes. They reproduce from the
local pinned archive and author metadata; the initial Windows CRLF receipts
remain local history and are not mixed into this canonical hash chain.

- [Complete archive intake](supporting/rgfca_monarda_archive_intake_v1.json):
  SHA-256 `1f5acb5af939b83c76b75d6513a0f404019ea88d4a4924d6d1322c8b8cde70f8`.
- [All-image source mapping](supporting/rgfca_monarda_source_mapping_v1.json):
  SHA-256 `fd0c58eb2657b8a7f5c620dca851e195ed341aa6f6781e7c704645b14a308106`.
- `scripts/analysis/audit_rgfca_monarda_archive.py` checks the exact received ZIP.
- `scripts/analysis/map_rgfca_monarda_reference_sources.py` checks pinned inputs,
  maps all 110 sources and compares both complete immutable FCP frames.
- `tests/test_rgfca_monarda_reference.py` reconstructs the saved receipts and
  checks artificial malformed archives, unsupported geometry, incomplete
  frames, source-ID parsing and the photo-versus-observation join distinction.

The reference workflow is read-only and runs offline tests only for this
checkpoint. It does not acquire photos, run a model or reopen the completed
six-document FlowerMask acquisition. Passing CI verifies the implementation
and retained evidence, not segmentation accuracy or ecological support.
