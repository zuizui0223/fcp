# FCP: source-image availability and original-byte pilot, 2026-10-10

## Where this sits

The frozen *New Phytologist* paper and its evidence ledger show recurrent 50-km photographic colour-label geographic allocation in the discovery, validation and third high-depth cohorts. Constructive post-outcome label swaps in PR #150 erased the **mean** positive local-depletion statistic with 54/46/40 distinct strategically chosen edited photos, without altering any original source photo. The resulting 405-case blinded review queue is an actionable source-photo target panel, **not a human review result**.

This stage tests whether original image bytes are retrievable *today* on the iNaturalist open-licensed photo dataset and match the original measurement-time SHA256. It neither relabels photos nor quantifies empirical class-specific error rates.

## Frozen 24-photo pilot

- Balanced by the processor, using original source-record cohort and panel membership **only before** exposing blind review IDs: each discovery, validation and third source supplies 3 high-leverage, 3 same-species/same-original-colour matched controls, 2 photo-random controls.
- Reviewer/exportable input has ONLY `audit_case_id`, exact original `photo_id`, and original measurement-time `image_sha256`. No previous morph, panel membership, species, spatial position, observer identity or reviewed colour. Filename: `data/audit/fcp_blinded_source_photo_pilot_24_20261010.csv`.
- Input SHA256 `06191b36eca9862f6241b5110ae13aee0ff7e7414a5caaaa1025daffca33eadd`; 24 exact distinct original photo IDs.
- Public original-photo retrieval is restricted to paths inside `inaturalist-open-data.s3.amazonaws.com/photos/{original_photo_id}/`, in `large/original` sizes and `jpeg/jpg/png` variants. Never query arbitrary user-supplied URLs or switch to `static.inaturalist.org` which may contain non-open-licensed source images.
- Each retrieved image is decoded as an RGB photo; bytes are checked against the frozen original `image_sha256`. A valid photo ID with different bytes is **not** called historically byte-identical. A photograph unavailable on the open-data host remains `unavailable_on_open_data_host`, not a biologically non-flowering record.
- Only image-level technical quantities (dimensions, whole-frame over/under-exposure proxies and RGB channel spread) are calculated. These **cannot** identify petal colour or flower visibility in a non-segmented photo.
- Raw pixels are processed transiently and **not written, retained or uploaded to GitHub**. Outputs are a per-ID byte-provenance/technical metadata CSV and a cohort-independent numeric JSON.
- Unavailable image links may reflect a licence change or missing alternate file format; do not equate absence on the CC open-data bucket with actual deletion from iNaturalist.
- Exact iNaturalist license terms and author attribution must be respected if any future human-review package distributes images; no such package is constructed here.

## Readout rules

Full success of the CI's *technical* guards is distinct from scientific coverage. Only a count of actually SHA256-matched retrieved images demonstrates **original-byte** recovery. Different SHA256, even with identical photo ID, means the downloaded bytes differ from the measured source, potentially due to available size/rendering; this alone is NOT evidence of a flower colour change or image substitution.

The 24-photo preselected pilot is not a representative assessment of all 405 images, and original images may have new URLs/rights. No unblinded matching of adjudicated labels is undertaken. The original 42,111-species global opportunity frame, all previously frozen H1/H2/50-km scientific statements and unopened 2,000+730 future taxa remain unmodified.

The only practical next inferential advance after byte-level availability is **independent blinded human botanical review**, with reviewer agreement, floral-organ identification, source photo license/attribution, and colour-conditional error audit against the unblinding key.
