# Technical photo-byte preflight for FCP Validation 115 fixed photo IDs

**Status:** post-outcome technical check of images selected by a frozen metadata-only process; **not** a colour, pigment, genotype, floral organ or fitness estimate.

Original 10-km high-depth historical source-only climate/site-year opportunity remains HOLD. The previously source-verified metadata queue contains 25 species and 115 distinct iNaturalist photo IDs: 15 gap-one species (39 photos), 10 gap-two species (76 photos), with each year of the missing 2-slot species represented.

## Imaging contract

- Input is the GitHub-versioned `candidate_photo_quality_review.csv`, Git blob `467f1a1f51a8e14f243d25e7f62550521bee0b40`. Every row is UNREVIEWED for flower visibility, species identity, organ, exposure, colour and plant individuality.
- Fetch **only** the 115 numeric, previously source-frozen photo IDs, from the licensed iNaturalist Open Dataset. Try `medium.jpg` and then `medium.jpeg`, 500 px maximum side. This JPEG resolution is sufficient for initial technical accessibility, **not** all floral-quality questions. Do not query new observations or substitute different photo IDs. No colour-based selection.
- Enforce CC0/approved Creative Commons licences from the metadata source, numeric IDs, a fixed URL host, max response bytes and per-request timeout. Each in-memory image has its dimensions, SHA256, whole-image near-highlight/near-white and darkness fractions, plus a simple edge-energy *diagnostic* recorded. These **are not flower-mask metrics**. Technical failures remain explicit non-retrieval or decode errors.
- Do not upload, commit or redistribute photo pixels. Store numerical technical diagnostics, links to the original photo/observation page, source licence, photo identity and untouched botanical-review fields only. Photographer attribution and present-day licence must be revalidated before any image redistribution outside original source.
- The photo technical worker may report 0–115 image decode successes, but **all 115 remain botanically UNREVIEWED**. No photo earns botanical or colour eligibility automatically from resolution, brightness or edge metrics. Do not infer adaptation.
- Once technical outputs are fixed, the 115 photos require actual focal-flower/organ/exposure review. The 25 species qualify only when every deficient year has its original number of independent *usable* photos. No post-result swapping, and no use of the untouched new 2,000+730 species.

## Reproducibility

`scripts/analysis/audit_fcp_validation_photo_pixels_technical_20261009.py` runs with only numpy and Pillow. Unit tests use generated synthetic image bytes; no external URLs are fetched during testing.

GitHub Actions `.github/workflows/fcp-validation-photo-pixel-technical-20261009.yml` runs the full 115-photo public licensed-source technical inspection exactly once per new workflow version and commits only the derived output. It checks the frozen queue hash, all 115 identities and an intact HOLD claim. Original image pixels exist only transiently in runner memory, not in version control or uploaded Actions artifacts. External photo links should not be mistaken for user-owned redistributable images.

The Open Dataset describes file URL structure and photograph licensing at https://github.com/inaturalist/inaturalist-open-data and the official observation API documents sizes/domains at https://github.com/inaturalist/iNaturalistAPI .
