# FCP image-first atlas execution status

## Current decision

The active mainline is the image-first global flower-colour atlas. The former PR #21 target of 12 literature-selected species × 200 photographs is retained as a completed metadata-method precursor, not as the terminal cohort or confirmatory design.

Current gate: **metadata and geometry frozen before image pixels; flower-ROI measurement contract is next**.

No atlas candidate image pixel, flower ROI or continuous colour value has been opened. Therefore no atlas colour field, transition boundary or shared-boundary result exists yet.

## Completed and immutable

The six-species Chapter 1 result is closed:

- 1,200 photographs, 480 calibration and 720 held-out evaluation;
- Stage A primary lower-tail `p = 0.0113`;
- Stage B primary upper-tail `p = 0.0906`;
- no retuning or reinterpretation in the atlas pilot.

The 34-species comparative paper is also retained intact. Its literature extraction and spatial classification remain useful method/evidence infrastructure, but those labels do not select the atlas cohort.

## Implemented for the atlas pilot

1. 50-species metadata-only cohort interface with 300/400/500-photo tiers.
2. Positional accuracy `<= 5 km`, public coordinates and licensed-photo admission.
3. 0.25° primary and 0.5° sensitivity thinning.
4. Observer, cell, month, spatial-coverage and seasonal-coverage gates.
5. Geometry-only 100/250/500-km scale selection, finest passing scale first.
6. Species-free display-table interface with species-conditioned inference retained.
7. GitHub Actions contract tests, live feasibility, artifact upload, validation and passing-freeze commit.

## Completed live metadata gate

GitHub Actions run `33355972418` completed successfully on 2026-08-31. Under the unchanged pre-image rules it:

- received 200 species-count records, retained 199 after static filters and audited 61 species;
- admitted exactly 50 species and selected 20,200 observations;
- assigned 19 species to the 500-photo tier, 14 to the 400-photo tier and 17 to the 300-photo tier;
- passed the geometry criteria at 100, 250 and 500 km and selected the finest passing scale, 100 km, as primary;
- retained 250 and 500 km as mandatory sensitivities;
- kept candidate image pixels, flower ROI and continuous colour closed throughout.

At 100 km, all 50 species were geometry-evaluable, all 50 had shared-opportunity support, and 418 cells had opportunity from at least three species. The frozen outputs are the 50-species cohort, the 20,200-observation manifest, the feasibility audit, all-scale geometry diagnostics and their manifest. Text hashes use the declared `sha256_lf_canonical_v1` mode so the same freeze verifies on Linux and Windows without changing scientific content.

Focused local validation: `tests/test_image_first_atlas.py` — **8 passed**. Atlas plus retained scale-up/manuscript regression selection — **40 passed**. The committed complete freeze also passes its independent validator on Windows.

The repository-wide pytest collection additionally requires the pre-existing optional Florence stack (`transformers==5.16.1`); it is validated in its dedicated workflow and is not imported by the metadata gate.

## Next allowed gate

Write and freeze the flower-ROI measurement contract: model revision, input image derivative, ROI acceptance/completeness rules, colour calibration and failure states. Do not open the 20,200 candidate image pixels until that contract and its validator pass.

The present freeze establishes metadata and geometry feasibility only. It is not evidence for a flower-colour field, a transition boundary or cross-species shared-boundary concentration.
