# FCP image-first atlas execution status

## Current decision

The active mainline is the image-first global flower-colour atlas. The former PR #21 target of 12 literature-selected species × 200 photographs is retained as a completed metadata-method precursor, not as the terminal cohort or confirmatory design.

Current gate: **pre-image contract implemented; live 50-species metadata feasibility pending GitHub Actions**.

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

Focused local validation: `tests/test_image_first_atlas.py` — **7 passed**. Atlas plus retained scale-up/manuscript regression selection — **39 passed**.

The repository-wide pytest collection additionally requires the pre-existing optional Florence stack (`transformers==5.16.1`); it is validated in its dedicated workflow and is not imported by the metadata gate.

## Next allowed gate

Run `.github/workflows/jbi-image-first-atlas-metadata.yml` on the PR #21 head. The workflow must either:

- freeze exactly 50 admitted species, 300–500 observations per species and one geometry-selected primary scale; or
- stop as `not_evaluable` without downloading images or changing the rules.

Only a passing, hash-validated metadata and geometry freeze authorizes a separately frozen flower-ROI measurement phase.
