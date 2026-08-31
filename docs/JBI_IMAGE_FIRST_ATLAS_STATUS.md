# FCP image-first atlas execution status

## Current decision

The active mainline is the image-first global flower-colour atlas. The former PR #21 target of 12 literature-selected species × 200 photographs is retained as a completed metadata-method precursor, not as the terminal cohort or confirmatory design.

Current gate: **metadata and geometry remain frozen; bulk atlas image opening is stopped pending estimator validation**.

A separate, preregistered location-blind image-measurement admissibility study is now complete. It used six development species and an independent locked partition for the three species that passed every image gate. All three locked spatial tests were non-significant under the frozen combined rule. The atlas candidate pixels themselves remain unopened, so no 50-species colour field, transition boundary or shared-boundary result exists.

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

## Completed automated-colour admissibility study

The frozen `fcp-inaturalist-automated-colour-state-v2` workflow used a pinned CLIPSeg revision, three positive and two negative prompts, continuous CIELAB features, prompt/reflection stability gates and a coordinate firewall.

- Development: 480 encounters and 886 photographs across six species; three species passed and three became `not_evaluable`.
- Locked: 360 encounters and 717 photographs across the three passing species.
- Cache audit: 717/717 valid records, with zero missing, partial or unexpected records.
- Locked admission: 101, 97 and 108 encounters for *Erythranthe lewisii*, *Hesperis matronalis* and *Orchis mascula*, respectively.
- Spatial confirmation: no species rejected the 9,999-permutation species-conditioned random-mark null after BH correction; primary q values were 0.59685, 0.9637 and 0.59685.
- Observer-removal and flower-minus-background requirements did not rescue any species.

The frozen decision is `spatial_organization_not_detected` for all three species. This does not prove spatial randomness and does not revise the six-species Chapter 1 result, which used different species, image measurement and statistic.

## Next allowed gate

Freeze and pass an independent flower-tissue localization benchmark and signal-recovery simulation for the exact atlas estimator. The benchmark must be scored without coordinates and must quantify localization completeness, contamination and attenuation under the already declared 100/250/500-km analysis scales. Its pass/fail rule must be committed before the atlas benchmark images are scored.

Until that gate passes, do not open the 20,200 atlas candidate images, reconstruct atlas colour fields or run shared-boundary concentration. The metadata freeze and the negative three-species validation are publishable results; neither is evidence for a 50-species colour boundary.
