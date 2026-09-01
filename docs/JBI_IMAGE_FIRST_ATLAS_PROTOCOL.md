# FCP image-first global flower-colour atlas — prospective protocol

## Mainline and frozen boundaries

The active FCP mainline is now an **image-first global flower-colour atlas**. The ordered analysis is:

`iNaturalist image -> automated flower ROI -> continuous colour -> within-species spatial field -> species-conditioned transition boundary -> prequalified cross-species boundary concordance`

Species identity is hidden only in the atlas display. It remains mandatory in every inferential graph, standardization step and permutation null.

Two completed analyses are retained without reopening them:

1. the six-species, 1,200-photograph Chapter 1 development result (Stage A `p = 0.0113`; Stage B `p = 0.0906`);
2. the literature-derived 34-species comparative analysis, retained as a classification method and supplementary comparison rather than the atlas admission frame.

## Gate 0 — pre-image contract

The v1 contract is frozen at `docs/supporting/jbi_image_first_atlas_contract_v1.json`. It governs the unopened 50-species sentinel cohort and, before any candidate image pixel is opened, fixes:

- the iNaturalist metadata root and flowering annotation;
- the 50-species sentinel size;
- 300/400/500-photo sample tiers;
- licensing and positional-accuracy rules;
- observer, spatial-cell and month caps;
- 0.25° primary and 0.5° sensitivity thinning;
- geometry-only candidate scales of 100, 250 and 500 km;
- season/phenology sensitivities;
- display and inferential separation;
- fail-closed stopping rules.

Literature class, flower colour, Stage A/B effects, environmental layers and candidate image pixels are prohibited admission inputs.

The later v2 expansion contract prospectively supersedes 50 species as the terminal experiment. It fixes eight disjoint cohorts of 25 species x 300 observations: 200 species and 60,000 observations in total. The v1 sentinel remains frozen and unopened; it is not pooled with, substituted into or used to retune v2.

## Gate 1 — metadata feasibility and cohort admission

Both the sentinel and terminal experiment begin from research-grade, flowering-annotated iNaturalist records under the frozen Angiospermae root (`taxon_id = 47125`). Records require a licensed photograph, public coordinates, a numeric positional accuracy no greater than 5 km, a date and a non-captive status.

For the v1 sentinel, species were ranked by the count returned by that metadata query and the deterministic selector admitted the first 50 species passing all rules. For terminal v2, the complete 500-species candidate pool is audited before a stable SHA-256 permutation selects eight consecutive, species-disjoint panels. At most one species per iNaturalist genus is admitted. The completed six development species are excluded so their frozen evaluation is not recycled as a new confirmatory cohort.

Within species, metadata rows are balanced without opening photographs:

- one licensed photograph per observation;
- 0.25° spatial round-robin selection;
- simultaneous cap of 5 records per 0.25° cell;
- cap of 10 records per 0.5° cell;
- cap of 10 records per observer;
- cap of 100 records per calendar month;
- minimum 50 observers, 60 primary cells, 40 sensitivity cells, four months and three hemisphere-adjusted local-solar quarters.

The v1 sentinel retained the largest passing tier in the fixed order 500, 400, then 300 photographs. Terminal v2 instead requires exactly 300 observations for each of 200 genus-distinct species; all eight panels are required. A failure is recorded without weakening gates or using colour to choose replacements.

The completed live feasibility audit in GitHub Actions run `33405153936` examined all 500 predeclared candidates. It found 358 geometry-eligible species and froze 200 genus-distinct species, 60,000 unique observations and 60,000 unique photos across the eight panels. Maximum retained counts were 10 observations per species-observer, 5 per species x 0.25-degree cell and 10 per species x 0.5-degree cell. This is selection and feasibility evidence only; the dated-source and ROI gates still control image access.

The read-only API supplies flowering-annotation eligibility and the exact bounded pre-colour selection. The official monthly Open Data schema does not include annotation term/value records, so it cannot independently reconstruct that eligible universe. Every exact selected photo, observation, observer and taxon must therefore reconcile to the fixed 2026-08-27 Open Data snapshot before acquisition. The first resolver, `jbi_atlas_dated_source_amendment_v1.json`, stopped `not_evaluable` when selected `photo_id = 950871` appeared in multiple association rows. Independent iNaturalist documentation establishes that photos and observations are many-to-many and that `photo_uuid + observation_uuid` is the association key. Before any selected association row or image was inspected, `jbi_atlas_dated_source_m2m_amendment_v2.json` froze a successor rule: retain every link for each selected photo asset and require exactly one association whose observer, taxon, date, quality, licence, coordinates and positional accuracy all match the already frozen API row. Zero matches, multiple matches, duplicate composite keys or conflicting asset fields stop without replacement. The v1 STOP remains part of the evidence chain.

## Gate 2 — geometry-only scale freeze

Only species identity and frozen coordinates enter this gate. A spherical `k = 5` nearest-neighbour graph is built within each admitted species. Three scales are evaluated in the declared order:

| Candidate resolution | Equal-area grid | Approximate cell area |
|---:|---:|---:|
| 100 km | 320 × 160 | 10,000 km² |
| 250 km | 128 × 64 | 62,500 km² |
| 500 km | 64 × 32 | 250,000 km² |

The first scale satisfying the frozen retained-edge, detectable-cell and cross-species opportunity gates is primary. All three remain required sensitivities. A scale can never be selected because it produces a stronger colour pattern.

Species failing geometry at a candidate scale remain `not_evaluable` at that scale. Cells lacking the required species opportunity are `not_evaluable`, not zero.

## Gate 3 — dated source, environmental coverage and flower ROI

Image pixels may be opened only after the cohort, exact 60,000-row v2 dated-source manifest, geometry and environmental opportunity coverage are hash-frozen. The measurement firewall requires the many-to-many dated-source pass and macroclimate plus at least one other primary environmental family at all 100/250/500 km scales. Acquisition also requires the exact firewall-key hash and a passing locked ROI result. The freeze manifest declares `sha256_lf_canonical_v1` for applicable text parents, so validation is identical on Linux and Windows. Acquisition retains the snapshot photo licence and attribution inputs.

Every photo must resolve to one of:

- `roi_ok` with a reproducible flower mask and continuous features;
- `not_evaluable` with a machine-readable failure reason.

There is no manual biological morph label and no post-hoc threshold rescue. ROI completeness and background-contamination gates are evaluated before coordinates and colour fields are joined for inference.

The first admissibility implementation, `fcp-inaturalist-automated-colour-state-v2`, is frozen separately in `docs/JBI_INATURALIST_AUTOMATED_COLOUR_STATE_PROTOCOL.md`. It uses a pinned CLIPSeg revision to produce soft **flower-candidate** weights, not verified flower-tissue masks. In its independent locked test, three development-passing species all passed completeness but none rejected the fixed spatial random-mark null. Consequently, that implementation is not automatically promoted to the 50-species cohort.

The proposed geographic transition statistic failed its prospective exact-geometry signal-recovery gate and is frozen `not_evaluable`; it is never applied to real atlas colour. ROI v3 also failed its independent JRC development gate and its locked test remains sealed. ROI v4 prospectively fixes a YOLO11n flower detector, EfficientSAM mask selection and one shared runtime used by both JRC qualification and scale-out. Its final training evidence must be committed before development prediction; development must pass before the locked 100-image test is opened. This is a measurement-validity gate, not an opportunity to retune against the observed three-species spatial results.

## Gate 4 — continuous colour and within-species spatial fields

Flower ROI pixels are converted to continuous perceptual-colour features. Calibration-only transformations may be species-specific; held-out colour values cannot tune representation or admission.

Primary spatial fields and local transition scores are reconstructed within species. Complete colour vectors are always permuted as indivisible rows and only within species. Species receive equal inferential weight after their own field is reconstructed.

## Gate 5 — species-conditioned transition boundaries

Within each species and frozen scale, edge discontinuities are rank-transformed to continuous transition intensity. Geometry fixes which edges and cells are detectable before colour is scored. The boundary object is therefore a species-specific continuous surface, not a universal categorical line.

## Gate 6 — cross-species environmental-boundary concordance

Within each cohort the surface averages ranked transition intensity equally across detectable species and retains the opportunity denominator; the eight cohort surfaces then receive equal weight. Flower and multiplicity-protected background surfaces are tested against independently frozen macroclimate, land-cover and ecoregion boundary intensity. The complete null moves each three-component colour vector as one row within species, rebuilds every eligible species surface and uses 9,999 weighted Moran eigenvector sign randomizations under one joint maximum statistic.

Support can establish spatial concordance between automated flower-candidate colour transitions and a predeclared environmental boundary family among the sampled species. It cannot establish a universal biological boundary, causation, adaptation or a pollinator mechanism.

## Gate 7 — species-free display

The atlas has two public-facing products:

1. a global point map in which each point is coloured by continuous flower ROI colour;
2. a photo bar of automated flower ROI thumbnails ordered without species labels.

Taxon, genus and family fields are removed from the display table. They remain present in the protected inferential table. The governing rule is:

> Species may disappear from the map and photo bar, but species must never disappear from the null model.

## Season and phenology sensitivity

The primary analysis uses all admitted dates with species-conditioned permutations. Two predeclared sensitivities preserve date structure more tightly:

1. permutation within species × observed calendar month where strata are sufficient;
2. permutation within species × hemisphere-adjusted local-solar quarter.

Sparse strata follow only the frozen adjacent cyclic merge rule or become `not_evaluable`. These checks diagnose seasonal photographic confounding; they do not identify a phenological mechanism.

## Stop rules and claim ceiling

- Metadata feasibility is not a colour result.
- Geometry support is not a shared-boundary result.
- A visible species-free band on the map is not inferential evidence.
- A failed scale, species or ROI gate is retained as failure or `not_evaluable`.
- No colour-dependent retuning, silent data repair or literature-led replacement is allowed.
- Environmental and pollinator overlays cannot rescue or replace the geographic test. Under the separately frozen v2 expansion contract they are ordered, multiplicity-adjusted concordance branches evaluated from the same locked colour fields; no new overlay may be selected after colour is opened.
- The completed three-species non-detection stops bulk atlas pixel opening until the independent estimator-qualification subgate passes; it cannot be bypassed by adding species or relaxing image gates.

## Prospective repeated-cohort expansion and ordered branches

The v2 expansion contract was frozen after the three-species negative validation and before any atlas candidate pixel was opened. It adds eight disjoint metadata-random panels of 25 species × 300 observations. The exact API selection must fully reconcile to one dated iNaturalist Open Data snapshot. The 200 species, 60,000 observations, seed rule, genus cap, no-overlap rule and no-early-stop rule are one experiment, not eight chances to select a favourable result.

The geographic branch is already `not_evaluable` from pre-colour method qualification. WorldClim 2.1, ESA WorldCover 2021 and RESOLVE Ecoregions 2017 are the evaluable frozen environmental family. Copernicus terrain is `not_evaluable` because its registered-source route did not complete. The *Bombus* branch is also `not_evaluable` because a citable authenticated GBIF download could not be frozen. Neither failed branch may be replaced after colour.

The expansion contract, independent-source audit, validator, simulations and CI are:

- `docs/supporting/jbi_image_first_atlas_expansion_contract_v2.json`;
- `docs/supporting/jbi_atlas_dated_source_amendment_v1.json`;
- `docs/supporting/jbi_atlas_dated_source_v1_stop_result.json`;
- `docs/supporting/jbi_atlas_dated_source_m2m_amendment_v2.json`;
- `docs/supporting/jbi_atlas_colour_surface_contract_v1.json`;
- `docs/research/FCP_ATLAS_GLOBAL_BOUNDARY_DATA_SOURCES.md`;
- `scripts/data/validate_jbi_image_first_atlas_expansion.py`;
- `scripts/data/benchmark_jbi_atlas_flower_roi.py`;
- `scripts/data/run_jbi_atlas_signal_recovery.py`;
- `.github/workflows/jbi-image-first-atlas-qualification.yml`.

## Implemented interfaces

- deep module: `fcp_pipeline/image_first_atlas.py`;
- live metadata adapter and freeze builder: `scripts/data/build_jbi_image_first_atlas_metadata.py`;
- validator: `scripts/data/validate_jbi_image_first_atlas_freeze.py`;
- contract/tests CI: `.github/workflows/jbi-image-first-atlas-metadata.yml`;
- focused tests: `tests/test_image_first_atlas.py`.
