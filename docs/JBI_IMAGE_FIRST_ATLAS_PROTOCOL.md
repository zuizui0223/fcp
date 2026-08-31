# FCP image-first global flower-colour atlas — prospective protocol

## Mainline and frozen boundaries

The active FCP mainline is now an **image-first global flower-colour atlas**. The ordered analysis is:

`iNaturalist image -> automated flower ROI -> continuous colour -> within-species spatial field -> species-conditioned transition boundary -> cross-species shared-boundary concentration`

Species identity is hidden only in the atlas display. It remains mandatory in every inferential graph, standardization step and permutation null.

Two completed analyses are retained without reopening them:

1. the six-species, 1,200-photograph Chapter 1 development result (Stage A `p = 0.0113`; Stage B `p = 0.0906`);
2. the literature-derived 34-species comparative analysis, retained as a classification method and supplementary comparison rather than the atlas admission frame.

## Gate 0 — pre-image contract

The v1 contract is frozen at `docs/supporting/jbi_image_first_atlas_contract_v1.json`. Before any candidate image pixel is opened, it fixes:

- the iNaturalist metadata root and flowering annotation;
- the 50-species pilot size;
- 300/400/500-photo sample tiers;
- licensing and positional-accuracy rules;
- observer, spatial-cell and month caps;
- 0.25° primary and 0.5° sensitivity thinning;
- geometry-only candidate scales of 100, 250 and 500 km;
- season/phenology sensitivities;
- display and inferential separation;
- fail-closed stopping rules.

Literature class, flower colour, Stage A/B effects, environmental layers and candidate image pixels are prohibited admission inputs.

## Gate 1 — metadata feasibility and cohort admission

The pilot begins from research-grade, flowering-annotated iNaturalist records under the frozen Angiospermae root (`taxon_id = 47125`). Records require a licensed photograph, public coordinates, a numeric positional accuracy no greater than 5 km, a date and a non-captive status.

Species are ranked by the count returned by that metadata query. The deterministic selector traverses that ranking and admits the first 50 species passing all rules. At most one species per iNaturalist genus is admitted. The completed six development species are excluded so their frozen evaluation is not recycled as a new confirmatory cohort.

Within species, metadata rows are balanced without opening photographs:

- one licensed photograph per observation;
- 0.25° spatial round-robin selection;
- simultaneous cap of 5 records per 0.25° cell;
- cap of 10 records per 0.5° cell;
- cap of 10 records per observer;
- cap of 100 records per calendar month;
- minimum 50 observers, 60 primary cells, 40 sensitivity cells, four months and three hemisphere-adjusted local-solar quarters.

The largest passing tier is retained in the fixed order 500, 400, then 300 photographs. A failure is recorded and the selector continues through the already ranked metadata pool. If 50 species cannot be admitted, the workflow stops. It does not weaken gates or use colour to choose replacements.

The pilot uses the read-only API for bounded feasibility. Subsequent bulk image acquisition should resolve the frozen photo IDs through the monthly iNaturalist Open Data snapshot rather than scrape the public API.

## Gate 2 — geometry-only scale freeze

Only species identity and frozen coordinates enter this gate. A spherical `k = 5` nearest-neighbour graph is built within each admitted species. Three scales are evaluated in the declared order:

| Candidate resolution | Equal-area grid | Approximate cell area |
|---:|---:|---:|
| 100 km | 320 × 160 | 10,000 km² |
| 250 km | 128 × 64 | 62,500 km² |
| 500 km | 64 × 32 | 250,000 km² |

The first scale satisfying the frozen retained-edge, detectable-cell and cross-species opportunity gates is primary. All three remain required sensitivities. A scale can never be selected because it produces a stronger colour pattern.

Species failing geometry at a candidate scale remain `not_evaluable` at that scale. Cells lacking the required species opportunity are `not_evaluable`, not zero.

## Gate 3 — image and flower ROI

Image pixels may be opened only after the cohort, observation manifest and geometry scale are hash-frozen. The freeze manifest declares `sha256_lf_canonical_v1`, which canonicalizes text newlines to LF before hashing so validation is identical on Linux and Windows. Acquisition retains the photo licence and attribution. Flower localization is automated and uses a separately frozen model revision and thresholds.

Every photo must resolve to one of:

- `roi_ok` with a reproducible flower mask and continuous features;
- `not_evaluable` with a machine-readable failure reason.

There is no manual biological morph label and no post-hoc threshold rescue. ROI completeness and background-contamination gates are evaluated before coordinates and colour fields are joined for inference.

The first admissibility implementation, `fcp-inaturalist-automated-colour-state-v2`, is frozen separately in `docs/JBI_INATURALIST_AUTOMATED_COLOUR_STATE_PROTOCOL.md`. It uses a pinned CLIPSeg revision to produce soft **flower-candidate** weights, not verified flower-tissue masks. In its independent locked test, three development-passing species all passed completeness but none rejected the fixed spatial random-mark null. Consequently, that implementation is not automatically promoted to the 50-species cohort.

Before bulk atlas pixels can be opened, a new estimator-qualification subgate must freeze and pass (i) an independent flower-tissue localization benchmark and (ii) signal-recovery simulations for the exact downstream transition statistic. This is a measurement-validity gate, not an opportunity to retune against the observed three-species spatial results.

## Gate 4 — continuous colour and within-species spatial fields

Flower ROI pixels are converted to continuous perceptual-colour features. Calibration-only transformations may be species-specific; held-out colour values cannot tune representation or admission.

Primary spatial fields and local transition scores are reconstructed within species. Complete colour vectors are always permuted as indivisible rows and only within species. Species receive equal inferential weight after their own field is reconstructed.

## Gate 5 — species-conditioned transition boundaries

Within each species and frozen scale, edge discontinuities are rank-transformed to continuous transition intensity. Geometry fixes which edges and cells are detectable before colour is scored. The boundary object is therefore a species-specific continuous surface, not a universal categorical line.

## Gate 6 — cross-species shared-boundary concentration

The shared surface averages transition intensity only across species detectable in a cell and retains the opportunity denominator. The complete null pipeline permutes colour within species, then recomputes edges scores, species-cell intensities, the shared surface and concentration.

Support can establish repeated geographic concentration among the sampled species under the frozen scale. It cannot by itself establish a universal biological boundary, a shared mechanism, climatic cause or adaptation.

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

The v2 expansion contract was frozen after the three-species negative validation and before any atlas candidate pixel was opened. It adds eight disjoint metadata-random panels of 25 species × 300 observations from one dated iNaturalist Open Data snapshot. The 200 species, 60,000 observations, seed rule, genus cap, no-overlap rule and no-early-stop rule are one experiment, not eight chances to select a favourable result.

All evaluable geographic, environmental and pollinator branch statistics enter one complete nested maximum-statistic null. WorldClim 2.1, Copernicus DEM GLO-30, ESA WorldCover 2021 and RESOLVE Ecoregions 2017 are frozen before colour joins. A *Bombus* layer must be derived without flower colour from one citable GBIF occurrence download and pass its own realm, opportunity-cell, source-balance and stability gates; otherwise that branch is `not_evaluable`.

The expansion contract, independent-source audit, validator, simulations and CI are:

- `docs/supporting/jbi_image_first_atlas_expansion_contract_v2.json`;
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
