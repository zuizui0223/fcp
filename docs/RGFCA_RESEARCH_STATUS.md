# RGFCA research status

Updated 2026-09-07 (Japan). RGFCA is the active FCP mainline. The six-species
Chapter 1 and 34-species literature study are legacy and are not input to the
current ecological-signal search.

## Current evidence

| Stage | Evidence | Interpretation |
|---|---|---|
| Measurement | 500 species; 50,000 terminal records; 25,377 classifiable records | Automated public-photo colour measurements |
| Eligible RGFCA frame | 369 species; 21,424 photos; at least 40 classifiable photos per species | Conditional inference frame, not a census of all plants |
| Within-species spatial omnibus | Mean rho = 0.0270213; 999-permutation upper-tail p = 0.001; all 369 species included | Small positive exploratory photo-derived spatial association; independent replication pending |
| G1 repeated field | Primary p = 0.070; fine-scale sensitivity p = 0.006 | Primary not supported; scale-sensitive candidate structure |
| Species-disjoint commonness | p = 0.856; median fold correlation = -0.0880 | No supported transfer of boundary geography across held-out species |
| Five environmental process blocks | None passes the fixed five-block Holm gate; thermal partial rho = 0.00836, adjusted p = 0.050 | Weak exploratory thermal candidate, no confirmed mechanism |
| Environmental effect heterogeneity | No corrected main-effect variance, interaction variance or syndrome support | Species/context-specific environmental mechanisms are hypotheses |
| Real-climate synthetic qualification | Failed; moderate full-sharing recovery 0%, 3.2%, 6.8% across three blocks | Insufficient recovery under the tested design; no biological result |
| Sharedness-specific predictive qualification | Failed; moderate full-sharing recovery 0.8% | Method development, not evidence of absence in nature |

Result files are `docs/supporting/global_monte_carlo_measurement_result_v1.json`,
`global_rgfca_g1_result_v1.json`, `global_rgfca_prespecified_robustness_result_v1.json`,
`global_rgfca_species_disjoint_commonness_result_v1.json`,
`global_rgfca_expanded_environmental_panel_final_result_v1.json`,
`hypervolume_real_climate_synthetic_qualification_result_v1.json`, and
`global_rgfca_sharedness_specific_predictive_result_v1.json` in that same directory.
The heterogeneity inference result is on the parent branch, not at a local path:
[immutable result](https://github.com/zuizui0223/fcp/blob/f2f9c58e1d857a5d4b5b35a78adab8e8324b25eb/docs/supporting/global_rgfca_environmental_species_heterogeneity_inference_result_v1.json).
The [supporting evidence index](RGFCA_SUPPORTING_EVIDENCE.md) records its verified
artifact, exact file hashes and complete 5/10/2-family reconstruction checks.

## Completed: species-level spatial information

Completed the frozen 369-species within-species spatial omnibus under
`global_rgfca_within_species_spatial_omnibus_contract_v1.json`. This tests an
equal-species mean association between geographic distance and colour
dissimilarity. It is exploratory because earlier descriptive G3 outcomes are
known. Run [34088925008](https://github.com/zuizui0223/fcp/actions/runs/34088925008)
completed all 20 shards, all 999 global null statistics, eight focused tests and
12 direct SciPy equivalence checks. Result commit: `29584f3ad7ae0cd99a1d8f43459252af38f615da`.
The observed mean rho is **0.0270213**, null mean 0.0000396, upper-tail **p = 0.001**.
This is a small effect, not an independent replication, shared boundary,
environmental mechanism or estimate of the prevalence of spatial organization.
The null quantiles are not a confidence interval on the observed effect.
Full result: `supporting/global_rgfca_within_species_spatial_omnibus_result_v1.json`.

Run `34085861343` stopped at `ModuleNotFoundError: No module named 'fcp_pipeline'`
in all species shards, before numerical execution. Technical recovery exposes
the checkout through `PYTHONPATH` in all workflow jobs and smoke-tested both real
entry points. The data, runner, fixed seeds, statistics and thresholds are unchanged.

## In progress: prospective reserve-species replication

The metadata-only audit identifies the entire complementary **500 taxa / 50,000
photos**, 100 per taxon, with no observation/photo-ID overlap with RGFCA discovery,
the older photo-first measurements, H9 fresh metadata or the H9 exclusion ledger.
All coordinates have stated accuracy at most 5 km; no observer contributes more
than two photos within a species. The audit does not open images or authorize
measurement. See [the replication protocol](RGFCA_RESERVE_REPLICATION.md).

The same sampling platform, measurement model, regions and potentially observers
are shared: disjoint IDs do not eliminate systematic observational error.
Some species occupy few spatial cells; this is not uniform global coverage.
The newer upstream matched-background diagnostic is preserved separately. The
reserve measurement adds symmetric twelve-anchor background counts alongside
flower counts in the first decode; flower-specific interpretation requires the
prospectively specified differential control as well as the replication and
observer/season checks. Flower-only p = 0.001 alone is not a biological discovery.

Execution receipts:

- Reserve protocol and metadata freeze: `9fd4ae98632c74caa9e66dd24e7390ec221efb4c`.
- [Pre-pixel measurement preflight 34090865980](https://github.com/zuizui0223/fcp/actions/runs/34090865980):
  18 tests passed, full metadata audit and blind worker census reproduced.
- [Extended inference implementation preflight 34091922722](https://github.com/zuizui0223/fcp/actions/runs/34091922722):
  27 tests passed at `708b305e314087047a39d27956a76c9327585497`.
  This tested code without reading reserve outcomes while blind measurement
  was underway; it is not a second claim that no other process had opened pixels.
- [Reserve measurement 34091091640](https://github.com/zuizui0223/fcp/actions/runs/34091091640):
  authorized at `1f80af7f5db81d61f28ae2818c130ef058a9f850`, running.
  First partitions have completed; no partial colour inference is permitted.
- Fixed primary and all three controls are implemented. Inference execution
  remains unopened until all 256 terminal partitions and the exact cohort/ID
  census pass. The final inference authorization is not yet issued.

The separate discovery-background recovery encountered a numeric serialization
error (`1818.0` parsed as an integer string); upstream `af36e98` fixes exact-integer
parsing without changing masks, colours or statistics. An incomplete/failed
recovery is not a negative ecological result and cannot certify the flower-only
signal. Original run `34091174522` has now terminated with 128 failed recovery
partitions and no finalized background result. Preserve the original run and
track its properly recorded continuation using the already-fixed parser.
The complete failed-log audit found 128 tracebacks and 128 matching integer-string
errors; the sole artifact is the pre-acquisition firewall, not a recovered colour
result. [Technical continuation v2b](supporting/global_rgfca_background_control_numeric_recovery_v2b.json)
retains the same 21,424 photos, exact reproduction gates, finalizer and permutation
specification. Eleven focused parser regression tests pass and are required again
before the continuation opens images. Prior failed workers had decoded pixels;
this is explicitly not another pre-pixel freeze claim.

## Discovery manuscript and figures

The active [RGFCA manuscript](RGFCA_MANUSCRIPT.md) now separates the completed
exploratory distance-colour signal from unsupported shared geography, unresolved
mechanisms and prospective reserve validation. It is explicitly not submission-ready.
[Figure 1](figures/rgfca_figure1_discovery_atlas.png) displays all 21,424 eligible
discovery photographs without a species legend; [Figure 2](figures/rgfca_figure2_discovery_omnibus.png)
shows all 369 species effects and all 999 global null means. The map's colours
are display mixtures, not calibrated reflectance, and the null interval is not
an effect confidence interval. No reserve outcomes were read for these products.

The [figure contracts](RGFCA_PUBLICATION_FIGURES.md), reproducible builder and
source/output manifest retain exact discovery input hashes. Local verification
passed 11 tests, including complete-census guards, numerical reconstruction and
identical two-render PNG/PDF hashes. Both PNGs were visually inspected; an
observed-statistic label collision was corrected. The matching read-only CI
workflow reproduced these checks successfully in
[run 34094422559](https://github.com/zuizui0223/fcp/actions/runs/34094422559), with
artifact `10008033327` and archive SHA-256
`480ceaa25e0e5844dd482eff9e7e4bc4a2e695b5c9fb63a92f0f353f44a623f6`.
Independent validation, submission-wide references and final manuscript/SI
closure remain outstanding. The separate photo-bar update below now provides
licensed ROI crops without altering these first two figures.

## Completed: bounded core literature audit

The [image/ecology audit](RGFCA_IMAGE_ECOLOGY_LITERATURE_AUDIT.md) and
[statistical audit](RGFCA_STATISTICAL_LITERATURE_AUDIT.md) document twelve primary
papers, the sections actually read and explicit transfer/access limits. Eleven
are cited in the active manuscript; the additional survey-bias paper informs
the method-positioning note. This is a targeted reading, not a systematic review
or a claim that all submission references have been checked.

Direct floral-image precedents in *Erysimum*, *Monarda fistulosa* and North
American colour/phenology are acknowledged. Photo volume, image-first analysis
and pollinator overlays alone are not priority claims. The distinction between
one field's spatial structure and association between two spatial fields is
explicit; neither prior photographic validation nor the plus-one permutation
formula validates RGFCA's measurement or ecological interpretation.

Nine citation regression checks plus the eleven existing figure/evidence checks
passed locally and in [publication CI 34096712215](https://github.com/zuizui0223/fcp/actions/runs/34096712215)
(20 tests total). They reject uncited/duplicate references, stale
preprint metadata and missing reading records while preserving the incomplete
replication/submission labels. The read-only publication workflow now runs this
combined suite and preserves the audited draft/notes as a separate artifact.
Its four-document artifact `10008879104` was downloaded and every file matched
its committed Git object at `730612e4af6c64178b1c079d23bb6dc802e1dac0`.
The supporting-evidence extension adds nine checks: complete environmental Holm
families, all 17 heterogeneity tail probabilities, all 15 species variances and
preservation of failed synthetic qualifications. All 29 tests pass locally and
in [publication CI 34097657746](https://github.com/zuizui0223/fcp/actions/runs/34097657746).
Its five-document artifact `10009230975` was downloaded and each file matched
its committed Git object at `21d13586777eb90324b2ebfc54eeca8659e38b42`.
The same head's legacy 34-species CI `34097662439` also completed successfully.
No reserve outcome was read,
and no frozen measurement, inference, result or legacy file was changed by this
literature work.

## Completed: actual licensed discovery photo bar

The [display protocol](RGFCA_PHOTO_BAR.md) fixes 24 species/24 observers from the
659 CC0 photographs with credits in the discovery frame. Selection uses only
metadata after frozen eligibility, longitude-rank bins and a fixed hash rank;
neither colours nor observed effects choose the examples. The saved plan and
30 metadata/synthetic-pixel tests passed locally and in
[CI 34100227172](https://github.com/zuizui0223/fcp/actions/runs/34100227172).
All 24 current photo-level CC0 declarations, source-image hashes and original
ROI/palette summaries reproduced, followed by exact two-render PNG/PDF identity.
The downloaded artifact `10010284387` passed complete local crop/hash checks;
all 24 display slots were visually inspected. Sparse/fragmented masks are
retained without aesthetic replacement and do not establish segmentation accuracy.
The [release receipt](supporting/rgfca_photo_bar_release_v1.json) records provenance;
ongoing publication checks render committed crops without reacquisition.
This illustrative figure does not estimate global colour frequencies or count
as independent evidence; no reserve outcomes or legacy results are involved.

Four project-owned scientific-software citation records have also been checked
and added to the draft (15 manuscript references in total). This does not close
the remaining release, model, data-provider or environmental-layer citation and
reuse audits. See [the software reading record](RGFCA_SCIENTIFIC_SOFTWARE_AUDIT.md).
The [measurement/provider audit](RGFCA_MEASUREMENT_PROVIDER_AUDIT.md) distinguishes
YOLO11 detection from the installed package's later YOLO26 citation, custom
weights from the base checkpoint, and EfficientSAM's pinned export from generic
segmentation claims. All five basemap files independently matched their official
tagged source. Mixed-licence metadata and model/training-data distribution still
need their own release review; the 24 CC0 crops do not clear the whole data pool.

## Route to an ecological result

1. Freeze the reserve cohort, measurement, primary test and observer/quarter
   sensitivity before pixels, pass metadata and implementation checks, then
   measure the complete reserve once and retain every terminal outcome.
2. Develop bounded, registered exploratory questions using the existing atlas
   and qualification work. Local patchiness, monotone distance effects, shared
   boundaries and environmental mechanisms are distinct targets.
3. Register promising candidates for independent validation. Every existing
   outcome set is discovery data for hypotheses motivated by it; re-splitting
   previously explored species does not create an untouched validation set.
4. Produce species-free maps, species-conditioned inference, ecological
   interpretation, manuscript and reproducibility artifacts. Existing legacy
   positives cannot stand in for a new RGFCA result.

The active goal is in `RGFCA_RESEARCH_GOAL.md`. Technical and synthetic successes
alone do not complete it. All completed negative and non-evaluable outcomes stay
in the research record.
