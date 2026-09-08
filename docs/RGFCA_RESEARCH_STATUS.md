# RGFCA research status

Updated 2026-09-08 (Japan). RGFCA is the active FCP mainline. The six-species
Chapter 1 and 34-species literature study are legacy and are not input to the
current ecological-signal search.

## Current evidence

| Stage | Evidence | Interpretation |
|---|---|---|
| Measurement, two separate cohorts | 1,000 species; 100,000 terminal records; discovery/reserve classifiable counts 25,377 / 24,885 | Acquisition total, not a retroactively pooled inference frame |
| Eligible RGFCA frame | 369 species; 21,424 photos; at least 40 classifiable photos per species | Conditional inference frame, not a census of all plants |
| Discovery within-species spatial omnibus | Mean rho = 0.0270213; 999-permutation upper-tail p = 0.001; all 369 species included | Small positive exploratory photo-derived spatial association |
| Independent reserve primary | 363 species / 20,903 photos; mean rho = 0.0254826; p = 0.001 | Weak directional photo association replicated; not by itself flower-specific biology |
| Reserve observer / calendar controls | Mean rho = 0.0252180 / 0.0254826; both p = 0.001 | Both fixed controls passed; residual observational error is not excluded |
| Reserve matched-background differential | Mean rho = 0.0044773; p = 0.087 | Evaluable but not supported; flower-specific robust-replication gate did not pass |
| Discovery matched-background recovery | All 128 partitions complete; 21,339 exact rows and 85 reproduction failures among 21,424 | Not evaluable; no background-adjusted statistic or p-value; not an ecological negative |
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

Submission preparation also now includes a [JRC training-source audit](RGFCA_TRAINING_SOURCE_AUDIT.md)
and [saved ROI qualification reconstruction](RGFCA_ROI_QUALIFICATION_AUDIT.md).
The 100-image JRC test admitted 85 images; all eight frozen gate decisions
reproduce from the saved rows. This is European grassland box-validation
evidence, not global petal-mask accuracy: large-object recall is only 2/4 and
mask containment must not be read as petal segmentation IoU. Training-data
licence/citation checks do not complete scientific or package-wide release gates.

The frozen estimator also pools all retained flower instances without receiving
the observation's taxon label. Species-conditioning is therefore by labelled
photograph, not by a verified focal-species mask. Exact-source and artificial-mask
checks now document that scope; the contamination rate from co-photographed
flowers is unknown. Independent replication and background contrasts do not
themselves validate taxon-to-mask attribution. No measurement or admission rule
was changed to address this limitation after outcomes became available.

Run `34085861343` stopped at `ModuleNotFoundError: No module named 'fcp_pipeline'`
in all species shards, before numerical execution. Technical recovery exposes
the checkout through `PYTHONPATH` in all workflow jobs and smoke-tested both real
entry points. The data, runner, fixed seeds, statistics and thresholds are unchanged.

## Completed: prospective reserve-species replication, limited claim

The metadata-only audit identifies the entire complementary **500 taxa / 50,000
photos**, 100 per taxon, with no observation/photo-ID overlap with RGFCA discovery,
the older photo-first measurements, H9 fresh metadata or the H9 exclusion ledger.
All coordinates have stated accuracy at most 5 km; no observer contributes more
than two photos within a species. The original audit did not open images or
authorize measurement. Measurement and inference subsequently passed their
separate execution gates. See [the replication protocol](RGFCA_RESERVE_REPLICATION.md)
and [all completed results](RGFCA_RESERVE_REPLICATION_RESULTS.md).

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
  authorized at `1f80af7f5db81d61f28ae2818c130ef058a9f850`, now complete.
  Result commit `9f5abe7b45fcdc20ba83adf75d1a8d4a640f622c` retains all 50,000 rows
  from 500 species. All 256 CSV/receipt pairs were downloaded and hash-verified,
  then reassembled with the frozen completeness checks. The complete photo and
  observation census matches the metadata freeze and has no overlap with the
  recorded discovery/prior sets.
- 24,885 photographs are classifiable; the fixed at-least-40 gate admits **363
  species / 20,903 photographs**, passing the minimum-250-species rule. The other
  terminal records remain: 22,681 ROI/flip failures, 2,415 ambiguous palettes and
  19 records without biological-palette mass. These counts are measurement
  outcomes, not spatial inference or prevalence estimates.
- [Complete census audit](supporting/rgfca_reserve_complete_census_audit_v1.json)
  verifies all 512 terminal files, the completed output hashes and all 13 files
  from the successful inference implementation preflight. Exact Git bytes were
  used for frozen text; no hash rule was weakened for Windows newline conversion.
- [Fixed inference authorization](supporting/rgfca_reserve_inference_authorization_v1.json)
  at `391caecaa1a6d3c3e2407f5f19d0bad7057922e0` binds those 13 tested files and
  all three completed measurement files. No inferential outcome was opened to
  prepare that authorization.
- [Inference 34178957447](https://github.com/zuizui0223/fcp/actions/runs/34178957447)
  completed all 20 shards and all four tests once. Result commit
  `ef00a78a2eda7bc79f7e6b88f719a8a696dc8b43` records primary rho = 0.0254826,
  p = 0.001; observer and quarter p = 0.001; matched-background rho = 0.0044773,
  p = 0.087. All nulls are evaluable. The directional photo association replicated,
  but the required flower-specific robust-replication conjunction did not pass.
- [Full inference artifact audit](supporting/rgfca_reserve_inference_artifact_audit_v1.json)
  verifies all 363 x 4 x 999 stored null values, complete species/photo allocation,
  exact file hashes, 2,541 direct SciPy checks and the fixed bootstrap interval
  [0.0170055, 0.0343321]. That interval is conditional on the admitted species,
  not spatially or phylogenetically independent uncertainty for all plants.
  No new permutations were generated for the audit.

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

The authorized continuation [34094784607](https://github.com/zuizui0223/fcp/actions/runs/34094784607)
has now completed all 128 partitions and finalization. Result commit
`2f00847ccc4d15be8637a1e3589239ce2ab5c114` records
`not_evaluable_incomplete_exact_background_recovery`: 21,339 exact and 85 failed
rows, with no statistic, replacement or denominator adaptation. The failures
comprise 79 flower-pixel/palette mismatches and six separate background-pixel
mismatches. All original-image and admission flags pass; root cause and mismatch
magnitudes remain unknown. This completed technical limitation does not overturn
the discovery coefficient and does not validate it as flower-specific biology.
The [artifact completion audit](RGFCA_BACKGROUND_RECOVERY_COMPLETION.md) reconciles
every partition and unique measurement ID. Do not restart this closed recovery
or analyse its successful subset. The independent reserve has completed its fixed
inference; its evaluable but unsupported background differential is not a repair
of this recovery. A p-value of 0.087 is not proof of no effect or background causation.

## Discovery manuscript and figures

The active [RGFCA manuscript](RGFCA_MANUSCRIPT.md) now separates the completed
exploratory distance-colour signal and its weak reserve replication from
unsupported shared geography, unresolved mechanisms and the failed reserve
flower-specific gate. It is explicitly not submission-ready.
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
The subsequent reserve tests are now complete with the limited result above;
flower-specific validation, submission-wide references and final manuscript/SI
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
and added to the draft. The later JRC article/dataset audit brings the active
manuscript to 17 references. This does not close
the remaining release, model, data-provider or environmental-layer citation and
reuse audits. See [the software reading record](RGFCA_SCIENTIFIC_SOFTWARE_AUDIT.md).
The [measurement/provider audit](RGFCA_MEASUREMENT_PROVIDER_AUDIT.md) distinguishes
YOLO11 detection from the installed package's later YOLO26 citation, custom
weights from the base checkpoint, and EfficientSAM's pinned export from generic
segmentation claims. All five basemap files independently matched their official
tagged source. Mixed-licence metadata and model/training-data distribution still
need their own release review; the 24 CC0 crops do not clear the whole data pool.

## Technical maintenance: bounded legacy figure updates

The frozen legacy pipeline passed in runs
[34105526167](https://github.com/zuizui0223/fcp/actions/runs/34105526167) and
[34106658236](https://github.com/zuizui0223/fcp/actions/runs/34106658236), including
all five models, finite-sample checks, numerical regression and two-render figure
identity. Its automatic write step nevertheless produced consecutive bot commits
`a5c209a473f763b30d7d901d8f9ecab16e9e11b8` and
`fb5565cb587ef37dfb48283d7cb53c2dc8d97aec`, each changing only
`figure5_inference_method_sensitivity.pdf`; the second returns exactly to the
PDF bytes at `b3370a88d39ba90ce150bdccbc3817ccf818f45a`. Thus same-run rendering
identity did not prevent cross-run byte oscillation and repeated PR checks.

The final write step now stops for the original `github-actions[bot]` actor,
including after human approval/rerun. All analysis, numerical and rendering
checks still run; a human-originated update can still commit changed figures.
The distinction between original and rerun actors follows the
[GitHub variable contract](https://docs.github.com/en/actions/reference/workflows-and-actions/variables).
Six shell-replay regression cases intercept all Git operations: the unguarded
step failed four cases and passed two; the guarded step passed all six.
Together with the legacy unit tests, 18 tests passed using the exact stored LF
input. The Windows checkout and initial export had CRLF and failed the existing
exact-byte gate; exporting with `core.autocrlf=false` reproduced the required
`bdc06dd671f41ce062ebf4ba687437909d9617b268657504c1c6c5e991d417ed`
SHA-256 without changing the dataset or validator. This bounds write recursion;
it does not claim to remove the underlying cross-run floating-point PDF drift.
Reserve measurement, background recovery and all ecological outcomes are unchanged.

## Route to an ecological result

The [measurement-reference feasibility audit](RGFCA_TARGET_DOMAIN_VALIDATION_FEASIBILITY.md)
now covers six source families. A public FlowerMask metadata snapshot exposed
300 annotation-file entries, without proving exhaustive access to the advertised
3,600 photographs. A prospectively fixed six-document structure check completed
in [34182715550](https://github.com/zuizui0223/fcp/actions/runs/34182715550) at
`df0838c15f01eb98eaee6b620c3bb4fd5c9bc78e`, after the offline preflight passed.
All six exact file hashes and sizes matched. Each contained one polygon labelled
`Flower`, `flower` or `rose`, with no populated shape group ID; these fields do
not establish focal-taxon or petal-only truth. Embedded image strings were
transferred and hashed but not decoded, displayed or measured. The saved
[terminal ledger](supporting/rgfca_reference_schema_result_v1.json) is an access
and format result, not a measurement-performance result or ecological support.
The next useful step is reference image/annotation correspondence, provenance,
overlap and ontology qualification for a clearly limited flower-region diagnostic;
the independent iNaturalist focal-taxon measurement gate remains unresolved.
Do not rerun the six completed schema examples or admit all 300 from this sample.

**Monarda archive received and audited:** the completed user-supplied v1 ZIP
contains 110 paired images and 788 polygon annotations. All 110 source photo
and observation IDs are mapped through pinned author metadata, with zero ID
overlaps against either full 50,000-photo frame. This resolves the prior access
blocker. One image has no annotations and remains unknown; 98 photos have
historical noncommercial licence conditions. At that intake checkpoint, no pixels
were decoded, model run or measurement score calculated. See the
[complete intake and next prospective gate](RGFCA_MONARDA_REFERENCE_INTAKE.md).
ID disjointness is not model-training, event or observer independence, and
generic flower-region annotations are not verified focal-petal truth.

**Monarda execution clarification, before outcomes:** initial eight-test
qualification succeeded at `846b2ad` / run `34189809578`, but its runner called
the full CIELAB-computing runtime despite saying it never measured colour.
Both tasks confirmed no reference-image/model execution and assigned Monarda
to one desktop owner. The [explicit amendment](RGFCA_MONARDA_EXECUTION_AMENDMENT.md)
preserves the original contract, permits only incidental internal colour
computation, and requires exact source/environment checks before decoding,
a full 110-row durable ledger, and no automatic rerun. At this pre-outcome
checkpoint the hardened runner awaited qualification and separate authorization.
Neither CI nor this clarification alters
the reserve non-support or establishes botanical measurement accuracy.

**Subsequent completed Monarda result:** qualified code `1b6ccc7` and CI
`34193780128` preceded single-execution authorization `887b4da`. The run
completed on 8 September 2026, 06:32 UTC: **110/110 aligned, 110 terminal rows,
109 annotated + 1 unknown, zero runtime failures**. The fixed gate did not pass:
pooled precision **0.56824250** (< 0.70), pooled recall **0.42608787** (>= 0.35),
median image precision **0.51480059** (< 0.70). Fifteen annotated images had
empty predictions, retained under the original rules. An independent verifier
recomputed every count-derived metric, split summary and gate and checked
ledger/events/hashes. See [the retained result](RGFCA_MONARDA_REGION_AGREEMENT_RESULTS.md).
This is measurement-agreement non-support, not technical failure or ecological
absence. Localization error, annotation incompleteness and ontology mismatch
remain unresolved alternatives; no root cause is established. Monarda is now
outcome-opened and cannot be used as a fresh holdout for model selection.

The result head `6a9ab216766d99f4444d33041f5d4d085db593ce` passed actual
retained-result CI `34195682122`, reference `34195682118`, publication
`34195682181`, manuscript `34195686426`, boundary `34195686251` and legacy
`34195686254`. The legacy run completed both reproduction and double-render
PNG/PDF identity checks; generated figures and manuscript integration were
unchanged. CI completion does not turn the failed measurement gate into a pass.

**Next independent-reference gate:** the
[three-family source audit](RGFCA_INDEPENDENT_REGION_REFERENCE_AUDIT_20260908.md)
keeps all candidates unadmitted for focal-petal validation. USDA's public
100/30 filename lists share `IMG_0339.JPG` (union 129); exact original text
receipts and an offline verifier preserve the conflict, not a corrected split.
The [measurement-target design](RGFCA_MEASUREMENT_TARGET_AND_VALIDATION_V2.md)
separates generic flower-region engineering from focal-taxon and colour
validation. Next freeze a metadata-only paired-manifest check of the already
exposed 300 FlowerMask annotation IDs; no new image/model run occurred here.

1. Preserve the now-completed reserve measurement and all four fixed tests.
   The weak photo association replicated, but the flower-specific gate did not
   pass. Complete its bounded manuscript/reproducibility package, then address
   target-domain measurement and focal-species attribution before stronger claims.
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
