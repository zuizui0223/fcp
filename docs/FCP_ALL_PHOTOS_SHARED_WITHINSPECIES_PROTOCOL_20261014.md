# FCP: uncapped photographic opportunity → shared structure of *within-species* colour variation

**FRESH STUDY revision: 14 October 2026 — source photo collection and all ROI-v4 flower-colour measurements start anew.** The immutable old 149,900 measured-photo CSVs are **NOT biological inputs, not an exclusion list and not mixed with the new results**. They remain technical regression fixtures only. The input source is the frozen outcome-blind 42,111 taxon opportunity registry, not an old colour table. No biological fresh-photo results are claimed until the census/measurement completion receipts pass. The frozen original New Phytologist discovery/validation/prospective H1/H2 results and their source SHA identities remain unchanged.

## The scientific question, in one sentence

> **Across plant species, is within-species flower-colour variation repeatedly organized around the **same white–chromatic phenotype axis**, and does the *strength and organization of that intraspecific axis* vary with climatic and soil heterogeneity?**

A species' **average colour composition** regressed on world rainfall is a different estimand; it belongs in [`mois`](https://github.com/zuizui0223/mois), **not** here. Do not aggregate one photograph per species or turn within-species FCP into interspecific colour-frequency modelling.

## 1. New target frame: all photographable species and ALL available usable photographs

The first immutable sampling frame remains 42,111 **iNaturalist metadata-discovered named taxa**, not the world's plant inventory. The original post-observer-cap metadata census distinguished 8,753 species with ≥40 candidate images and 4,730 with ≥100, but those are *old opportunity counts*, not new inclusion ceilings. A species with fewer original candidates is **retained** in the full census, because a fresh source download may reveal more records. Only after image measurement do we decide which within-species statistics are supported.

- **No limit of 100 photographs/species. No global 473,000-photo ceiling. No selection using apparent colour, polymorphism, geography, environment, W or p-value.** Photo opportunity is metadata-only and auditable.
- Collect **all available distinct licensed photo IDs** from qualifying research-grade, georeferenced iNaturalist observations of the frozen taxa. Preserve **all** photos even when multiple photographs belong to one observation, but retain observation ID and observer ID. Distinguish the photographed trait (photo) from independent observation units and breeding plants.
- Freeze a `created_d2` upper bound and API/source version; traverse sorted observation IDs using `id_above` and record every cursor. Do **not** stop at page 50, at 10,000 observations, or once the species reaches a photo target. Document per-species `complete`, `in_progress`, inaccessible, or provider-failure status. For truly massive harvesting prefer source-provider bulk exports and metadata snapshots; follow [iNaturalist API recommended practices](https://www.inaturalist.org/pages/api+recommended+practices).
- No requery until a desirable amount of colour variation is detected, no replacement of hard-to-classify photographs. An unavailable URL or failed ROI is a **terminal unclassified photo**, not an observation to silently replace.
- Each rediscovered photo ID is downloaded and measured *afresh*, even if that photo ID occurred in the historical 149,900-row study. Do NOT concatenate historical colour/measurement rows with fresh rows. Count every distinct photo ID only once in this new study, and note that redownloading an earlier ID is not an independent biological replicate. Canonical taxon IDs and names require a logged taxonomy crosswalk if they drift.

Existing immutable sources: the full 42,111 frame at historical commit `7e538e5c51c05a7cc47b2fcf53eea92634c8a863`, `data/frozen/global_monte_carlo_species_discovery_combined_v2_species_v1.csv`; original three measured photo CSVs at `5142f7951af0dde5364bb047a566d67e8c479e51` and `7e538e5c51c05a7cc47b2fcf53eea92634c8a863`. U100=4,730 is retained as a historical subset, no longer the expansion population.

## 2. Same source-colour measurement as the 1,499-species study

Reuse the frozen **ROI-v4 detector → EfficientSAM flower mask → fixed nine-colour sRGB/CIELAB palette → four coarse colour states**, with every failed photo preserved as an unclassified outcome. The exact historical measurement runtime is at Git commit `9fae6ccdf684a46026f72ba12e98de2c5c54bf2a`; its original `scripts/analysis/measure_random_photo_first_blind_partition.py` and `fcp_pipeline/photo_first_measurement.py` control the image pipeline. The original detector SHA256 is `f1aaeec4664fe2c178e5cf2bc1f508977bef3e4aa7b40613026cb8ae3de789d5`; original EfficientSAM encoder/decoder hashes remain as frozen in the New Phytologist provenance contract.

The fresh source-audit code writes the **same three-field location-blind worker interface** `measurement_id,image_filename,photo_license` and keeps `photo_url_large`, species/taxon, observer, date and coordinates outside the worker. A complete sealed worker partition/terminal-result/metadata-join receipt is required before environmental information is available to colour measurement. Old H1/H2 numbers must not be recomputed using the later generic `disttrait` two-means implementation: near-tied k-means initializations have documented numerical differences.

This new acquisition code establishes the *metadata and measurement-firewall interfaces*, **not** an authorization to run all 42,111 taxa or a claim that downloaded image pixels were retained. An expensive full-source collection/measurement run should be separately budgeted, monitored, checkpointed and data-source-compliant.

## 3. Keep all photo information, avoid pseudoreplication

Measurement-grade distinction for every species:

- `not_measured` or `insufficient_classifiable` is missing support, never monomorphism.
- `n_classifiable ≥ 40` = unchanged paper-defined species-level H1/D and potential H2 **measurement** gate; a photographed species with one observed state is `single_photo_state_detected`, **not** genetically monomorphic.
- Only species with frozen second coarse colour-state fraction ≥0.10 and subsequently frozen continuous minor-cluster fraction ≥0.10 become eligible for per-species W; strict 0.20 is a sensitivity. A **candidate** passing the coarse gate has **not yet passed** the continuous H2 gate.
- Three descriptive D views: (1) photo-equal; (2) **observation-equal** (all photographs retained but each observation contributes one total weight); (3) **observer-equal**. The latter two are observation-process stress tests, not substitutes for the frozen photographic-D estimand. Preserve effective distinct observers, observations, locations and blooming seasons as audit fields.
- **Each species gets one vote in between-species inference**, regardless of how many photos are available. Species' sampling depth affects uncertainty and detection power; it must not act as a regression weight that makes popular species decide the ecological pattern.

## 4. Exactly three biological layers; keep their tests separate

| Structural feature within each species | Target | Evidence required across species |
|---|---|---|
| FCP amount | Four-state `D`, corrected pair diversity, uncertainty / observer-equal sensitivity | Observer-disjoint reproducibility; finite-sample, geography and missingness audits. Never infer population allele frequencies. |
| **Where FCP resides** | `D_specieswide − D_local(50 km)` and continuous palette turnover | Same original 50-km spatial vertex-permutation null (species-specific colour count and geography fixed) and controls for shared observer, location and time. Neighbourhood ≠ interbreeding population. |
| **Common form of FCP** | Unlabelled two-mode axis `W_i=(u_i·q_white)^2`, with `q_white=normalize([1,−1/8,...])` | The original **999-world coarse-state-preserving fixed-axis null** for aggregate W, 0.10 primary and 0.20 sensitivity. W sign is arbitrary; do not claim white arose from pigmented flowers. |

The **single main comparative endpoint is species-level fixed white-axis alignment `W_i`**, evaluated across species relative to the **frozen coarse-state-preserving null**, with **environmental breadth/heterogeneity inside each species** as primary ecological covariates. Species-wide D is a phenotype-opportunity gate and 50-km geographic allocation is a necessary diagnostic and ecological modifier, not competing headline outcomes. Do not maximize nominal significance by adding a new colour axis, selecting geographical radii post hoc, or cherry-picking species with apparent variation.

## 5. Environmental mechanism: directly annotate *each FCP photo*, not another study's nearby photo

Adopt the named and provenance-aware environmental axes from `mois` (WorldClim: BIO1–BIO19, mean/CV of solar, wind and vapour pressure, elevation; SoilGrids: pH, SOC, N, clay, available-water proxy, CEC, sand, silt, bulk density and coarse fragments). **Do not transfer 36 covariates from a nearby *mois* photograph and call them observations of the target flower.** The original strict photo-ID `mois` overlap was too small for a defensible full FCP soil result. New source-photo coordinates require fresh raster sampling, matching versions, units, pixel resolution, ocean/urban/no-data handling, georeferencing uncertainty and coverage receipts.

Separate two environmental estimands rather than treating all tests as the same:

1. **Across species:** does climatic/soil *heterogeneity within a species' sampled localities* explain variation in D, W or spatial partitioning across species? Model FCP metrics as the outcomes, species as replicates; account for taxonomic similarity, photograph opportunity, sampled geographic span and landscape/background availability.
2. **Within species:** are *differences between locality colour distributions* greater between environmentally different places even after accounting for geographic distance and observation history? Use species-specific matched-colour label permutations on fixed photos/sites, preserving the full colour-state composition; require observer-disjoint, exposure/flower-vs-background and matched-location controls. This is phenotypic, not proof of adaptation.

Signed white-minus-nonwhite mean environmental deltas are **secondary** because they test local colour-environment sorting, not a species-level FCP feature. The primary species comparison is `W_i ~ heterogeneity_i(environment)` with original-coarse-state conditioning and confounder/observation-opportunity checks; `D` and 50-km differentiation are only supporting diagnostics. Annual-average floral whiteness from `mois` is never the response. Neither the old BIO5 association nor any exploratory BIO7 or `mois` annual precipitation/soil coefficient becomes causal by rerunning more photographs. **Missing environments remain missing; never zero-fill.** Multiple environment-axis comparisons must be familywise controlled and accompanied by null/negative checks.

## 6. Exact NEW measurement sequence — ZERO reused historical photo colours

The historical 149,900 photos, H2 vectors and old environment tables are prohibited as **input** to the new biological estimation. Historic immutable scripts and trained weights are permitted **as fixed methods**, plus in separately named technical-regression jobs only. The pipeline MUST NOT use `--already-measured` and MUST NOT supply historical measured CSVs to either `fcp_allphotos` or `fcp_axis_replay` for the new study.

The following commands are a work plan (not evidence that the global fetch has run):

```bash
# 1. Original OUTCOME-BLIND 42,111 taxon identities; no species outcome selection.
# File SHA256 = 5d871bd1f1190c68fd513bb527e6c93de1606dc35ff8cfff5f6187850dc382cc
FRAME=/trusted/historical/outcome_blind_42111_species.csv
DATE=2026-10-14

# 2. Acquire one resumable shard of the full photo-source census.
# Compute page budget=10 is a checkpoint, NOT a terminal photo limit.
python scripts/acquisition/run_fcp_uncapped_shard.py \
  --frame "$FRAME" --output-root /trusted/fcp_fresh_metadata \
  --shards 128 --shard 0 --created-before "$DATE" \
  --pages-per-species-per-run 10
# Re-run the SAME shard/date until every species checkpoint is 'complete'.
# Repeat for all 128 shards; HTTP errors remain incomplete. Do not scale
# independent API clients beyond the provider's allowed overall rate.

# 3. Only after ALL 42,111 metadata taxa are complete:
python -m fcp_pipeline.fcp_fresh_study freeze \
  --source-species "$FRAME" --metadata-root /trusted/fcp_fresh_metadata \
  --shards 128 --worker-partitions 128 --outdir /trusted/fcp_fresh_firewall

# 4. Reacquire ALL source images. The original source/binary checksums are
# checked before opening pixels; no historical flower-colour labels are read.
# Run one job per (partition, work-batch). The 64-row value bounds memory
# only, NOT source photos or species. Repeat until every job is completed.
python scripts/acquisition/run_fcp_fresh_image_partition.py \
  --firewall /trusted/fcp_fresh_firewall \
  --partition 0 --batch 0 --batch-size 64 \
  --original-measure-script /trusted/original/measure_random_photo_first_blind_partition.py \
  --original-detector /trusted/original/jrc_yolo11n_last_v4.pt \
  --original-efficient-sam /trusted/original/efficient_sam \
  --output /trusted/fcp_fresh_sealed

# 5. A fresh biological measured table is authorized ONLY once every source
# photo ID has exactly one attested, SHA-checked terminal record.
python -m fcp_pipeline.fcp_fresh_study reassemble \
  --firewall /trusted/fcp_fresh_firewall \
  --sealed-results /trusted/fcp_fresh_sealed \
  --method-attestation /trusted/fcp_fresh_sealed/frozen_method_attestation.json \
  --outdir /trusted/fcp_fresh_complete

FRESH=/trusted/fcp_fresh_complete/newly_remeasured_photos.csv

# 6. NEW photo-only D/observer/observation support; no old measured inputs:
python -m fcp_pipeline.fcp_allphotos --opportunity "$FRAME" \
  --measured fresh_remeasured="$FRESH" \
  --outdir results/fcp_new_full_FCP

# 7. Original bit-identical frozen H2 geometry and 999-world coarse-state
# preserving null, with zero axis refitting on environmental predictors:
python -m fcp_pipeline.fcp_axis_replay --measured "$FRESH" \
  --opportunity "$FRAME" \
  --frozen-h2-source /trusted/original/run_polymorphism_delta_geometry_validation_20260912.py \
  --frozen-null-source /trusted/original/run_polymorphism_delta_geometry_structured_null_20260912.py \
  --cohort new_full_source --outdir results/fcp_new_W_axis

# 8. Annotate THIS photograph's own position with 36 provenance-pinned
# environmental raster fields, not nearby mois photo attributes:
python -m fcp_pipeline.fcp_direct_environment --measured "$FRESH" \
  --manifest /trusted/provenance_verified_mois36_rasters.json \
  --require-all36 --out /trusted/new_source_photo_environment.csv
python -m fcp_pipeline.fcp_environment_profiles \
  --measured-complete-species-partition "$FRESH" \
  --environment /trusted/new_source_photo_environment.csv \
  --outdir results/fcp_new_species_environment

# 9. Species-equal comparisons of W vs species-internal environmental spread.
# D and 50-km spatial FCP are descriptive/control metrics, not primary response.
python -m fcp_pipeline.fcp_shared_structure \
  --ledger results/fcp_new_full_FCP/all_photo_FCP_species_ledger.csv \
  --w results/fcp_new_W_axis/species_W_fixed_axis.csv \
  --environment-profile results/fcp_new_species_environment/species_environment_heterogeneity.csv \
  --environment results/fcp_new_species_environment/species_white_nonwhite_environment_sorting.csv \
  --outdir results/fcp_new_shared_FCP
```

**Attention:** Command 9 returns preliminary descriptive species-equal associations. A signed environmental selection result is **not** demonstrated by Spearman correlations. Species-level environmental opportunity nulls, observer/quality, photograph-background, geographic confounding, phylogenetic structure, raster nodata and multiple-testing adjustments are required for ecological inference. The original 999-world H2 null tests recurrent axis concentration, *not* a null distribution for the environmental slope; ecological associations need their own construction-preserving null. A taxon with insufficient photos still counts in the source opportunity denominator but is excluded from W estimation under unchanged H2 gates.

This study deliberately **does not claim** that all source metadata, original image pixels or 36 raster assets are now downloaded. Actual iNaturalist harvesting must respect provider guidance and may require an approved bulk-export route. A completed one-page pilot is never a completed world photograph census. The source taxon frame is a frozen outcome-blind historical identity layer; a genuinely new taxon discovery frame would require an additional provenance gate.

## 7. Current evidence ceilings (unchanged)

The **original 149,900 measured photos / 1,499 named species** produced **1,109 species with >=40 classifiable images**; **439 W-eligible** (152 +129 +158) across three independent species cohorts. A frozen 158-species prospective H2 white-axis test passed its construction-preserving null (`W=0.517246`, null median 0.457143, `p=0.001`), but it shares the same iNaturalist acquisition/measurement system. The three-cohort photographic 50-km depletion is replicated; signed white–BIO5 and BIO7 environment links do **not** pass strict observer/geographical transport. 4,730 is still a valid **old depth-opportunity census**, not the upper limit of this revised study. No new number of taxon-wise FCPs or direct soil/fitness explanation is asserted from unmeasured photos.

## 8. Explicit all-photo and soil limitations

- **Available photos is not a known finite constant.** A new iNaturalist crawl or bulk source export is required for a frozen census. `473,000` is neither the target nor an upper bound. All eligible original photo IDs may still include multiple images of the same observed plant; observer- and observation-equal versions are required.
- The historical sample is **149,900 measured photo rows**, not independently verified 149,900 breeding plants, and the historical 42,111-taxon metadata frame is an initial candidate frame rather than the biological total of possible taxa. New taxa can enter only through a separately audited outcome-blind frame amendment.
- Never call a taxon with <40 classifiable photographs monomorphic, or call species-wide white-versus-red geographic turnover local genetic FCP. Site-colour coexistence has to be analysed separately from between-site turnover.
- The source-matched `fcp_direct_environment` raster manifest enumerates exactly the `mois` 26 nonsoil and 10 soil **variable names**; **it does not invent or rederive a missing SoilGrids layer, raw unit conversion, or the available-water proxy**. Missing soil rasters remain unmeasured until the exact provenance-controlled 0–30cm/source-defined raster layers are available. Values are sampled at the *target photo's own location*, with location precision retained, and neither neighbouring photos nor missing-value interpolation count as direct observations.
- Structural white-axis `W` is sign invariant and its **between-species environmental regression cannot be called a confirmed pattern from ordinary Spearman p-values alone**. The H2 structured null, ecology-specific fixed-panel null, observer/flower-background controls, phylogenetic nonindependence and geographic opportunity must precede ecological inference. The old New Phytologist manuscript result is not overwritten by an expansion.