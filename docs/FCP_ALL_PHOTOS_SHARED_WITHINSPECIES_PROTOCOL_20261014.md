# FCP: uncapped photographic opportunity → shared structure of *within-species* colour variation

**Revised scope: 14 October 2026, following explicit user decision to remove the 4,730 species × 100 photo ceiling.** This protocol is an *exploratory expansion plan and executable acquisition/support framework*, **not** a claim that unlimited photographs have already been measured. The frozen original New Phytologist discovery/validation/prospective H1/H2 results and their source SHA identities remain unchanged.

## The scientific question, in one sentence

> **Across plant species, is within-species flower-colour variation repeatedly organized around the **same white–chromatic phenotype axis**, and does the *strength and organization of that intraspecific axis* vary with climatic and soil heterogeneity?**

A species' **average colour composition** regressed on world rainfall is a different estimand; it belongs in [`mois`](https://github.com/zuizui0223/mois), **not** here. Do not aggregate one photograph per species or turn within-species FCP into interspecific colour-frequency modelling.

## 1. New target frame: all photographable species and ALL available usable photographs

The first immutable sampling frame remains 42,111 **iNaturalist metadata-discovered named taxa**, not the world's plant inventory. The original post-observer-cap metadata census distinguished 8,753 species with ≥40 candidate images and 4,730 with ≥100, but those are *old opportunity counts*, not new inclusion ceilings. A species with fewer original candidates is **retained** in the full census, because a fresh source download may reveal more records. Only after image measurement do we decide which within-species statistics are supported.

- **No limit of 100 photographs/species. No global 473,000-photo ceiling. No selection using apparent colour, polymorphism, geography, environment, W or p-value.** Photo opportunity is metadata-only and auditable.
- Collect **all available distinct licensed photo IDs** from qualifying research-grade, georeferenced iNaturalist observations of the frozen taxa. Preserve **all** photos even when multiple photographs belong to one observation, but retain observation ID and observer ID. Distinguish the photographed trait (photo) from independent observation units and breeding plants.
- Freeze a `created_d2` upper bound and API/source version; traverse sorted observation IDs using `id_above` and record every cursor. Do **not** stop at page 50, at 10,000 observations, or once the species reaches a photo target. Document per-species `complete`, `in_progress`, inaccessible, or provider-failure status. For truly massive harvesting prefer source-provider bulk exports and metadata snapshots; follow [iNaturalist API recommended practices](https://www.inaturalist.org/pages/api+recommended+practices).
- No requery until a desirable amount of colour variation is detected, no replacement of hard-to-classify photographs. An unavailable URL or failed ROI is a **terminal unclassified photo**, not an observation to silently replace.
- Confirm the same original photo IDs are not counted twice when the 149,900 original high-depth rows are combined with new metadata/measurements. Canonical taxon IDs and names require a logged taxonomy crosswalk if they drift.

Existing immutable sources: the full 42,111 frame at historical commit `7e538e5c51c05a7cc47b2fcf53eea92634c8a863`, `data/frozen/global_monte_carlo_species_discovery_combined_v2_species_v1.csv`; original three measured photo CSVs at `5142f7951af0dde5364bb047a566d67e8c479e51` and `7e538e5c51c05a7cc47b2fcf53eea92634c8a863`. U100=4,730 is retained as a historical subset, no longer the expansion population.

## 2. Same source-colour measurement as the 1,499-species study

Reuse the frozen **ROI-v4 detector → EfficientSAM flower mask → fixed nine-colour sRGB/CIELAB palette → four coarse colour states**, with every failed photo preserved as an unclassified outcome. The exact historical measurement runtime is at Git commit `9fae6ccdf684a46026f72ba12e98de2c5c54bf2a`; its original `scripts/analysis/measure_random_photo_first_blind_partition.py` and `fcp_pipeline/photo_first_measurement.py` control the image pipeline. The original detector SHA256 is `f1aaeec4664fe2c178e5cf2bc1f508977bef3e4aa7b40613026cb8ae3de789d5`; original EfficientSAM encoder/decoder hashes remain as frozen in the New Phytologist provenance contract.

The new source-audit code writes the **same three-field location-blind worker interface** `measurement_id,image_filename,photo_license` and keeps `photo_url_large`, species/taxon, observer, date and coordinates outside the worker. A complete sealed worker partition/terminal-result/metadata-join receipt is required before environmental information is available to colour measurement. Old H1/H2 numbers must not be recomputed using the later generic `disttrait` two-means implementation: near-tied k-means initializations have documented numerical differences.

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

## 6. Reproducible program and what is not yet run

```bash
# Outcome-blind inventory: every discovered species, no >=100 gate.
python -m fcp_pipeline.fcp_allphotos \
  --opportunity /trusted/42111_species_frame.csv \
  --outdir results/fcp_allphoto_inventory

# Complete one species metadata enumeration (no max batches and no per-species cap):
python scripts/acquisition/collect_fcp_uncapped_metadata.py \
  --taxon-id 12345 --species 'Example species' --created-before 2026-10-14 \
  --out-dir /frozen/metadata/12345

# Freeze completed metadata without opening pixels; source URL stays sealed:
python scripts/acquisition/build_fcp_uncapped_blind_batches.py \
  --metadata-dir /frozen/metadata/12345 \
  --output /frozen/blind_workers/12345

# Once frozen ROI-v4 acquisition/measurement and terminal reassembly exist,
# ALL photos are accepted, not a 100-photo sample:
python -m fcp_pipeline.fcp_allphotos \
  --opportunity /trusted/42111_species_frame.csv \
  --measured original_discovery=/trusted/discovery.csv \
  --measured original_validation=/trusted/validation.csv \
  --measured original_prospective=/trusted/prospective.csv \
  --measured all_photos_new=/frozen/all_terminal_measured.csv \
  --outdir results/fcp_allphoto_FCP

# Original H2 *mathematical implementation*, exact historical blob,
# on every classifiable measured photo; does not refit any novel axis.
python -m fcp_pipeline.fcp_axis_replay \
  --measured /complete/species_partition/photo_colour.csv \
  --opportunity /trusted/42111_species_frame.csv \
  --frozen-h2-source /historical/7e538e/scripts/analysis/run_polymorphism_delta_geometry_validation_20260912.py \
  --cohort all_photo_expansion --outdir /frozen/W_expansion
# Inferential W claim requires --frozen-null-source from the same original
# Git commit to execute 999 coarse-state-preserving worlds.

# Direct photo-coordinate raster sampling, same 36 variable names as mois.
# Raster versions, units/depth, monthly reduction and SHA256 must be supplied.
python -m fcp_pipeline.fcp_direct_environment \
  --measured /complete/species_partition/photo_colour.csv \
  --manifest /trusted/fcp_mois36_rasters.json --require-all36 \
  --out /frozen/photo_environment.csv
python -m fcp_pipeline.fcp_environment_profiles \
  --measured-complete-species-partition /complete/species_partition/photo_colour.csv \
  --environment /frozen/photo_environment.csv --outdir /frozen/species_environment

# Evidence-only cross-species analysis, W primary and D/spatial diagnostic.
python -m fcp_pipeline.fcp_shared_structure \
  --ledger results/fcp_allphoto_FCP/all_photo_FCP_species_ledger.csv \
  --w /frozen/study_specific_H2_species_W.csv \
  --spatial /frozen/50km_species_metrics.csv \
  --environment /frozen/species_environment/species_white_nonwhite_environment_sorting.csv \
  --environment-profile /frozen/species_environment/species_environment_heterogeneity.csv \
  --outdir results/fcp_uncapped_shared_FCP
```

**Original measured photo IDs remain frozen and authoritative.** The blind batching script accepts `--already-measured discovery=...` (similarly validation and prospective), checks the original three data SHA256 values, and excludes those photo IDs from the new acquisition queue without dropping their historical measurements. Exact duplicate retry rows are reconciled; conflicting photo metadata halt the run. A completed metadata manifest is not a fully measured photo dataset.

The archive-based CLI must take *actual* photo/image sources; a synthetic or incomplete metadata result cannot be labelled `full`. Metadata and source URL acquisition, pixel processing, 999-world H2 null, 50-km label-permutation null and environmental raster joins are distinct gates. Input photos already used in the 1,499-species historic sets cannot be counted twice. The new module's own output **is an inventory, not a fully performed 42,111-species biological analysis**.

## 7. Current evidence ceilings (unchanged)

The **original 149,900 measured photos / 1,499 named species** produced **1,109 species with >=40 classifiable images**; **439 W-eligible** (152 +129 +158) across three independent species cohorts. A frozen 158-species prospective H2 white-axis test passed its construction-preserving null (`W=0.517246`, null median 0.457143, `p=0.001`), but it shares the same iNaturalist acquisition/measurement system. The three-cohort photographic 50-km depletion is replicated; signed white–BIO5 and BIO7 environment links do **not** pass strict observer/geographical transport. 4,730 is still a valid **old depth-opportunity census**, not the upper limit of this revised study. No new number of taxon-wise FCPs or direct soil/fitness explanation is asserted from unmeasured photos.

## 8. Explicit all-photo and soil limitations

- **Available photos is not a known finite constant.** A new iNaturalist crawl or bulk source export is required for a frozen census. `473,000` is neither the target nor an upper bound. All eligible original photo IDs may still include multiple images of the same observed plant; observer- and observation-equal versions are required.
- The historical sample is **149,900 measured photo rows**, not independently verified 149,900 breeding plants, and the historical 42,111-taxon metadata frame is an initial candidate frame rather than the biological total of possible taxa. New taxa can enter only through a separately audited outcome-blind frame amendment.
- Never call a taxon with <40 classifiable photographs monomorphic, or call species-wide white-versus-red geographic turnover local genetic FCP. Site-colour coexistence has to be analysed separately from between-site turnover.
- The source-matched `fcp_direct_environment` raster manifest enumerates exactly the `mois` 26 nonsoil and 10 soil **variable names**; **it does not invent or rederive a missing SoilGrids layer, raw unit conversion, or the available-water proxy**. Missing soil rasters remain unmeasured until the exact provenance-controlled 0–30cm/source-defined raster layers are available. Values are sampled at the *target photo's own location*, with location precision retained, and neither neighbouring photos nor missing-value interpolation count as direct observations.
- Structural white-axis `W` is sign invariant and its **between-species environmental regression cannot be called a confirmed pattern from ordinary Spearman p-values alone**. The H2 structured null, ecology-specific fixed-panel null, observer/flower-background controls, phylogenetic nonindependence and geographic opportunity must precede ecological inference. The old New Phytologist manuscript result is not overwritten by an expansion.
