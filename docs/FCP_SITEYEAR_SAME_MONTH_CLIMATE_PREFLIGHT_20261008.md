# Preflight: same-species, same-month, repeated-site-year photographic opportunity — 2026-10-08

## Scientific question and constrained immediate objective

The three fixed historical FCP cohorts show local photographic colour clustering even after month and year×month label nulls. **Those conditionings do not test real interannual climatic variability.** Before retrieving ERA5-Land or similar year-specific climate data, measure whether the historical photographs permit *comparisons of the same species at a physically bounded site in the same calendar month across at least two different years*.

This is a new **retrospective, post-outcome coverage audit**. It is not a flower-colour result, not an environmental regression and not a revision to the frozen New Phytologist, H1/H2 or pending prospective 2,000+730-taxon programme. It cannot remove the already-known imaging/exposure coupling of photographic white.

## Immutable photo sources

Same files, identities and SHA256 as the completed species-wide spatial-season analysis:
- Discovery: `global_monte_carlo_measured_photos_v1.csv`; SHA256 `ee854126eed2cfe23e333abe2c28d14df24895a5c52cbc389c060a5a4d6f91f4`.
- Validation: `rgfca_reserve_replication_measured_photos_v1.csv`; SHA256 `0e2ed349122739eecfc725fb2d5e313d284cf30752da91f9a0429cff0eeaa5e6`.
- Third cohort: `polymorphism_h2_third_cohort_measured_photos_v1.csv`; SHA256 `57630fc9f281bce94a0c40a70aaf7bce879dde93d6154175adcd021e8f5c1186`.

Keep species with >=40 already classifiable flower photographs (existing 369/363/377 denominator), valid coordinates and observed calendar dates 1990–2026. Only dates, observer identifiers, locations and photo identities are exposed to site/year selection: **the fixed four-category colour label is dropped after the historical classifiability gate**. Classifiability itself may depend on phenotype/lighting, which this audit cannot undo. No new pixels or climate data are retrieved.

## Frozen site-year-month opportunity rule

1. **Primary physical locality**: photographic point set with a maximal pairwise great-circle *diameter* <=10 km. For reproducibility, every candidate is anchored on a source photograph and includes only photographs within 5 km of that anchor (a conservative radius/2 construction); the actual maximum pairwise distance is checked. This construction can miss some legitimate <=10-km sets and must not be interpreted as exhaustive co-locality identification.
2. **Same flowering season**: all candidate photos come from the *same exact calendar month* in every qualified year; pooling March and June does not create repeated-month support.
3. **Repeat years**: >=2 different calendar years. Each qualifying year×month must contain >=2 classifiable source photographs by >=2 distinct identified observers, all within the same anchored locality. This avoids pseudo-repeat from a single user's seasonal photos, but different photographs/observers need not imply distinct plant genotypes.
4. **Outcome-independent opportunity ranking**: among eligible anchor×month candidates, prefer more qualifying years, then more source photos, then more observers; deterministic photo ID/month tie-breaks. Report one best qualifying locality per species and each radius without extrapolating it to total sites per species.
5. **Sensitivities**: 25 and 50 km *diameter* using the same anchored half-radius algorithm. Their stronger coverage does not rescue a failed primary 10-km analysis; a 50-km photographic locality cannot be equated to a mating population.
6. **Feasibility threshold**: >=30 eligible source species in *each* disjoint cohort at 10 km before planning a comparative year-anomaly model on these historical sources. This is a **resource/coverage gate**, not a statistical significance rule. Source photographs were geographically maximin sampled and may not resolve same-site years even when the underlying natural populations are biologically stable.

Output is a cohort/scale coverage summary and an outcome-blind per-species metadata CSV. Unqualified species are marked *photo opportunity insufficient*, never monomorphic, temporally invariant or climate-insensitive.

## Interpretation and next step

If 10-km coverage clears 30 eligible species per cohort, **do not immediately claim biological identification**. Before climate model fitting:
- Freeze the year-specific source, monthly value definition (e.g. flowering-month temperature/drought anomaly versus historical same-month site baseline), temporal coverage/lat-lon resolution, source hashes, request cost and extraction rules **without using photo-colour effects**.
- Compare **within the same source species and bounded locality across years**, with fixed season, source/observer checks, matched alternative photo exposure/background metrics and no extrapolation of unobserved sites.
- Check temperature, precipitation, water-deficit and radiation anomalies as competing *predeclared* environmental mechanisms, while recognizing that observer/camera, phenology, garden/wild and genotype turnover can still confound.
- Require ecological replication across disjoint taxa, and independent pigment/organism-level measurements before interpreting colour change as adaptive rather than plastic or photographic.

If 10-km coverage is below 30 in any cohort, record a frozen HOLD on *historical same-site climate inference*. Do not manufacture pseudo-time series by attaching static WorldClim means to observations with different years. Do not opportunistically increase site radius until a nominally significant result appears.

## Run

```bash
python -m pytest tests/test_fcp_siteyear_same_month_opportunity_20261008.py -q
python scripts/analysis/audit_fcp_siteyear_same_month_opportunity_20261008.py \
  --discovery /tmp/discovery.csv \
  --validation /tmp/validation.csv \
  --third /tmp/third.csv \
  --outdir /tmp/fcp_siteyear_opportunity
```

The workflow re-fetches exact historical Git commits and verifies each input's SHA256. Results become verified observational coverage *only after* tests, input hashes and receipt checks succeed.


## Terminal original anchor-based result (source-verified, 2026-10-08)

GitHub Actions run [37787901739](https://github.com/zuizui0223/fcp/actions/runs/37787901739), job `113347305008`, passed 8 synthetic tests, all 3 original photo SHA256 checks and full result validation; artifact `11556301126`, SHA256 `c959bc3ced3921a8a9d63289cb190523b8bf5aa777b31c1827f13cb868d35e8a`. Frozen human- and machine-readable results are preserved at `results/fcp_siteyear_same_month_opportunity_20261008/terminal_coverage_receipt.json`.

| Cohort | Original >=40 photo species | 10-km photo-centred half-radius anchor | 25-km diameter | 50-km diameter |
|---|---:|---:|---:|---:|
| Discovery | 369 | 10 | 43 | 66 |
| Validation | 363 | 10 | 22 | 50 |
| Third | 377 | 20 | 40 | 81 |

All original 10-km precommitted **30 eligible species per cohort** gates **FAILED**. `HOLD_10KM_REPEATED_SITE_YEAR_MONTH_PHOTO_OPPORTUNITY` is the terminal decision for the originally frozen conservative anchored definition. The 25-km and 50-km findings cannot rescue it. This is an observation-coverage HOLD, not proof of zero interannual climate effects or biologically stable flower colour.

### Separate exploratory geometry check, not an amended primary threshold

Photo-anchored circles with radius 5 km are a conservative implementation of a 10-km site *diameter*. Four photographed locations can have maximum pairwise distance <10 km yet no photograph within 5 km of every other photograph. This could lead to a method-driven false-negative coverage classification for an otherwise valid same-month, different-year, two-observer candidate site.

After seeing the first HOLD, a separate method-sensitivity script, `scripts/analysis/audit_fcp_exact_10km_siteyear_opportunity_20261008.py`, was introduced. It exhaustively tests **existence** of a four-photo minimal witness: two observer-distinct photograph pairs, one per year, with all six between-photo distances <=10 km and exactly the same month. The existence test is mathematically sufficient and necessary for the *minimum four-photo* year-month requirement, but does not identify a sampled genetic population or optimize for >2 years. It neither changes the original sampling rule nor serves as a prospectively registered confirmation. It must be evaluated from checked Actions outputs, not asserted before execution.

Regardless of geometric sensitivity, a retrospective year-anomaly flower-colour association still lacks field-verified genetic morph identity, observation-level camera/illumination matching, and a genotype-aware local change design.


### Terminal exact-diameter diagnostic — verified, does not rescue coverage

After the original anchor-HOLD, the four-photo minimal-witness exact 10-km geometric sensitivity passed its seven additional synthetic tests and source-identity replay in [Actions 37792424151](https://github.com/zuizui0223/fcp/actions/runs/37792424151); checked artifact `11557705349` (SHA256 `245c3cea172fd1fb97e3ca863752470416a62079464fe2cebaf34d2c9968bbf5`). Permanent [numerical receipt](../results/fcp_exact_10km_siteyear_opportunity_20261008/terminal_geometry_receipt.json).

| Cohort | Original conservative anchor at 10km | Exhaustive exact FOUR-photo witness at 10km | Geometric undercount corrected | Meets 30-species exploratory threshold? |
|---|---:|---:|---:|---|
| Discovery | 10 | 19 | +9 | No |
| Validation | 10 | 12 | +2 | No |
| Third | 20 | 23 | +3 | No |

This exposes a real **geometry-induced undercount of 14 species** in the original conservative anchor screen; the corrected minimum-witness existence calculation still falls short of 30 eligible species in **every** cohort. The original anchor-based precommitted HOLD is unmodified. Even under the retrospectively strengthened exact test, the species-pool and within-site photographs cannot sustain the proposed cross-cohort comparative same-site yearly anomaly analysis.

**Terminal ecological verdict:** the historical maximin-style, high-depth iNaturalist photograph resource is strong for comparing flower-colour organization across a sampled range, but too sparse in independently observed, strictly repeated local year×month windows for this distinct causal-environmental question. The extent to which the sampling algorithm versus underlying photographer opportunity caused low coverage was not itself decomposed. No absence of climatic effects, genetic polymorphism or selection is established.

**Next data-design boundary:** any climate anomaly test must begin from a separately collected or predesignated repeated-*site* sample, not enlarge the locality diameter after seeing 10-km failures, or repurpose old same-provider classes as untouched prospective replication. Required are site/year/month/observer windows with enough photos, camera exposure controls, repeated source plants or genotype identities, and matched year-specific climate covariates. With the current frozen photo sample, **do not download ERA5-Land for the blocked three-cohort same-site test**.


## Follow-on data-recruitment gap audit (outcome blind; NOT a re-test of spatial or climate effects)

After both 10-km tests remained HOLD (conservative anchors 10/10/20, exact four-photo witness 19/12/23), the study next asks a different *feasibility* question: **for the existing frozen high-depth species, where are there already two photographed years of the same month in a 10km neighbourhood, but too few observer-disjoint photos per year?**

The additional `scripts/analysis/audit_fcp_siteyear_replenishment_gap_20261008.py` searches only source-photo-centered 5-km anchor circles and only **years that already have actual dated photographs**. For two years in the same exact month, it computes a theoretical observer-photo slot shortfall
`max(0, 2 - n_observers_in_year_A) + max(0, 2 - n_observers_in_year_B)`.
The minimum over anchors, calendar months and source-observed year pairs is reported, while species without a two-year/month photographic anchor are separately coded *not observable in the current source*. All photo-colour labels have been dropped before selecting locations and dates; individual photo and observation IDs remain historically source-verified.

A gap of one is a **mathematical recruitment target**, NOT evidence that the missing observation exists, that the corresponding species colour is genetically variable, or that climate explains its variation. The chosen site and date may be biased by source participation and photography. Because the source photograph-selection design maximized geographic coverage, this is specifically a *supplemental site-first observational* requirement for existing taxa, not an independent confirmation cohort.

The original ≥30-per-cohort 10-km source-only gate remains HOLD, regardless of results. No new iNaturalist requests or images are fetched, and no future 2,000+730 taxon allocation is opened by this audit. Results are to be interpreted only after the additional pinned CI succeeds and commits a count receipt.


## Verified replenishment shortfall results (2026-10-08)

The source-SHA checked [GitHub Actions run 37793952964](https://github.com/zuizui0223/fcp/actions/runs/37793952964) completed the original eight, exact-geometry seven, and missing-observer-slot eight synthetic tests (**23 passed**), and emitted the label-blind site/year/month gap ledger as artifact `11558485103` (SHA256 `d2a2627d92243350bdc11cce8afc034d5bd77aeac7de2b332d89767a087898ee`). The durable outcome-limited [JSON receipt](../results/fcp_same_siteyear_replenishment_gap_20261008/terminal_gap_receipt.json) is extracted verbatim from the successful source-verified workflow output.

| Species cohort | >=40 classified-photo source taxa | Already anchored 10km complete (gap 0) | Same-month/year-local gap 1 observer-photo | Gap 2 | No repeated year/month 10km source-photograph anchor |
|---|---:|---:|---:|---:|---:|
| Discovery | 369 | 10 | 37 | 49 | 273 |
| Validation | 363 | 10 | 22 | 60 | 271 |
| Third | 377 | 20 | 33 | 62 | 262 |
| **Total** | **1,109** | **40** | **92** | **171** | **806** |

The **303** species with an already observed repeat-year/month 10km photo anchor are **historical observation opportunities**, not 303 naturally polymorphic species. The gap-1 class contributes 37/22/33 candidate taxa across the three source-disjoint cohorts. Starting from the original anchored 10/10/20, the minimum hypothetical supplement to reach 30/30/30 is **20/20/10 = 50 additional independently observed photograph slots**, respectively, provided their source-year/month/locality and photo-classifiability constraints are satisfied. It is only a lower bound on future acquisition workload, **not** an observed count of unused images. In particular, Validation needs 20 successful gap-1 completions out of just 22 candidate species (**90.9%**) unless it also recovers gap-2 candidates, making that cohort the immediate limiting case. The three cohorts cannot be combined to bypass the per-cohort gate.

Photo-colour outcomes were excluded from this **metadata-only locality selection** after the historically fixed classifiability stage. New historic photographs have not been fetched. The precommitted original sample sufficiency decision remains **HOLD**, because candidate slots have not been filled and repeated photograph neighbourhoods do not establish common genotypes, pollinator selection, plant fitness or adaptive maintenance.

### Data collection order, explicitly not a biological confirmation

(1) Retain the three source-disjoint original taxon cohorts, years, months, observation opportunity and anchored source-only localities. (2) First check availability of *distinct unused observation/photo IDs* from other observers for all 92 gap-1 species (not only the 50 easiest hindsight successes). (3) If Validation in particular remains <30, assess the 60 pre-existing gap-2 candidates under an explicitly labelled secondary recruitment plan. (4) Validate actual photo licence, coordinate accuracy, plant/flower visibility, imaging and wild/population identity independently, and do not replace failures based on new colour. (5) Only after measured photos and year-specific climate are separately qualified may an environmental association be attempted, with the original source-dependent design and noncausal limitations preserved. No unexposed 2,000+730-taxon prospective species or photo identities should be reused as this retrospective recruitment pool.
