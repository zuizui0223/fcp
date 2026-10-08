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
