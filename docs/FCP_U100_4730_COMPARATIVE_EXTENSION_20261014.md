# FCP U100 / 4,730-species within-species comparative extension

**Status (2026-10-14): implemented input/analysis bridge with synthetic regression tests; no new 4,730-species colour results.** This is a follow-up to the frozen New Phytologist FCP pipeline. It does not reopen, overwrite or rerun the original discovery/validation/prospective H1/H2 conclusions. The global *between-species* moisture–flower-colour question belongs to `zuizui0223/mois`, not here.

## Scientific target

Compare **properties of variation within each species** across many species, without mistaking the observed range-wide trait for a local natural-population genetic morph frequency.

| Level | Primary outcome | Comparisons | Important limit |
| --- | --- | --- | --- |
| Photographic phenotype geometry | Study-specific fixed `W_i=(u_i · q_white)^2` | Is achromatic–chromatic displacement recurrent beyond coarse-state construction? | W is sign-invariant; no white-origin or pigment-loss inference |
| Species-wide diversity | Four-state `D=1-Σp²`, pairwise-corrected `D_pair` | How much visible colour variation is sampled? | D=0 means uniform *sampled photos*, not a proven monomorphic natural population |
| Geographic allocation | Frozen 50-km `D_pair−D_local` | Is sample-wide diversity depleted within local neighbourhoods? | A 50-km neighbourhood is not a breeding population |
| Environmental sorting | Within-species z-scored `white−nonwhite` BIO5/BIO14/SRAD contrasts | Do colour states occupy different long-term environments within a species? | Observer/exposure/background/missingness and spatial opportunity can create spurious directions |
| Cross-species commonality | Equal-species distributions, D–W–spatial and environment contrast associations | Does a recurrent relationship transport across cohorts and families? | Family leave-one-out is not a phylogenetic model or causal test |

**Frozen geometry:** nine palette coordinates are `[white,yellow,orange,red,pink,magenta,purple,blue,bronze]`; `q_white = normalize([1,−1/8,...,−1/8])`. The previously frozen Hellinger two-means eligibility (second coarse state and minor continuous cluster each >=0.10; strict >=0.20) and within-coarse-state, species-count-preserving 999-draw null remain the authority. **Do not recompute official W using `disttrait`**: documented near-tie normalization differences prevent exact numerical identity.

## Coverage / inference firewall

- **U0=42,111:** outcome-blind species metadata frame; not global colour prevalence.
- **U100=4,730:** species with `after_observer_cap>=100`; opportunity, not a 4,730-species measured colour result.
- **Allocated 1,500:** frozen discovery 500, validation 500, prospectively selected 500; the last contributed 499 measured species. **Actual fixed measured photographs: 149,900, on 1,499 species.** A third-cohort selected taxon that yielded 98 fresh metadata rows is not treated as measured.
- **Unallocated 3,230:** remaining opportunity species before any new opening. Accordingly, 3,231 U100 species have no completed high-depth measurement in the three fixed cohorts. Do not replace or reuse original cohort IDs in a new prospective test.
- **D-evaluable:** 369 discovery, 363 validation, 377 prospective (=1,109 species with >=40 classifiable photos). **Official primary W-eligible = 152/129/158 (=439),** not 4,730.

Input `--u100` must be a *real* identity table derived from the checksum-frozen metadata/capacity census. Columns: `inat_taxon_id,species` (strongly recommended: `family,after_observer_cap`). It can be reconstructed from the provenance snapshot described in `CURRENT_PAPER_REPRODUCIBILITY.md`; the snapshot contains the complete frozen opportunity frame and P100/P500 selection provenance. **Never generate the missing U100 IDs from colour outcomes, inferred names or a hardcoded list.** The script enforces the expected 4,730 count by default and optionally the observer-cap >=100 rule. When those input bytes have not been supplied, outputs remain *not computed*.

## Steps (do not skip evidence gates)

1. **Census/partition before opening photos.** Validate all 4,730 taxon IDs, exact cohort disjointness and immutable SHA256 of the three frozen measured CSVs. Partition species with `FCP_U100_FOLLOWUP_V1|inat_taxon_id|species` SHA256 modulo number of shards. Shards are computation units, *not biological replicates*. No automatic downloading or expanding is triggered.
2. **Colour measurement with frozen method.** For new species, use the study-specific ROI-v4 / fixed nine-palette measurement source and the existing global classifiability/morph assignment. Store exact source-photo IDs, observer identity, quality/near-clipping variables, audit manifests and a content hash. Do not use generic `disttrait` mode refitting as a substitute for official FCP H2 replication. Report `not_measured`, `insufficient_classifiable` and `D_evaluable` separately, including `D=0` if classifiable photos happen to show one state.
3. **Spatial structure and W.** Evaluate 50-km local depletion with the original species-composition-preserving label-permutation procedure, and W with the frozen coarse-state-preserving geometry construction. The accompanying script **only imports** per-species `W_i` and per-species spatial metrics exported from those authoritative study-specific procedures, checking ranges/identity. Do not interpret a species' W_i alone as individually significant: inference requires the cohort-wide structured null.
4. **Environmental join and checks.** Use an **exact `(inat_taxon_id,photo_id)` join** to site-matched WorldClim baseline BIO5/BIO14/mean solar radiation and the **frozen, response-blind** high-clip classification. The comparator retains >=5 white and >=5 nonwhite, within-species standardizes each variable, and records observer-paired diagnostics (minimum 3 observers with both colour states). Missing environment rows are excluded, not imputed as zero. The `near_clip_fraction` variable is required but this bridge does *not* replace the original conditional logistic near-clip adjustment, matched-background falsification or formal detection-bias model.
5. **Compare across species without pseudo-replication.** Report equal-species effect distributions by cohort, D–50-km depletion, D–W, W–depletion and each of those versus signed environment contrasts. Missing W/space/environment values are `not_estimable`, never zero. Provide family leave-one-out ranges when family information permits; full phylogenetic analysis, sampling-opportunity controls and observer/background negative controls are still necessary for a robust ecological mechanism claim.

### Inputs and execution

```bash
python -m fcp_pipeline.u100_comparison \
  --u100 /trusted/snapshot/u100_species.csv \
  --measured discovery=/trusted/global_monte_carlo_measured_photos_v1.csv \
  --measured validation=/trusted/rgfca_reserve_replication_measured_photos_v1.csv \
  --measured prospective=/trusted/polymorphism_h2_third_cohort_measured_photos_v1.csv \
  --expected-u100 4730 --shards 24 \
  --outdir /trusted/output/u100

# Once independent outcome/quality measurement files are available:
python -m fcp_pipeline.u100_comparison \
  --u100 /trusted/snapshot/u100_species.csv \
  --measured expansion=/trusted/new_measured.csv \
  --environment /trusted/photo_environment_with_frozen_clip.csv \
  --w-metrics /trusted/frozen_H2_species_W.csv \
  --spatial-metrics /trusted/frozen_50km_species.csv \
  --outdir /trusted/output/new_cohort

python -m pytest tests/test_u100_comparison.py -q
```

Environment CSV must have `inat_taxon_id,photo_id,bio5,bio14,srad_mean,near_clip_fraction,high_clip_excluded`; `observer_id` in measured rows enables observer-paired estimates. `high_clip_excluded` **must be imported from the previously frozen response-blind clipping decision, not fit on white results**. W CSV needs `inat_taxon_id,cohort,w_alignment_sq,method`, with `method=frozen_study_specific_H2`; spatial CSV needs `inat_taxon_id,cohort,D_pair,local_pairs_50km,depletion_50km`. Exact taxon ID and cohort keys cannot be duplicated across any input.

Outputs: `u100_species_registry.csv` (4,730 species or explicit test expectation, with reason-coded missingness and deterministic shard); `white_environment_species_contrasts.csv` (one row per estimable species × variable), `comparison_summary.json` (coverage, cross-species descriptive statistics, family leave-one-out correlations and SHA256 fingerprints). No image pixels are downloaded or written by this bridge.

## Claim/evidence boundaries and next biological test

The current prospective H2 W result is 158 eligible species, observed W=0.517246 vs structured-null median 0.457143, p=0.001. The distributed pattern at 50 km replicates across three photographed cohorts. By contrast the signed white–BIO5 result, although passing the original third-cohort mechanism gate (281 species; median +0.069 within-species SD), is **not supported after observer-based transport checks**. BIO14 and solar did not pass that gate. These are existing results, **not new estimates for U100**.

The strongest new target is **heterogeneity of ecological sorting**, rather than a global positive white–heat coefficient: test whether variation in `W_i` and in 50-km depletion across species predicts the *magnitude* (not forced sign) and repeatability of within-species temperature/moisture differences. Separate mean environment from spatial environmental turnover; match accessible habitat and observer effort, and compare against species-preserving nulls. Pre-register any untouched 3,230-species environmental transport analysis before opening its colour outcomes. Family leave-one-out is a sensitivity check, not proof against phylogenetic confounding.

**Not allowed:** interpreting 4,730 as measured or polymorphic species; calling photo categories alleles; pooling old retrospective and prospective p-values; calling BIO5 causation/local adaptation; promoting this synthetic-tested code bridge into a 4,730-species empirical finding before real input coverage and independent QC exist. The global interspecific colour × moisture paper remains in `mois`.
