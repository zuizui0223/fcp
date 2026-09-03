# Dynamic TerraClimate H1 — compute-only partition amendment v2

Freeze date: 2026-09-03  
Status: **frozen before any v2 dynamic model output exists**.  
Scope: compute topology and transport robustness only; no scientific definition is changed.

## Reason for the amendment

The first exact dynamic TerraClimate execution, GitHub Actions run `33599249288`, ended `not_evaluable` under the pre-result decision rule because the remote-data extraction did not complete within the six-hour job limit. The frozen 1,320-point key passed, but the extraction job reached only 1,000 / 6,545 point-variable requests before cancellation; the prospective model, validation and result-artifact steps never ran.

No dynamic ecological metric or model result exists from that run. The v1 `not_evaluable` decision therefore reflects incomplete computation/data transport, not biological non-support.

## Scientific objects that remain exactly frozen

The v2 recovery **must not change** any of the following:

- core protocol: `display-core-v6-focal-consistent`;
- climate-eligible species: 66;
- fixed spatial support: exactly 20 deterministic locations per species = 1,320 rows;
- fixed-point key SHA-256: `f538d563a667f16b7707060ce53962c7e37c7f6b70cd6d8ecc94aca6a267c279`;
- spatial selection rule: 95% radial core -> TerraClimate-cell deduplication -> deterministic spherical maximin;
- climate source: TerraClimate v1.1/current aggregated service;
- time interval: 1958–2025;
- scaling/climatology baseline: 1991–2020;
- variables: `tmin`, `tmax`, `ppt`, `def`, `vpd`, with `tmean` derived exactly as in v1;
- seven prospective metrics and their directions;
- five direct H1 temporal metrics;
- response-state definitions C-only / S-only / C+S / unresolved;
- documentation-propensity model;
- stabilized-IPW L2 multinomial design;
- geographic-extent covariate;
- pure C-vs-S family-clustered sensitivity;
- 499 family bootstraps;
- seven-metric Holm family;
- seed `20260827`;
- pre-result H1-dynamic supported / unsupported / not-evaluable decision rule frozen in PR #17 comment `5506045557`.

No additional metric, covariate, source, species, location, time window, threshold or favourable subset may be introduced after the v2 result is observed.

## Compute-only change

The 6,545 frozen point-variable requests are assigned deterministically to **64 disjoint compute partitions** using their canonical sorted request rank modulo 64. Each frozen request therefore belongs to exactly one partition.

A partition retrieves only its assigned point-variable histories from the same TerraClimate NCSS source and writes a long-form transport shard. Partitions may execute concurrently. Retry count and worker concurrency are implementation parameters and may differ from v1 because they do not change the requested coordinate, variable, time interval or returned scientific value.

After all 64 partitions complete, a separate reassembly stage must prove:

1. exactly 64 partition indices are present once each;
2. their request-key union equals the exact frozen 6,545-request universe;
3. no request key occurs in more than one partition;
4. every request has the required monthly history for 1958–2025;
5. no coordinate, species, location, variable, month or value is replaced or imputed;
6. the reconstructed wide fixed-location time series contains all 1,320 frozen species-location rows;
7. the model input is produced only after this exact reassembly passes.

The partition topology cannot be used to select favourable completed shards. **Every partition is required.** Any missing or failed partition makes the recovery `not_evaluable`; partial data never advance to modelling.

## H1-dynamic decision remains unchanged

The direct H1 temporal family remains:

1. `seasonal_centroid_cycle_mean`;
2. `interannual_centroid_drift_mean`;
3. `interannual_overlap_loss_mean`;
4. `annual_hypervolume_log_sd`;
5. `temporal_variance_component`.

The two S-predicted metrics, `spatial_variance_component` and `space_time_variance_ratio`, remain contrast/diagnostic evidence and cannot alone promote H1-dynamic.

`SUPPORTED` still requires the complete strict analysis plus at least one direct temporal metric satisfying the already-frozen conjunction:

- `pure_S_vs_C_OR < 1`;
- seven-metric Holm-adjusted `pure_S_vs_C_cluster_p_holm <= 0.05`;
- family-bootstrap `S_vs_C_OR_ci_high < 1`.

`UNSUPPORTED` requires a fully evaluable exact analysis in which all five direct temporal metrics are estimable and none satisfy that conjunction.

`NOT_EVALUABLE` applies if exact partition coverage/reassembly fails, required TerraClimate histories are incomplete, the model set is incomplete or the prespecified decision outputs are not estimable.

## Claim ceiling

Even a supported v2 result would mean comparative concordance between documented C/S organization and dynamic temporal environmental heterogeneity. It would not establish temporal selection, a temperature-specific causal mechanism, morph-level fitness differences, pollinator turnover or gene-flow causation.

A v2 unsupported result would mean only that the seven frozen TerraClimate space-time metrics did not distinguish C/S organization under this design. It would not establish the absence of temporal or fine-scale environmental selection.
