# Calendar-conditioning attenuation audit of geographic flower-colour organization — 2026-10-08

## Immediate decision

**The existing spatial signal survives stronger date-preserving photo-label nulls, but its estimated magnitude is sensitive to the temporal exchangeability definition.** This is a retrospective *null-sensitivity comparison*, not evidence that time causally accounts for a fixed fraction of spatial colour variation.

This note does not change the frozen New Phytologist manuscript, the historical H1/H2 results, or prospective 2,000+730-taxon selection. It accompanies the already source-verified [species-wide space/season result](../results/polymorphism_specieswide_space_season_posthoc_20261008/result.json).

## Numerical readout (same species denominator within each observer policy)

Each effect is the species-equal mean observed local photo-colour depletion minus the corresponding within-species photograph-label permutation-null expectation. The geographic statistic, photo sample and geographically evaluable species denominator do not change between the three columns; *only the admissible label permutations* change.

| Cohort | Pair policy | N geo species | N exchangeable unconditional/month/year×month | Unconditional excess | Month excess | Year×month excess | Attenuation vs unconditional: month / year×month |
|---|---|---:|---|---:|---:|---:|---|
| Discovery | All | 165 | 148 / 140 / 92 | 0.020935 | 0.017455 | 0.012092 | 16.6% / 42.2% |
| Validation | All | 180 | 163 / 152 / 108 | 0.017843 | 0.015287 | 0.008996 | 14.3% / 49.6% |
| Third | All | 203 | 191 / 174 / 122 | 0.015012 | 0.013324 | 0.008761 | 11.2% / 41.6% |
| Discovery | Different observer | 165 | 148 / 140 / 92 | 0.020472 | 0.016038 | 0.011847 | 21.7% / 42.1% |
| Validation | Different observer | 177 | 161 / 150 / 106 | 0.016371 | 0.014641 | 0.009201 | 10.6% / 43.8% |
| Third | Different observer | 202 | 191 / 175 / 122 | 0.015295 | 0.014102 | 0.008610 | 7.8% / 43.7% |

For every year×month result above, the 199-permutation plus-one upper-tail P is **0.005**. Its 95% species-bootstrap intervals for all pairs are **[0.00505, 0.01998]**, **[0.00252, 0.01562]** and **[0.00370, 0.01390]**, respectively. Different-observer-only year×month intervals are **[0.00476, 0.01977]**, **[0.00289, 0.01595]** and **[0.00335, 0.01388]**.

All reductions are calculated as `(effect_unconditional - effect_conditional) / effect_unconditional`; percentages are a *change in null-corrected effect magnitude*, not a variance partition, proportion mediated, causal attribution to climate or estimate of seasonal selection.

## Why the result is informative without establishing mechanism

1. **Conserved phenomenon:** Within-species floral photo states remain geographically clustered even after photo labels are randomized only within the same species' year×month groups. The same-source, species-disjoint cohorts and different-observer sensitivity agree on the positive direction.
2. **Important estimand sensitivity:** The stricter year×month restriction reduces conditional label exchangeability: the identifiable species count in the three all-pair cohorts falls from 148/163/191 without date conditioning to 92/108/122 under year×month. Non-exchangeable species stay in the geographical denominator with exactly zero *identified additional effect*. Neither that zero nor the smaller effect proves absence of biological geographic structure.
3. **Not a temporal climatic test:** The date strata constrain photographic opportunity and the observed colour-label composition. They do not assign site-specific *year-dependent* climatic anomalies, identify independent flowering genotypes, or remove camera exposure, observer–site associations, land use, neutral population history or plasticity.
4. **No new replication claim:** Discovery, validation and third consist of disjoint taxa but share the same iNaturalist/flower-colour measurement framework. Outcome-exposed null comparisons cannot replace a genuinely new prospective experiment.

## Reproducible calculation and fail-closed conditions

Run from repository root:

```bash
python -m pytest tests/test_polymorphism_space_time_null_sensitivity_20261008.py -q
python scripts/analysis/audit_polymorphism_space_time_null_sensitivity_20261008.py \
  --source results/polymorphism_specieswide_space_season_posthoc_20261008/result.json \
  --output /tmp/fcp_space_time_null_sensitivity_20261008.json
```

The audit checks the original source schema and 199-permutation protocol, keeps both observer-pair policies separate, verifies identical geographic species denominators and observed local statistics across nulls, enforces monotonically shrinking *identifiable* strata and reports the source file SHA256. If the original result is edited, the result's byte identity and diagnostics must be reconsidered.

## Next causal-discrimination gate (not measured here)

The next distinct empirical test would require **site-and-year matched environmental anomalies**, not reusing static WorldClim values with a different occurrence year. To ask whether abiotic conditions generate geographic sorting, match occurrences to flowering-month/year temperature, water-deficit, precipitation and radiation anomalies, preserve geography, taxon and observer opportunity, and test whether the *same species and sites* exhibit reproducible phenotype change with anomaly direction beyond site fixed effects. Photo-exposure calibration and organism-level genotype/pigment measurements are necessary before interpretation as adaptive colour-morph maintenance.

**Paper-safe claim:** Repeated, species-balanced photo-colour geography persists under calendar-month and year×month nulls, while its effect size is sensitive to temporal permutation restrictions; this supports nonrandom allocation of photographed states, not a resolved spatial versus temporal selective cause.
