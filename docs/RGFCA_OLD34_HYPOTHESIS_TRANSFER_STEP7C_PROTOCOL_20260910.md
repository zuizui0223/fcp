# RGFCA Step 7C — transfer of the historical 34-species hypotheses

Date: 2026-09-10 JST
Status: frozen before opening transferred-hypothesis outcomes.

## Aim

Test which hypotheses from the historical 34-species C/S programme can be transferred to the photo-derived RGFCA organization states without silently changing either the old hypothesis or the evidential meaning of the new data.

The old literature states and the new photo states are related but not identical:

- historical `C`: documented within-population/local coexistence;
- historical `S`: documented spatial segregation / among-place differentiation;
- RGFCA `C*`: repeated <=100-km observer-disjoint cross-morph photo co-occurrence;
- RGFCA `S*`: significant morph-label spatial segregation in the fixed photo coordinates.

C* is not proof of population coexistence. S* is not evidence of adaptation or genetic differentiation.

## Transferability registry

| ID | Historical mechanism | Historical observable | RGFCA status | Frozen RGFCA observable / reason |
|---|---|---|---|---|
| H1a | temporal thermal heterogeneity | mean occupied BIO4 | external environmental acquisition required | preserve BIO4 as climatological seasonality only; dynamic climate remains a distinct future test |
| H1b | temporal precipitation heterogeneity | mean occupied BIO15 | external environmental acquisition required | preserve BIO15; no replacement variable chosen after outcome |
| H2 | geographic fragmentation / restricted connectivity | multiscale disconnected-range index | **testable now from frozen photo coordinates** | exact historical index: mean of `1 - largest_component_fraction` over DBSCAN connectivity thresholds 100, 250, 500 km |
| H3a | spatial environmental turnover | long-vs-short geographic-distance climate turnover | external environmental acquisition required | preserve historical PC1-PC3 environmental-distance statistic |
| H3b | regional niche displacement | normalized two-sector climate centroid distance | external environmental acquisition required | preserve two geographic sectors and normalized PC1-PC3 centroid separation |
| H3c | regional niche overlap | Gaussian/Bhattacharyya overlap | external environmental acquisition required | preserve historical Gaussian overlap definition |
| H4 | pollinator turnover | pollinator diversity/turnover by place and season | not identifiable from present RGFCA photos | requires external pollinator data; prior global pollination trait gate is inadequate |
| H5 | opposing mutualist-antagonist selection | interaction-community evidence | not identifiable from present RGFCA photos | requires interaction data |
| H6 | breeding system / clonality | self-compatibility, autonomous selfing, clonality | external trait acquisition required | no outcome-dependent backfill |
| H7 | drift / gene flow | neutral differentiation / dispersal proxies | not identifiable from present RGFCA photos | requires genetics or defensible independent dispersal proxies |
| H8 | pigment-mediated abiotic trade-offs | drought, UV, soil, extremes | external environmental acquisition required | requires explicit environmental variables; do not infer from S* alone |
| H9 | phenological partitioning | morph-specific flowering phenology | **proxy-testable now** | local morph-specific photo-date partitioning among <=100-km pairs; not true flowering phenology |
| H10 | NFDS / overdominance | reward/FDS/genetic architecture evidence | not identifiable from present RGFCA photos | photo frequencies alone cannot diagnose the mechanism |
| H0 | total niche-size null | rarefied PC1-PC3 hull volume | external environmental acquisition required | preserve rarefied climate-space volume as negative-control estimand |

## Frozen analysis family for photo-only transfer

Run discovery and reserve independently on the Step 7A species universe (second morph >=10%, >=40 classifiable photos), using the already frozen photo rows.

### H2 — sampled fragmentation

For each species, compute the historical multiscale index exactly:

1. DBSCAN on coordinates with haversine distance, `min_samples=1`, at 100, 250, 500 km;
2. at each threshold calculate `1 - largest component size / n`;
3. average the three values to obtain `frag_multiscale`.

Primary response: `S_star` (binary). H2 predicts greater fragmentation in S*=1 species.

Because C* itself requires close pairs, H2 is **not** tested against C*; doing so would partially reuse the C* definition and create circularity.

Primary statistic: difference in median `frag_multiscale`, S*=1 minus S*=0. Significance: one-sided label permutation of S* labels, 20,000 permutations.

Span-adjusted robustness: regress ranked fragmentation on ranked log1p(maximum pairwise span), use residual fragmentation, and repeat the same S* median contrast/permutation. This is a robustness diagnostic, not a redefinition of H2.

### H9 — local photo-date partitioning proxy

Parse `observed_on`; require a valid calendar date. Convert to day-of-year on a 365-day circular scale (Feb 29 mapped to Feb 28 equivalent after day 59). Circular date distance is `min(|d1-d2|, 365-|d1-d2|)`.

Within each species use only photo pairs:

- geographic distance <=100 km;
- different observers;
- both dates valid.

Species gate: >=30 eligible local pairs, with >=5 same-morph and >=5 cross-morph pairs.

Statistic:

`phenology_delta_days = median(circular date distance | cross morph) - median(circular date distance | same morph)`.

Positive values mean locally co-occurring different morphs tend to be recorded farther apart in seasonal time than same-morph pairs.

Within-species significance: permute morph labels among the fixed dated photos 999 times and recompute the pair statistic. No date threshold is tuned.

Cross-species transferred test: compare `phenology_delta_days` between C*=1 and C*=0 using a two-sided 20,000-label permutation. Direction is deliberately two-sided because the historical H9 allows phenological partitioning to stabilize local coexistence but also allows geographic phenological turnover to contribute to S; the old programme did not predeclare a single C-vs-S sign for this photo proxy.

Also report Spearman association of `phenology_delta_days` with continuous D and with S*, descriptively.

## Multiplicity and replication

The two photo-transfer families are H2 and H9. Within each tranche apply Holm correction across the two primary transferred tests. A transferred hypothesis is called **two-tranche recurrent** only if the frozen primary test has the predicted/appropriate direction and Holm-adjusted p <=0.05 independently in discovery and reserve.

If one tranche passes and the other does not, report `single_tranche_only`; do not pool species across tranches to rescue significance.

## Environmental follow-on firewall

H0/H1a/H1b/H3a/H3b/H3c/H8 are not opened in this photo-only run. Their variables must be acquired for the RGFCA species/coordinates under a separate source freeze, preserving the historical WorldClim/PC-space estimands where possible. Dynamic year-specific climate is a separate extension and must not be substituted for BIO4/BIO15 in the historical replication family.

## Claim boundary

Step 7C is a transfer test of historical organizational hypotheses into a new photo-derived state system. It is not a confirmatory replication of the original 34-species study because the response definitions differ and the Step 7A states were developed after the historical programme.