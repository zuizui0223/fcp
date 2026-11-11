# Fifty-km FCP photographic geography after removing shared-photo pair endpoints (2026-10-10)

## Source-verified decision

**The main observed geographic colour-concentration pattern survives when no photograph is used in more than one 50-km pair.** This remains a post hoc **visible photograph-colour** result, not evidence of inherited population polymorphism or evolutionary selection. It directly closes the narrow concern that repeated use of an identical photograph within many local edges is *necessary* for the original result.

Exact source-state, image-colour categories and 50-km distances were recovered without alteration from the frozen measured-image sources. All four graph/pseudoreplication tests and original historical numerical baseline checks passed in [GitHub Actions 38036033533](https://github.com/zuizui0223/fcp/actions/runs/38036033533), SHA-verified. The original per-species full execution and source data are archived in artifact 11663587388, ZIP SHA256 `51bb22cb0f04cb59f3822b494f74b3fb5533b6deab95c99e2362c4c5d7efbc32`.

## Primary: nearest-geodesic label-blind disjoint matching, at least 10 pairs

| Cohort | Distinct high-depth species passing disjoint gate | Mean original four-colour local depletion | One-sided original matched-null p | Species bootstrap 95% CI |
|---|---:|---:|---:|---|
| Discovery | 189 | +0.021215 | 0.005 | [+0.009844,+0.032314] |
| Validation | 204 | +0.016800 | 0.005 | [+0.004782,+0.028486] |
| Third | 239 | +0.019393 | 0.005 | [+0.008829,+0.030274] |

These are not the identical species sets as the historical >=30 **all-pairs** source gate; the new >=10 photo-disjoint pairs can include extra high-depth species. That distinction is scientifically consequential.

## Secondary: exactly preserve historical 30-local-edge eligible source species

This is an explicitly **post-exposure denominator-comparability check**, added after realizing the first support gate expanded the species sets. It keeps only the historical 166/181/204 source species first, then requires >=10 nonoverlapping photo pairs.

| Cohort | Original eligible species remaining | Disjoint-local mean depletion | One-sided p | 95% species bootstrap CI |
|---|---:|---:|---:|---|
| Discovery | 165 | +0.020630 | 0.005 | [+0.010015,+0.032018] |
| Validation | 181 | +0.012820 | 0.010 | [+0.001580,+0.025538] |
| Third | 204 | +0.020011 | 0.005 | [+0.008295,+0.032439] |

Under the more restrictive **different-observer** within-pair policy, the corresponding historical species subsets still show:

- Discovery n=165, depletion +0.021523, p=0.005, bootstrap CI [+0.010320,+0.033270]
- Validation n=181, depletion +0.013010, p=0.010, CI [+0.002049,+0.025378]
- Third n=203, depletion +0.019021, p=0.005, CI [+0.008270,+0.030168]

The positive result is therefore not dependent solely on pairing observations by the same photographer.

## Counterexamples and no-overstatement boundaries

The alternative **label-blind fixed-seed random edge-order** matching generally gives positive estimates at 10 disjoint pairs, but **Validation, restricted to historical source species, loses bootstrap-lower-bound support for unrestricted-observer pairs**: n=181, depletion +0.009009, permutation p=0.025, species bootstrap CI **[−0.003465,+0.020369]**. This is a real qualification, not a code failure. The same fixed-random different-observer Validation condition has a weaker but just-positive lower bound ([+0.000149,+0.023822]; p=0.030). A claim that *every* matching strategy is robust would be false.

A requirement of **>=30 nonoverlapping pair sets per species** falls below the minimum 30 species/cohort support gate in all combinations. Do not call the resulting small samples null or decisive. Other thresholds 5, 15 and 20 and both selection policies, including negative/unresolved outcomes, appear in the complete machine-readable receipt.

The greedy algorithm is **maximal**, not necessarily a global maximum-cardinality matching. Each original photo can occur only once **within** a matching, but photographs may still be from the same botanical individual, identical locality, observer, or correlated source site. Same-location photo pseudoreplication, true genetic morph segregation and selective fitness remain unidentified. Reported permutation p values have not been familywise adjusted across 5 support thresholds × 2 match policies × 2 match selection orders. This follow-up is retrospective and cannot become independent confirmation.

## Paper integration rule

Keep the current New Phytologist main-story claim at its existing ceiling: *species-wide variation in photographed floral colour is geographically allocated.* Include this photo-endpoint-disjoint analysis as supplementary observation-support robustness and expose the Validation alternative-matching CI that crosses zero. It strengthens the assertion that the pattern is not merely a graph-edge reuse artefact, but not inherited morph maintenance, adaptation, pollinator-mediated selection or floral biochemistry.

The separate 405-photo expert botanical image-colour review remains completely **unscored by independent reviewers**. Existing H1/H2 primary decision, immutable original image labels, original manuscript main and unopened future 2,000+730 photo-taxa are unchanged.

[Permanent numerical receipt](../results/polymorphism_photo_disjoint_50km_posthoc_20261010/reporting_receipt.json) · [Exact source Actions execution](https://github.com/zuizui0223/fcp/actions/runs/38036033533)
