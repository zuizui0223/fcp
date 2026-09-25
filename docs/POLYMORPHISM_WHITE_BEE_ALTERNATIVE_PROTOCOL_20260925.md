# FCP independent bee-alternative mechanism protocol — 2026-09-25

Status: **prospectively frozen before opening plant-specific bee-interaction breadth from the 2026 curated GloBI v3.1 dataset for the third-cohort FCP species.**

This analysis is separate from the frozen New Phytologist manuscript. It asks whether a broad bee-pollination correlate can explain the measured white/non-white pattern after accounting for the heat association already identified.

## External source

Noori et al. 2026 curated GloBI bee–plant interaction dataset, Zenodo v3.1:
- DOI: 10.5281/zenodo.18303036
- file: `GloBI_Curated.csv`
- expected MD5: `a8e280aad62f2f6d99ff00f84ef30cd1`

The source contains harmonized bee and plant taxonomy and is strongly geographically/taxonomically biased. Missing plant species are therefore **missing coverage**, never zero bee use.

## Biological prediction

Primary hypothesis B1:

> If broad bee-mediated ecological filtering contributes to the repeated white/non-white contrast, plants interacting with a broader set of bee species should have a lower measured white fraction after controlling for interaction-record effort and long-term heat.

Frozen coefficient direction:
- `bee_species_richness`: **negative**.

This is deliberately narrower than “pollinators cause flower colour.” It tests only breadth of documented bee associations.

## FCP outcome source

Third-cohort prospective measurement artifact:
- workflow run `35177668182`;
- only globally classifiable rows in the four frozen biological colour states;
- direct technical seal run `35687421896` is used;
- exact response-blind high-clip set is removed;
- rows lacking near-clip information are excluded.

Species-level outcome:
- `n_white`;
- `n_nonwhite`;
- `white_fraction = n_white / (n_white+n_nonwhite)`;
- empirical logit `log((n_white+0.5)/(n_nonwhite+0.5))`.

A species requires >=40 retained classifiable rows after these technical filters.

## GloBI plant matching and coverage

Only exact string matches to harmonized `plant_species` are used.
No genus-level fallback, fuzzy matching, synonym search, or unmatched-as-zero recoding is allowed after opening.

For each exact matched plant species derive:
- `n_interaction_records`;
- `bee_species_richness = n_distinct(scientificName)`;
- `plant_family` from GloBI.

Primary analysis eligibility:
- >=5 bee–plant interaction records;
- >=2 unique bee species;
- one unambiguous plant family;
- valid FCP outcome and BIO5.

Analysis is `not_estimable` unless:
- >=100 plant species;
- >=20 plant families;
- both white fraction <=0.10 and >=0.90 are represented by >=20 species each (coverage diagnostic only; the continuous outcome remains primary).

## Heat covariate

WorldClim 2.1 10-arc-minute BIO5 is extracted at every eligible FCP photograph coordinate before species aggregation.

Species heat covariate:
- mean BIO5 across retained FCP photographs.

BIO5 is included because its association with white has already been opened. The purpose is to ask whether bee breadth explains additional between-species variation rather than rediscovering the heat result.

## Primary model

Species is the inferential unit.

Fit unweighted OLS:
`
empirical_logit_white ~ z(log1p(bee_species_richness))
                      + z(log1p(n_interaction_records))
                      + z(mean_BIO5)
`

Inference uses plant-family clustered sandwich standard errors.

Primary B1 support requires all:
1. coverage gate passes;
2. coefficient on bee richness < 0;
3. one-sided p for the pre-specified negative direction < 0.05;
4. family-cluster 95% CI excludes 0 in the negative direction.

Because the directional hypothesis is frozen, one-sided p is calculated from the cluster-robust normal statistic. Two-sided p and CI are also reported.

## Effort-residual sensitivity

As a non-rescuing sensitivity:
1. regress `log1p(bee_species_richness)` on `log1p(n_interaction_records)`;
2. use the residual as effort-adjusted bee breadth;
3. fit `empirical_logit_white ~ bee_breadth_residual + z(mean_BIO5)` with family-clustered SE.

This sensitivity cannot rescue a failed primary gate.

## Heat-versus-bee interpretation

- Bee supported after BIO5: bee breadth remains a viable common ecological filter.
- Bee not supported while BIO5 remains positive: broad bee breadth does not explain the heat-associated white sorting under this dataset.
- Neither result establishes absence or presence of non-bee pollinator effects (hawkmoths, birds, flies, beetles, etc.).

## Hard nonclaims

This analysis does not establish:
- local pollinator preference;
- bee dominance in the pollinator community;
- successful pollination from every GloBI record;
- absence of bee effects for uncovered plants;
- effects of moths, birds, bats, flies, or beetles;
- causal selection;
- evolutionary transition direction.

The curated source itself documents strong geographic and taxonomic biases, so the inference is restricted to covered FCP species.
