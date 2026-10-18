# Independent scoped method review: GIFT trait modifiers of FCP environmental slopes

Date: 2026-10-18. Reviewed `scripts/data/fetch_fcp_gift_modifiers_20261018.py`, `scripts/analysis/run_fcp_gift_modifiers_20261018.py`, cached GIFT 3.2 API pages, and outputs in `results/fcp_gift_modifiers_20261018/`.

## Assessment

No blocking implementation defect was found. The extraction uses unique exact species-name matches, keeps unknown trait records distinct from biological absence, and preserves source identifiers, agreement, API metadata, and page hashes. The analysis estimates the nine declared contrasts consistently, with equal-species means and joint genus-cluster bootstrap resampling of both trait groups.

None of the nine nominal 95% bootstrap intervals excludes zero; none of the approximate nine-comparison family intervals excludes zero. The supported conclusion is that these data do not resolve trait-associated differences in the fitted standardized climate slopes. This is not evidence of equivalence, absence of biological tradeoffs, or absence of effects in particular species. The trait join provides independent broad biological descriptors, not molecular mechanism labels.

## Independent verification

All 72 cached public API page hashes matched the frozen extraction profile. I independently reconstructed the roster join from the 372,098 cached GIFT species rows and reproduced all 4,475 assignments, including 4,334 unique exact-name matches and 141 unmatched names. I then reconstructed each admitted trait mapping directly from its cached API pages using the agreement threshold of .66. All three reconstructed mappings matched the saved join exactly. There were no duplicate work_ID values in the combined species or individual trait tables, and no disagreement between matched GIFT genus and the first word of the FCP species name.

The joined environmental panel contains exactly 1,574 species. All nine analytical group counts and point estimates were independently recomputed from `species_slopes_with_traits.csv` and agreed to numerical precision. The saved source records and join preserve missing values rather than recoding them as an absent trait or reference level.

I separately checked joint bootstrap arithmetic with explicit expanded samples. Under prescribed genus multiplicities, the production function's first-group mean minus reference-group mean agreed exactly with expanding the sampled clusters into their constituent species. Shared genera correctly contribute to both groups under the same draw weight. The production data have zero invalid bootstrap draws for every contrast.

## Coverage and denominators

| Trait | Known in full 4,475 roster | Known in environmental 1,574 panel | Groups actually contrasted | Excluded from contrast |
|---|---:|---:|---|---:|
| Self fertilization | 1,024 | 377 | Present 177; absent 200 | 1,197 |
| Lifecycle | 3,832 | 1,347 | Annual 253; perennial 1,070 | 251 |
| Woodiness | 4,148 | 1,456 | Woody 277; non-woody 1,161 | 136 |

Lifecycle exclusions comprise 227 missing, 22 biennial, and two variable species. Woodiness exclusions comprise 118 missing and 18 variable species. Other categories are excluded explicitly rather than silently collapsed into a reference group. All specified minimum group-size conditions (20 species and five genera per level) are met.

The self-fertilization groups contain 104 and 152 genera, with 240 distinct genera across both groups. Annual/perennial groups contain 150 and 524 genera, with 626 in their union. Woody/non-woody groups contain 163 and 539 genera, with 679 in their union. Union cluster counts correctly differ from sums of the two group counts because some genera contain both levels.

Coverage is selective and strongly differs across traits, especially self-fertilization, for which only 377/1,574 panel species have a usable record. The source-backed coverage table includes cohort breakdowns. Conclusions concern available matched species, not all measured or eligible species. Exact string matching is an auditable operational join, not comprehensive taxonomic synonym reconciliation; unmatched names must remain unmatched in these outputs.

## Main results

Differences below are first level minus reference level, in fitted log odds per species-specific primary climate SD. Intervals are nominal genus-bootstrap 95% percentile intervals.

| Contrast | BIO5 difference [interval] | BIO14 difference [interval] | SRADmean difference [interval] |
|---|---|---|---|
| Self fertilization present − absent | −0.0315 [−0.1624, 0.0879] | 0.0035 [−0.1040, 0.1052] | 0.0610 [−0.0614, 0.1929] |
| Annual − perennial | −0.0198 [−0.1003, 0.0630] | 0.0642 [−0.0115, 0.1421] | 0.0458 [−0.0547, 0.1508] |
| Woody − non-woody | −0.0790 [−0.1677, 0.0074] | 0.0269 [−0.0486, 0.1013] | −0.0052 [−0.0993, 0.0869] |

No selected point estimate should be singled out as a demonstrated modifier. The woody/non-woody BIO5 result and annual/perennial BIO14 result remain compatible with zero even before family adjustment.

## Bootstrap and interpretation limits

Each of the three trait analyses resamples its union of genera with replacement and recomputes the two species-weighted means jointly, including random group denominators. This preserves cluster membership and cross-group covariance from shared genera. It estimates stability across an empirical genus-sampling process; it is not a phylogenetic model and does not control all taxonomic or ecological confounding.

The family comprises the declared three trait contrasts by three climate axes. Nominal family intervals use per-comparison alpha .05/9. Their percentile coverage is approximate. With 5,000 draws, the lower family-tail endpoint is near the 14th ordered draw, so tail Monte Carlo resolution is limited. The current results do not hinge on claiming a near-threshold exclusion of zero; all nine intervals cross zero. No additional multiplicity claim should be made for other FCP tests outside this family.

The bootstrap resamples estimated Firth slopes, without refitting photographs or explicitly propagating fitted covariance matrices. The empirical spread of estimated slopes contains sampling noise, but this is not an errors-in-variables hierarchical analysis separating biological heterogeneity from slope-estimation noise. It also inherits the environmental model's photo-independence assumptions, outcome-dependent eligibility, classification uncertainty, and residual geographic/observer confounding.

Because each slope uses its own species climate SD, group differences can reflect differences in sampled environmental range as well as biological sensitivity. Native SD medians appropriately expose this limitation. For example, annual/perennial median BIO5 SD is 3.192 versus 2.594 degrees, while woody/non-woody median SRADmean SD is 1,054 versus 1,312 in native radiation units. These summaries reveal a scale issue but do not adjust for it. The contrasts should not be described as effects of a common physical temperature, rainfall, or radiation increment.

The three traits are analyzed separately and can be associated with each other, geography, cohort, or other unmodeled features. A self-fertilization record is the GIFT present/absent category; it does not quantify realized selfing rate, pollinator dependence, or population-specific mating system. Woodiness and lifecycle are broad descriptors. None directly measures pigment chemistry, organ-specific pigment control, genetic mechanism, energetic costs, or fitness.

The analysis therefore extends the environmental work with an auditable, independently sourced exploratory trait comparison, while finding no resolved explanation of slope differences from these three broad traits. Molecular and experimental tradeoff questions remain open.
