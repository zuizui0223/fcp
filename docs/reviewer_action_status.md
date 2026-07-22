# Reviewer action status

## Completed in this branch

### Claim and manuscript safeguards

- Reframed moisture niche breadth as the focal reported association unless prospective documentation proves preregistration.
- Added explicit baseline-unambiguous inclusion and freeze rules.
- Replaced ambiguous `family-adjusted` wording with the exact model and covariance specification.
- Added multiplicity reporting requirements for all evaluated climatic niche metrics and thresholds.
- Restricted interpretation of fragmentation and turnover non-significance.
- Added a reviewer-preemption checklist and submission blockers.

### Implemented analytical outputs

- Exported the exact formula `among ~ metric_z + effort_z`, estimator, covariance method and `statsmodels` version.
- Exported species count, family count and within/among class counts.
- Exported coefficient, family-clustered SE, odds ratio and family-clustered 95% CI.
- Exported convergence, iteration and fitted-probability diagnostics.
- Retained 9,999 two-sided label permutations and leave-one-family-out refits.
- Generated a manuscript-facing model table from the same run as the robustness analysis.
- Generated the complete five-metric by four-threshold result matrix.
- Generated a deterministically ordered frozen baseline manifest with SHA-256 digest.
- Generated an explicit correction-log template.

### Real-pipeline validation

GitHub Actions run `29915531122` completed successfully on the literature, GBIF and WorldClim pipeline. The companion macroecology workflow run `29915531047` also completed successfully.

The verified baseline-unambiguous moisture result is:

- 34 species from 25 families;
- 19 within-population and 15 among-population cases;
- OR 0.403;
- family-clustered 95% CI 0.165–0.985;
- clustered Wald p 0.0463;
- two-sided permutation p 0.0439 from 9,999 valid permutations;
- leave-one-family-out OR range 0.284–0.456, all in the same direction;
- converged in four iterations.

These values now replace the older provisional 35-species / OR 0.368 result in `docs/current_data_research_results.md`.

### Figures

- Added `docs/figures/moisture_effect_strict_vs_enriched.svg`.
- Added `docs/figures/moisture_leave_one_family_out.svg`.

## Remaining tasks

1. Extend the frozen classification manifest so each strict species retains traceable source identifiers and an explicit decision note. The current digest freezes membership and classification, but the model dataset does not yet preserve the complete source-level audit trail.
2. Convert the manuscript outline into complete Introduction, Methods, Results and Discussion prose.
3. Add the final source-audit table and figure captions to the supplement structure.

## Submission blockers

- The frozen manifest must contain traceable source-level evidence before it is treated as the final supplementary classification audit.
- Complete manuscript text remains required before submission.
- Any correction to the frozen set must be recorded in the correction log and must trigger regeneration of the reported values and figures.