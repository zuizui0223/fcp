# Reviewer action status

## Completed in this branch

### Claim and manuscript safeguards

- Reframed moisture niche breadth as the focal reported association unless prospective documentation proves preregistration.
- Added explicit baseline-unambiguous inclusion and freeze rules.
- Replaced ambiguous `family-adjusted` wording with explicit reporting requirements.
- Added multiplicity reporting requirements for all evaluated climatic niche metrics and thresholds.
- Restricted interpretation of fragmentation and turnover nulls.
- Added a reviewer-preemption checklist and submission blockers.

### Implemented analytical outputs

- Added exact model formula, estimator, covariance method and `statsmodels` version to spatial-scale outputs.
- Added number of species, number of families and within/among class counts.
- Added the log-odds coefficient, family-clustered standard error, odds ratio and 95% confidence interval.
- Added convergence, iteration and fitted-probability diagnostics.
- Retained two-sided permutation inference and leave-one-family-out estimates.
- Added a manuscript-facing model table generated from the same run as the robustness analysis.
- Added a complete metric-by-threshold table covering five metrics and four occurrence thresholds.
- Added a frozen baseline-unambiguous manifest with deterministic ordering and a SHA-256 digest.
- Added an explicit empty correction-log template so later changes cannot be made silently.
- Preserved source title, DOI/OpenAlex identifier, evidence snippet and review note into the spatial-scale dataset.
- Normalized those fields to `evidence_source`, `source_id` and `decision_note` for the frozen audit.

### Verified numerical result

- Baseline-unambiguous set: 34 species, 25 families, 19 within-population and 15 among-population.
- Moisture-breadth OR 0.403, family-clustered 95% CI 0.165–0.985.
- Clustered Wald p = 0.0463; permutation p = 0.0439 from 9,999 valid permutations.
- Leave-one-family-out OR range 0.284–0.456, with the negative direction retained for every family omission.
- Broader evidence estimate: OR 0.566, 95% CI 0.293–1.093; permutation p = 0.0886.

### Manuscript progress

- Added `docs/jbi_manuscript_draft.md` with a complete first-pass Abstract, Introduction, Methods, Results and Discussion.
- Inserted the validated model values and explicit exploratory-status language.
- Marked unresolved citation, climate-resolution and journal-format details rather than inventing them.
- Added strict-versus-enriched and leave-one-family-out SVG figures.

## Validation in progress

The full climatic-niche workflow has been re-triggered after preserving source-level evidence. Completion requires confirming that all 34 strict-set rows in `baseline_unambiguous_frozen_manifest.csv` contain traceable source fields and that model outputs remain numerically unchanged.

## Remaining tasks in this PR

1. Read the new successful workflow artifact and verify source-field completeness for all strict-set rows.
2. Add workflow assertions for non-empty `evidence_source`, `source_id` and `decision_note` once completeness is confirmed.
3. Add literature citations and exact WorldClim variable/resolution details to the manuscript draft.
4. Add final figure captions and supplement cross-references.
5. Perform journal-format and language editing.

## Submission blockers

- The source-audit workflow must pass and the frozen manifest must be traceable row by row.
- Complete citations and exact climate-data documentation must be inserted.
- The manuscript draft must undergo final scientific and journal-format review.
- A permanent archived release and DOI are required before submission.
