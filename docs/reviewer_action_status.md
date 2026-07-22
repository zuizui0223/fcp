# Reviewer action status

## Completed in this branch

- Reframed moisture niche breadth as the focal reported association unless prospective documentation proves preregistration.
- Added explicit baseline-unambiguous inclusion and freeze rules.
- Required source-level audit and correction logging.
- Replaced ambiguous `family-adjusted` wording with explicit reporting requirements.
- Made the strict-set 95% confidence interval a submission blocker.
- Added multiplicity reporting requirements for all evaluated niche metrics.
- Restricted interpretation of fragmentation and turnover nulls.
- Added a reviewer-preemption checklist and submission blockers.

## Next implementation tasks

1. Extract and report the strict-set coefficient and 95% confidence interval from the model artifact.
2. Export the exact model formula, estimator, package versions, family count and class counts into a manuscript-facing table.
3. Create a frozen strict-classification manifest with source IDs and decision notes.
4. Generate a complete metric-by-threshold multiplicity table.
5. Produce strict-versus-enriched and leave-one-family-out figures.
6. Convert the manuscript outline into complete prose.

These tasks should be completed in this PR or in clearly linked follow-up PRs before submission.
