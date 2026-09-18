# Polymorphism main-figure visual QA — 2026-09-18

Status: **FINAL PASS — submission-size visual QA complete**

This QA is based on direct inspection of the freshly regenerated figure artifact after all recorded layout corrections. It supersedes the earlier failed visual inspection of artifact `10526330340`.

## Final inspected artifact

- workflow run: `35317018635`
- job: `105510790410`
- artifact ID: `10535687927`
- artifact name: `polymorphism-publication-figures-20260918`
- artifact digest: `sha256:94cbdeafcd6cb04dafa8a40376e2433f96df8b61f3bbda361c496a9081bd21c2`
- generating head: `b394c72b390b0057711fd1b0a1c0ae75bae49bc8`
- committed generated-figure head: `779446c2d2df6d6ade54bb18661c0d31f483a01b`

The artifact ZIP was opened and all five PNG figures were inspected directly at submission-readable scale. The corresponding PDF outputs were present in the same artifact.

## Correction history

The first inspected artifact revealed three presentation problems:

1. Figure 1 had incorrect cohort topology / arrow direction.
2. Figure 2 had a stress-test annotation collision.
3. Figure 3 had a legend/data collision.

Those findings were encoded in `tests/test_make_polymorphism_manuscript_figures.py` and corrected in the figure generator.

A later fresh artifact showed the intended Figure 1 branching and downward arrows, but the two middle lane boxes were still slightly crowded. A new regression contract was therefore added before the final correction:

- RED regression commit: `fe3381e5601a76a35f4e0610053a60699385300b`
- final lane-separation implementation: `b394c72b390b0057711fd1b0a1c0ae75bae49bc8`
- final `figure-tests`: PASS
- final `build-figures`: PASS

The final manifest requires `middle_lane_center_gap_axes >= 0.54`; the implemented lane centers are 0.22 and 0.78, giving a gap of 0.56.

## Figure-by-figure final audit

### Figure 1 — PASS

Panel A is readable.

Panel B now shows the correct sampling topology:

- the 42,111-species opportunity frame branches into two separate lanes;
- the original validation lane leads to the discovery/reserve high-depth source and D inference;
- the prospective-confirmation lane leads to the frozen third-cohort selection/fresh metadata and the 499-species terminal cohort;
- arrows point from source to target, top to bottom;
- the two lane boxes are visually separated with no material overlap at submission size.

The displayed third-cohort denominator remains 499 species x 100 rows, 377 measurement-evaluable species, and zero replacements.

### Figure 2 — PASS

The observer-disjoint H1 panels are readable.

The rho = 0.80 deterministic stress-test floor, confidence interval and rho/CCC annotations no longer collide. The figure continues to distinguish the first-frozen repeated-partition support from the later stricter deterministic stress-test miss.

### Figure 3 — PASS

The legacy H2 localization panels are readable.

The observed/null encoding legend is outside the data region and does not collide with the strict-reserve result. The panel retains the chronology that the named white-versus-nonwhite target was localized after the original broad geometry had been opened.

### Figure 4 — PASS

The prospective primary and strict structured-null panels are readable with no material annotation overlap.

Frozen values remain unchanged:

- primary: 158 species, W = 0.5172457461, p = 0.001;
- strict: 86 species, W = 0.5329282123, p = 0.001;
- measurement-evaluable species = 377;
- verdict = `H2_PROSPECTIVE_WHITE_AXIS_CONFIRMED`.

The figure retains the required boundary language: species-disjoint prospective confirmation within the same iNaturalist opportunity universe, not independent-source replication.

### Figure 5 — PASS

The H3a S1-S3 results and discovery/reserve H3b contrast are readable.

The figure retains the boundary that sampled photographic span is not true biological range size and does not overgeneralize the unsupported broad phylogenetic-signal test.

## Scientific-integrity check

The visual corrections changed layout only.

The final figure manifest continues to report:

- `status = generated_from_frozen_results`;
- `scientific_claims_changed = false`.

No frozen biological value, threshold, null, cohort, verdict or hard nonclaim was changed during figure QA.

## Final decision

All conditions for submission-size visual QA are satisfied.

**Final figure QA verdict: PASS.**

The New Phytologist checklist item “Figure files checked visually at submission size” may be marked complete.
