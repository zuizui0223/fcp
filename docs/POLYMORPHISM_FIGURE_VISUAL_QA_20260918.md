# Polymorphism main-figure visual QA — 2026-09-18

Status: **REGENERATION REQUIRED — final visual PASS not yet issued**

This QA inspects the actual generated figure artifact rather than inferring quality from plotting code.

## Inspected artifact

- workflow run: `35290714346`
- job: `105432735541`
- artifact ID: `10526330340`
- artifact name: `polymorphism-publication-figures-20260918`
- artifact digest: `sha256:21bb4004fdad68f3869009497c2ea18550d63368924d6820aeddf3fdf8388383`
- generating head: `84c9dc2621734023808663241403e462c1732d18`

The ZIP was opened and the five committed PNG figures were inspected directly.

## Figure-by-figure visual audit

### Figure 1 — FAIL, layout/topology

Panel A was readable.

Panel B had two substantive visual problems:

1. the arrowheads in the linear stack pointed upward, opposite to the intended sampling chronology;
2. the third prospective cohort was displayed downstream of the original D-inference cohort, although it is a separate later branch from the opportunity frame.

Required correction:

- branch the global frame top-down into:
  - original validation lane -> original high-depth source -> D inference;
  - prospective confirmation lane -> frozen third-cohort selection/fresh metadata -> 499 × 100 / 377 evaluable / zero replacements.

### Figure 2 — FAIL, annotation collision

Panel A was readable.

Panel B had overlapping reserve annotations near the rho = 0.80 stress-test floor. The floor legend, reserve rho/CCC label and CI region collided at publication size.

Required correction:

- remove the in-axis floor legend;
- label the floor directly;
- vertically offset discovery/reserve rho/CCC annotations from their CI lines;
- reserve additional y margin.

### Figure 3 — FAIL, legend collision

Panel A was readable.

In Panel B, the encoding legend overlapped the strict-reserve data/annotation region.

Required correction:

- move the observed/null encoding legend outside the plotting area below the right axis.

### Figure 4 — PASS on inspected version

The primary and strict prospective structured-null histograms were readable with no material overlap.

Scientific boundary text was visible:

- species-disjoint prospective confirmation;
- same iNaturalist opportunity universe;
- not independent-source replication.

No change requested from this visual audit.

### Figure 5 — PASS on inspected version

The three H3a scenarios and discovery/reserve H3b contrast were readable.

The sampled-span/non-range boundary was visible.

No change requested from this visual audit.

## Regression contract

The layout findings are encoded in:

- `tests/test_make_polymorphism_manuscript_figures.py`

The generator correction is:

- commit `4f2ebe0d227d9f38bcaf51d88c9f60a54ffea599`
- file: `scripts/analysis/make_polymorphism_manuscript_figures.py`

The corrected manifest must report:

- Figure 1 cohort topology = `global_frame_branches_to_original_and_third_cohort`;
- Figure 1 arrow direction = `top_to_bottom`;
- Figure 2 stress annotations = `offset_no_legend_overlap`;
- Figure 3 legend = `outside_below_axis`.

## Legend/documentation synchronization

The New Phytologist figure legends have been rewritten to match the implemented two-panel figures rather than the earlier aspirational multi-panel plan.

The canonical plan is now:

- `docs/POLYMORPHISM_FIGURE_PLAN_20260918.md`

## Final PASS rule

Do **not** check off submission-size visual QA until:

1. the corrected generator completes successfully;
2. a fresh figure artifact is uploaded;
3. its actual Figure 1–5 PNG/PDF outputs are opened;
4. Figure 1 has correct branch topology and downward arrows;
5. Figure 2 has no stress-test annotation collision;
6. Figure 3 has no legend/data collision;
7. Figures 4–5 remain readable and scientifically unchanged;
8. frozen numerical values and hard nonclaims are unchanged.

Until those conditions are met, the submission checklist item “Figure files checked visually at submission size” remains open.
