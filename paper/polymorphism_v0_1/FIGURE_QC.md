# Figure QC — FCP polymorphism paper v0.1

Date: 2026-09-10 JST

## Status

All four main figures render successfully as both PNG and PDF, pass the workflow file-size/existence checks, and are generated exclusively from the synchronized `paper_numbers.json` and frozen figure-data tables.

## Structural QC

- **Figure 1:** appropriate opening figure. The estimand is introduced before any map or mechanism claim; discovery-frame prevalence is explicitly labelled as in-frame rather than angiosperm-wide.
- **Figure 2:** carries the flagship result. Discovery and species-disjoint reserve are visually separated; the continuous D gradient, threshold contrast, and reserve nuisance controls are not collapsed into one statistic.
- **Figure 3:** carries the strongest reviewer-facing robustness evidence. Technical ROI/flip failure is separated from ambiguous palette observations, and Dmin4/observed/Dmax4 are shown as sensitivity variants rather than alternative biological labels.
- **Figure 4:** correctly remains secondary. The panel combines clustering-gain robustness with direct cross-tranche genus-mean concordance, while displaying the Dmin4 failure rather than hiding it.

## Layout risks identified from renderer audit

The binary image stream is not directly exposed to the current review surface, so this is a renderer/layout audit rather than a pixel-level visual inspection. Two minor risks should be checked when the figures are viewed at final journal column size:

1. Figure 2C currently uses separate below-axis tranche labels in addition to `<10%` / `>=10%` tick labels; replacing these with four self-contained tick labels would be cleaner if crowding is visible.
2. Figure 3C pointwise P-value labels and Figure 4B genus labels reach 5.8 pt in the renderer. They are legible in the full-size source but should be raised to at least ~6.5–7 pt if the journal applies strong down-scaling.

Neither issue affects the scientific content, number synchronization, or current manuscript assembly. They are classified as final-layout polish rather than evidence blockers.

## Claim-to-figure contract

- Figure 1: species-level polymorphism estimand and descriptive distribution only.
- Figure 2: replicated association between polymorphism and within-species geographic organization.
- Figure 3: sampled-opportunity, technical-failure, matched-background, and ambiguity-endpoint robustness of the spatial association.
- Figure 4: qualified genus-level taxonomic clustering; not formal phylogenetic signal.

No main figure is permitted to imply adaptation, a universal climate threshold, a shared global boundary, or that ambiguous observations are known hidden morphs.
