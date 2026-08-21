# Canonical figure plan for the JBI manuscript

The figures are selected from the paper question backward, not by choosing the most significant analysis output.

**Generation status:** all four canonical figure products are currently generated from versioned inputs and committed in both PNG and PDF form under `docs/figures/`.

## Selection rule

A main-text figure must answer one of three reader questions:

1. **What is the geographic scope of the comparison?**
2. **What is the primary cross-metric result?**
3. **Can the reader see the 34 species underlying that result?**

Robustness diagnostics that primarily answer reviewer/statistical questions remain in tables or Supporting Information.

## Figure 1 — Geographic context of the 34 focal species

**Input:** broader exact GBIF occurrence subset for the same 34 focal species + frozen category labels.

**Purpose:** establish the global biogeographic context and show that the focal taxa are geographically heterogeneous. The lower species strip shows the number of occupied climate cells in the checksum-locked primary analysis.

**Important boundary:** the mapped 58,455 occurrence records come from the broader exact GBIF citation bundle and are used as geographic context/QC. They are not represented as the exact primary occurrence sample that created the frozen climatic metrics, and they are not morph-labelled records.

## Figure 2 — Five-metric forest plot

**Input:** checksum-locked 34-species frozen dataset; exact production model `among ~ metric_z + effort_z` with family-clustered sandwich uncertainty.

**Purpose:** show the main result in one glance. All five odds ratios lie below one. The figure therefore emphasizes the cross-metric direction, rather than elevating the moisture p-value alone.

This is the central result figure.

## Figure 3 — Raw 34-species climatic metrics

**Input:** checksum-locked 34-species frozen dataset.

**Purpose:** expose the observations behind Figure 2. Every point represents one species; metrics are standardized only for plotting so the five dimensions can be shown on the same axis. The plot is not a new inferential analysis.

This figure lets readers judge overlap, extreme species and the limited sample size directly instead of seeing only model coefficients.

## Supporting Figure S1 — Per-species geographic occurrence context

**Input:** broader exact GBIF occurrence subset for the same 34 species.

**Purpose:** distribution/QC audit. One species is shown per map. These maps are supporting geographic context rather than a morph-specific range analysis.

## What is deliberately not promoted to a main figure

- Open Tree and V.PhyloMaker2 results: important robustness evidence, already summarized numerically; best retained in Supporting Information or a sensitivity forest if requested.
- CR2/Satterthwaite: finite-cluster robustness, not a separate biological result.
- power/precision simulation: design diagnostic, not evidence for the ecological hypothesis.
- collinearity: model diagnostic.
- historical control, fragmentation, paginated-sensitivity and expanded-set analyses: retired from the active paper pipeline.

## Narrative order

`Figure 1: where are the focal species?`

→ `Figure 2: what is the comparative effect across five niche metrics?`

→ `Figure 3: what do the 34 species themselves look like?`

→ `Supporting Figure S1: can each species' occurrence context be audited individually?`

This ordering makes the paper readable as a biogeographic result rather than as a sequence of statistical robustness checks.
