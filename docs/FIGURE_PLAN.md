# Canonical figure plan for the JBI manuscript

The figures are selected from the paper question backward, not by choosing the most significant analysis output.

**Canonical target set:** five main-text figures plus three Supporting Figures. The PR workflow regenerates PNG and PDF products from versioned inputs and the canonical analysis artifact.

## Selection rule

A main-text figure must answer one of five reader questions:

1. **What is the geographic scope of the comparison?**
2. **What is the primary cross-metric result?**
3. **Can the reader see the 34 species underlying that result?**
4. **Is the common direction concentrated in one represented plant family?**
5. **Does the direction persist when family clustering, finite-cluster correction and phylogenetic structure are handled differently?**

Diagnostics that do not change the biological interpretation remain in Supporting Information. A later display-core-v6 ecological-mechanism screen is retained as a separate supporting hypothesis-discrimination figure because it uses a different sample/state framework and must not be pooled with the frozen 34-species primary comparison.

## Figure 1 — Geographic context of the 34 focal species

**Input:** broader exact GBIF occurrence subset for the same 34 focal species + frozen category labels.

**Purpose:** establish the global biogeographic context and show that the focal taxa are geographically heterogeneous. The lower species strip shows the number of occupied climate cells in the checksum-locked primary analysis.

**Important boundary:** the mapped 58,455 occurrence records come from the broader exact GBIF citation bundle and are used as geographic context/QC. They are not represented as the exact primary occurrence sample that created the frozen climatic metrics, and they are not morph-labelled records.

## Figure 2 — Five-metric forest plot

**Input:** checksum-locked 34-species frozen dataset; exact production model `among ~ metric_z + effort_z` with family-clustered sandwich uncertainty.

**Purpose:** show the main ecological result in one glance. All five odds ratios lie below one. The figure therefore emphasizes the repeated cross-metric direction — geographically structured cases occur toward narrower sampled climatic niches — rather than elevating the moisture p-value alone.

This is the central result figure and the direct visual test of the broad-gradient-sorting versus broad-niche-coexistence expectations.

## Figure 3 — Raw 34-species climatic metrics

**Input:** checksum-locked 34-species frozen dataset.

**Purpose:** expose the observations behind Figure 2. Every point represents one species; metrics are standardized only for plotting so the five dimensions can be shown on the same axis. The plot is not a new inferential analysis.

This figure lets readers judge overlap, extreme species and the limited sample size directly instead of seeing only model coefficients.

## Figure 4 — Leave-one-family-out stability

**Input:** checksum-locked 34-species frozen dataset; 25 family-deletion refits for each of the five metrics.

**Purpose:** show the full leave-one-family-out result instead of reducing it to min–max numbers in Table 2. All 125 deletion estimates remain below OR = 1, making it visually clear that no single represented family creates the shared direction.

**Interpretive boundary:** this is a concentration/influence diagnostic. It is not a substitute for phylogenetic comparative analysis.

## Figure 5 — Inference-method sensitivity

**Input:** canonical workflow outputs for the primary family-clustered models, CR2/Satterthwaite analysis, Open Tree/Grafen models and dated V.PhyloMaker2 S1–S3 analyses.

**Purpose:** make the paper's central uncertainty visible rather than leaving it scattered across Tables 4–5. Point estimates remain below one under every treatment, while confidence intervals broaden under phylogenetic and finite-cluster analyses.

Open Tree values are medians across 100 polytomy resolutions. The dated-phylogeny point is the median across S1–S3 and the plotted interval is the envelope of the three scenario-specific intervals. These are sensitivity treatments of the same data, not independent tests.

## Supporting Figure S1 — Per-species geographic occurrence context

**Input:** broader exact GBIF occurrence subset for the same 34 species.

**Purpose:** distribution/QC audit. One species is shown per map. These maps are supporting geographic context rather than a morph-specific range analysis.

## Supporting Figure S2 — Finite-sample power/precision design diagnostic

**Input:** the canonical 3,000-replicate design simulations across specified odds-ratio scenarios.

**Purpose:** show why effect-direction recovery and crossing a conventional p < 0.05 threshold can behave differently in a 34-species, 25-family design. Both sign-recovery probability and p < 0.05 probability are shown across the same effect-size grid for all five metrics.

**Interpretive boundary:** this is a design/precision diagnostic only. It is not evidence for the ecological hypothesis and is not used as a post-hoc adequacy criterion.

## Supporting Figure S3 — Ecological hypothesis discrimination in display-core-v6

**Input:** `docs/supporting/jbi_space_time_hypothesis_screen_v1.csv` and `docs/supporting/jbi_bio4_robustness_v1.csv`.

**Purpose:** show which simple ecological explanations survive a later, separate 32-species C/S-informative mechanism screen. Plot S-vs-C odds ratios and family-bootstrap intervals for temperature seasonality, precipitation seasonality, multiscale fragmentation, spatial environmental turnover, regional climate-centroid separation, regional Gaussian climate overlap and total climatic hypervolume. A small inset or second panel shows how the BIO4 estimate changes after mean-temperature and/or absolute-latitude controls.

The intended visual result is a **hypothesis filter**, not a second primary forest plot: the initial BIO4 signal is the only resolved candidate and becomes uncertain after broad geographic controls; precipitation seasonality, fragmentation, large-scale spatial climate partitioning and total hypervolume do not provide clean discrimination.

**Interpretive boundary:** this screen is exploratory and uses a different evidence/state framework from the frozen 34-species primary comparison. It must not be pooled numerically with Figure 2, used to retroactively redefine the primary endpoint, or described as confirmatory evidence for temporal climate. The exact dynamic TerraClimate H1 test is reported separately according to its own frozen decision rule.

## What is deliberately not promoted to a main figure

- collinearity diagnostics: compact and already unproblematic (VIF 1.095–1.415), so Table 3 is sufficient;
- the display-core-v6 mechanism screen: ecologically informative but a different later sample, therefore Supporting Figure S3 rather than a primary result figure;
- paginated-sensitivity and unreviewed expanded-set analyses: not primary manuscript evidence;
- any morph-specific climate map: unavailable because GBIF records are not flower-colour-morph labelled;
- a causal mechanism schematic implying moisture selection: not justified by the species-level occupied-climate analysis;
- the terminal 200-species image atlas: its measurement gate was not evaluable before biological coordinate-colour inference, so it has no ecological result to plot.

## Narrative order

`Figure 1: where are the focal species?`

→ `Figure 2: which competing ecological expectation matches the five-metric comparative direction?`

→ `Figure 3: what do the 34 species themselves look like?`

→ `Figure 4: could one plant family be creating the common direction?`

→ `Figure 5: what survives when the inferential structure is changed?`

→ `Supporting Figure S1: can each species' occurrence context be audited individually?`

→ `Supporting Figure S2: what does the finite-sample design imply about sign recovery versus threshold crossing?`

→ `Supporting Figure S3: which simple ecological mechanisms remain plausible after an independent later hypothesis screen?`

This ordering keeps the biological result first, exposes the observations, then addresses robustness, and finally shows ecological mechanism discrimination without conflating the separate inferential lanes.
