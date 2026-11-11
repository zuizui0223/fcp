# Companion macroecology paper — figure and table implementation plan

The figures below are for **one** macroecology manuscript. Do not recycle the active New Phytologist within-species H1/H2 figures as independent macro evidence. Some preview panels can be plotted from the frozen tables now; maps and mechanism panels require distinct validation stages.

## Main Fig. 1 — Sources, units and inference (central to the paper)

A flow diagram with three independent columns: **photographic trait**, **long-term environment at iNaturalist sites**, **NHM text-derived species colour**. Show 100,543 photo IDs / 42,111 taxon names as the opportunity frame; 45,953 four-hue classified photographs; 45,836 climate26 complete (25,028 taxa); 31,925 soil36 complete (20,786 taxa); the directional eligible **20,694 species / 1,934 genera**; and independent botanical matching **2,476 species /1,144 genera**. Annotate 54,590 original four-classification failures explicitly. The NHM test reuses FCP photo-site climates and is **not** geographically independent. End the diagram before causal selection/genetic maintenance: those have no measured endpoint.

## Main Fig. 2 — Signed colours along climate gradients

For global full-image, genus-equal climate26 source, forest plots of **BIO12, vapour pressure, BIO5, solar** × **white/yellow-orange/red-pink/blue-purple** probability slopes. Axes are percentage points per global source SD, with colour-category label as text rather than pigment assay. CI bars use the original *within-cohort* six-climate-axis×four-photo-hue maxT envelope. Do not annotate individual p-values. All four BIO12 hue slopes algebraically sum to ~0.

**Frozen key BIO12 slopes**: white +2.255 [0.319,4.190]; yellow-orange −1.221 [−2.987,+0.545]; red-pink −0.446 [−1.513,+0.621]; blue-purple −0.588 [−1.429,+0.253]. These are not four independent effects.

## Main Fig. 3 — From global correlation to spatial support

Panel A: world + East BIO12-to-white same-scope intervals for broad geography, strict photographic quality, and within genus×coarse cell. **These are different samples**, explicitly label species and genus counts.  
Panel B: A/B/C matched support stages; same species for B and C. World B=+1.326, C=+0.808; east B=+1.102, C=−2.076.  
Panel C (optional): paired B−C effect intervals, show nominal95 and corrected24-comparison95 side by side. World B−C +0.518 [−1.253,+2.289] corrected; East +3.178 [−0.258,+6.613]. The East nominal CI excludes zero but **corrected does not**.

## Main Fig. 4 — What independent colour measurements do and do not replicate

A source schematic, not a direct four-hue vs binary effect-size forest plot on a common axis. Two response definitions: 4-hue photo composition and species-level NHM pure white vs blue/purple. Show external matched support 58 genera, 670 species and effect sign with genus-equal 95% interval [−15.34,+22.60]. Mark independent phenotype measurement but **shared FCP observation-site climate**. Do not depict natural population morph frequency or climate adaptation.

## Main Tables

**Table 1:** Study populations / photographic and text trait source / unit / inclusion / eligible N / inferential independence. Summarize 42,111/100,543 opportunities, classifiability, climate+soil completeness, global directional fit, 17,183/13,961 predictive test species, matched within-cell, external botanical.

**Table 2:** BIO12-white effects and simultaneous CIs by region/source/image and paired matched comparison. Include an explicit CI-family column and scale (global BIO12 source SD) so comparisons are not misread as identical estimands.

## Supplementary evidence (archive instead of adding 8 new main figures)

1. The full 1,280 directional contrast rows: [FCP #193 CSV](https://github.com/zuizui0223/fcp/blob/analysis/fcp-macro-abiotic-hue-direction-20261013/results/fcp_macro_direction_20261013/direction_ledger.csv).
2. The full 576 matching/spatial comparison rows: [FCP #194 CSV](https://github.com/zuizui0223/fcp/blob/analysis/fcp-m1-matched-spatial-support-20261014/results/fcp_m1_matched_spatial_support_20261014/matched_spatial_support_differences.csv).
3. Seven-block predictive-information Shapley fits (#166) for methods traceability only; not seven causal mechanisms.
4. Photo quality/missingness and taxon-weighting influence (#158, #190), presented as technical limits, not new adaptive findings.
5. Cross-region genus+500km transport (#169), interpretation as limited generalization with environmental support caveat.
6. Separate programme-level RGFCA world-boundary result: context only, no new P_hot map is a data product of this manuscript.

**Drawing rule:** Never generate or label a world flower-colour hotspot map from coordinates as if it were a validated true natural colour distribution, unless the explicit M2 null, measurement and coverage gates have been met.
