# FCP — within-species floral colour variation only (scope after the mois split)

**Date:** 14 October 2026. **Active repository:** [fcp](https://github.com/zuizui0223/fcp) for *intraspecific* photographed floral hue variation. The global *between-species* water–flower-colour environmental source and manuscript now live in [mois](https://github.com/zuizui0223/mois).

## Scientific question and central claims

**Among plant species with high-depth repeat photographs, does photographed intraspecific floral-colour variation show recurrent directions in phenotype space and repeated structure across their geographic distribution?**

The primary empirical question is about **variation within named species**, not a worldwide species-wise rainfall coefficient. The photo-derived trait is measured visible flower colour (four coarse labels and nine continuous palette fractions). It is not individually verified genotypic morph status, pigment biochemistry or relative reproductive fitness.

| FCP biological result | Exact high-depth evidence | Claim and boundary |
|---|---|---|
| Photo-colour diversity can be measured repeatedly | 149,900 fixed photos, 1,499 intentionally sampled species (500 discovery,500 species-disjoint validation,499 prospective); 369/363/377 ≥40-classifiable photographed species | Range-wide photographed four-state Gini–Simpson D; not locally breeding-population allele frequencies |
| Original observer-disjoint H1 reliability | Validation 200 original splits, median Spearman rho=0.789; fresh same-source photo IDs among 136 overlapping species rho=0.968 | Reproducibility under original qualification, not artefact-free or independent-source flora measurement |
| **Recurrent achromatic–chromatic direction** | 152 discovery + 129 validation +158 prospectively confirmed vector species = **439** already evaluated primary W vectors; source-weighted descriptive observed **W≈0.515557** versus weighted original *null means* **0.450362** | 439 is **descriptive** across source-exposed retrospective plus pre-frozen prospective cohorts; it is **not** 439 untouched tests or a newly constructed pooled null |
| **Pre-frozen confirmatory W** | Untouched 158 eligible prospective species, observed W=0.517246 vs original structured-null median=0.457143, p=0.001; strict prospective 86 pass | Genuine prior-frozen H2 species/photo-disjoint confirmation **within the same source system** |
| Spatial allocation of photographed colour variation | Fixed-photo 50-km exact species-composition label-preserving null, local depletion 0.02053/0.01867/0.01468 in 166/181/204 species across three cohorts, each p=.005 | Photographed colours in nearby conspecific records are less mixed than source species-wide sample; does not measure demographic local coexistence |
| IBD-like versus environmental colour differentiation | Phenotypic geographic-distance colour component more robust than BIO5-associated residual; stringent observer/background BIO5 tests unsupported | History/dispersal remains possible, but genetic IBD/IBE and adaptation unmeasured |
| Technical validity | White-colour classifier near-clipping OR=1.444; high-clip sensitivity retains 142/158 vectors (89.9%), just below pre-fixed 90% support | Official high-clip outcome **INDETERMINATE**; does not certify biochemical white |

The source H1/H2 conditions, pre-frozen thresholds and all original historical manuscript results remain **unchanged** on `main`. The new [source arithmetic-only W summary](../results/polymorphism_white_axis_three_cohorts_descriptive_20261014.json) includes all three cohort means *without pooling null p-values*. It is verified by `scripts/analysis/verify_fcp_withinspecies_axis_descriptive_20261014.py` using existing FCP main H2 receipt files.

## Scope separation (important for editing and merging)

- **Inside FCP main manuscript:** image measurement and observer-disjoint validation, species-wide D, prospectively confirmed white–nonwhite W, within-species geographic photo-colour allocation, process limitations and independent field-validation needs.
- **Moved into [mois](https://github.com/zuizui0223/mois):** 42,111-photo-species source universe as *measured interspecific atlas*, global 100,543-photo colour/abiotic joins, main congeneric 20,694-species climate slopes, soil36, 1,280 regression rows, 576 identical-species spatial-support coefficients, strict photo/geography/missingness and NHM alternative-measurement limitations, original macro source code and source-verified reports. The shared 42,111-species candidate **metadata** is still relevant for fcp sampling provenance and must not be deleted.
- **Historical PR #196:** this Draft integrated single-paper branch combined two different inferential levels. Its history is preserved for audit, but the *integrated submission proposal is superseded* by the two-repository decision. It should **not** be merged wholesale into FCP main (its branch is stacked on the macro PR #193/#194 code).
- **No inferred common causal chain:** a global positive among-congeneric white/BIO12 slope in mois does not automatically explain the within-species W direction in FCP. Neither project measures universal pigmentation costs, frequency-dependent morph maintenance, genotype or seed/pollen relative fitness.

## What remains before publication

The frozen [original New Phytologist within-species manuscript](POLYMORPHISM_MANUSCRIPT_NEW_PHYTOLOGIST.md) is the FCP paper, not the combined single manuscript in historical PR #196. Its figure/SI and archive status should be reviewed independently; journal wording must distinguish photographic flower-colour variation from genetically verified flower-colour polymorphism. Do not silently borrow the mois species-wise moisture coefficient as an FCP primary result.

**Read first:** [frozen source claim ledger](POLYMORPHISM_CURRENT_CLAIM_LEDGER_20260918.md), [source lineage](POLYMORPHISM_DATA_LINEAGE_MAP_20260925.md), [original manuscript](POLYMORPHISM_MANUSCRIPT_NEW_PHYTOLOGIST.md). Future design/biology can test pigmentation chemistry, local plant genotype identity and environmental morph-specific fitness, but no current support should be overstated as this causal mechanism.
