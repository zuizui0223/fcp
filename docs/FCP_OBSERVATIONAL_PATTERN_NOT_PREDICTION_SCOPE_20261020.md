# FCP analysis scope: observational patterns and explanatory ecological mechanisms, NO prediction

Date: 2026-10-20 (JST). Project: FCP within-species flower-colour variation.

## User-directed binding research scope

**Do not run, optimize, or center the manuscript on predictive modelling.** Specifically, no AUC/log-loss/Brier score campaigns, no held-out-region performance chase, no prediction engineering, and no ranking models by cross-validation accuracy. Earlier predictive diagnostics in PR #248 are archival qualifications and should not be presented as new science or as the next action.

The main goal is to discover and test **observed biogeographic patterns and explanatory associations**, not to forecast flower colours.

## Fixed scientific questions

1. **Global within-species patterns**: How are photographically measured colour diversity, recurrent achromatic/chromatic axes, and regional morph-frequency contrasts distributed among taxa and places? Always show ecological-opportunity and morphology-quality denominators.
2. **Distribution history versus spatial distance**: Does modern geographical distance associate with colour turnover, and is independent historical/phylogeographic evidence available? Describe distance as an IBD-*like phenotype* association, not as verified colonisation history or genetic IBD.
3. **Environmental associations, not prediction**: Estimate within-species ecological associations with contemporaneous vs long-term climate, soils and (where independently recorded) pollinator composition. Compare IBE-like environmental variation conditioned on geography, and explicitly retain confounding and accessible-area uncertainty. Association/effect sizes and matched spatial nulls are the deliverables, not predictive accuracy.
4. **Morph maintenance and molecular mechanisms**: Separate the descriptive white/chromatic colour direction from whether white is geographically overrepresented *after its own source-composition baseline*. Examine published independent anthocyanin-pathway/tissue-specific PAL/WAL evidence, and test real fitness/genetic mechanisms only where appropriate direct data exist.
5. **Local coexistence**: Keep photographed co-occurrence distinct from genetically stable within-population polymorphism; respect the negative same-month fine-local null and missing documented S-only cases in the literature/photo spatial crosswalk.

## Hard scientific decisions already established

- The real source contains 447,500 measured photos from 4,475 species; 227,107 are classifiable four-state photo labels.
- The original prospective white/chromatic direction (158 eligible species, W=0.517 vs structured-null 0.457, p=.001) remains supported and was NOT modified by subsequent exploration.
- Geographic photo-label sorting is reproducible but modest and is not automatically selection, migration history, or balancing polymorphism.
- White photo contrasts account for a large *absolute* part of geographic label turnover, but white-specific enrichment after both components' own opportunity baselines is **NOT SUPPORTED** (strict component-normalized p=.220 on n=174).
- Within-species BIO5 has a weak, positive species-equal photographic association in the 1,574-species climate panel, but it is not evidence of direct adaptive temperature selection; phenotype BIO5 IBE-like controls remain technically bounded.
- Current source photos do not identify neutral genetic isolation, LGM range trajectories, allele fitness, molecular causes or true local breeding-population morph frequencies.
- Full paper-level C/S comparative inference cannot be obtained from the present strict literature-photo intersection (only 5 taxa with spatial results; zero S-only).

## Next analysis permitted (read-only by scientific inference; no predictive validation)

Compare geographic colour differentiation to environmental dissimilarity using **identical species and observation units** with effect sizes and environmental-availability controls, and report unresolved confounding transparently. Use matched species/region and observer/season robust nulls, not train/test prediction. Do not call geographic distance a reconstruction of range history. For direct history tests, require independently supplied phylogeographic/genomic or paleodistribution evidence; label missing direct data as unidentifiable rather than substituting an AUC workflow.

## Code/repository boundaries

Work in the current FCP branch without altering the frozen 447,500-image source, 50-km original effects, or prospective W. The historical results and prediction diagnostics from PR #248 remain archived, but **are not the current research target**. Do not merge a branch or overwrite the active New Phytologist manuscript without an appropriate integration review.
