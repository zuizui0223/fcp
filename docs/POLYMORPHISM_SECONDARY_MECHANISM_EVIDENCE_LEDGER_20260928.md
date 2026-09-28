# Secondary mechanism evidence ledger — 2026-09-28

## Purpose

This ledger separates the evidential status of the secondary BIO5 and anthocyanin-pathway results used in the current New Phytologist manuscript. It prevents source-derived molecular context from being reported as a result newly estimated by this study.

## A. BIO5 environmental sorting

**Status:** study-derived secondary empirical result.

**Primary source artifacts**
- `results/polymorphism_white_environment_mechanism_20260925/result.json`
- `results/polymorphism_legacy_white_bio5_replication_20260925/result.json`
- `results/polymorphism_white_environment_observer_sensitivity_20260925/result.json`
- `docs/POLYMORPHISM_SECONDARY_REPLAY_ENVIRONMENT_20260928.md`

**Prospective cohort**
- eligible species = 281;
- median within-species white-minus-nonwhite BIO5 contrast = +0.0690113 SD;
- species-level Wilcoxon p = 0.0118162;
- Holm-adjusted p = 0.0354487;
- species-stratified conditional OR per one within-species SD BIO5 = 1.073475;
- 95% CI = 1.029395–1.119441;
- conditional p = 0.0009188;
- frozen within-cohort BIO5 gate = PASS.

**Observer sensitivities**
- observer-paired: 106 species / 144 paired species-observer strata; median delta = 0.000 SD; Wilcoxon p = 0.484962; conditional OR = 0.786636 (95% CI 0.543940–1.137618), p = 0.202310;
- observer-balanced: 352 species; median delta = +0.054098 SD; Wilcoxon p = 0.074964; sign-test p = 0.048443.
- status: post hoc after the primary environmental result; cannot upgrade or redefine the frozen gate.

**Transport**
- discovery species-level median contrast = -0.006777 SD, p = 0.743252;
- validation species-level median contrast = +0.066165 SD, p = 0.054067;
- frozen joint cross-cohort transport rule = NOT SUPPORTED.

**Allowed claim**
White states were associated with warmer BIO5 environments in the frozen primary prospective analysis, but the association is observer-sensitive and BIO5 is not supported as a common cross-cohort rule.

**Not demonstrated**
Causal heat selection, adaptive causation, or a universal temperature effect.

---

## B. *Silene littorea* molecular anchor

**Status:** source-derived quantitative molecular evidence structured by the current study; not a raw-read reanalysis.

**Machine-readable extraction**
- `results/polymorphism_silene_molecular_anchor_20260928/result.json`.

**Source**
- Casimiro-Soriguer et al. (2016), doi:10.3389/fpls.2016.00204.
- mRNA-seq of nine morph × developmental-stage samples.
- 29 anthocyanin-biosynthetic-pathway-related loci.

**Bud-stage expression**
- dark pink versus white: F3h1 = **49.0×**, p = **0.039**; C4h2 = **36.2×**, p = **0.013**; Myb1a = **5.1×**, p = **0.009**;
- light pink versus white: F3h1 = **42.2×**, p = **0.049**; F3′h = **4.5×**, p = **0.047**;
- F3h1 is the only locus significant in both pigmented-versus-white bud contrasts;
- dark versus light regulatory contrasts also include Myb1a = **4.2×**, p = **0.021** and Myb3 = 0.3×, p = 0.033.

**Sequence evidence**
- 622 SNPs reported across the 29 ABP-related loci;
- F3h1 has zero SNPs in the reported UTR/CDS table;
- nine initially colour-associated *Ans* SNPs are synonymous;
- expanded sequencing of 38 individuals found no SNP that consistently differentiated colour morphs.

**Petal biochemistry**
- cyanidin derivatives are the primary anthocyanins;
- source-reported white-versus-pigmented differences include rutin, quercetin and isovitexin;
- the source authors interpret transcriptomic and biochemical evidence as consistent with a blockage near F3h1, potentially mediated by Myb1a.

**Allowed claim**
One high-frequency *S. littorea* PAL system has quantitative transcriptomic and biochemical evidence consistent with a petal-specific regulatory blockage near F3h1/Myb1a.

**Not demonstrated**
A causal Myb1a mutation, raw-read independent reanalysis, cross-species generality, or molecular validation of H2.

---

## C. PAL/WAL natural-frequency contrast

**Status:** study-derived descriptive reanalysis of published supplementary frequency tables.

**Primary source artifacts**
- `results/polymorphism_silene_decoupling_persistence_20260925/result.json`
- `results/polymorphism_crossspecies_pal_wal_frequency_20260925/result.json`
- `docs/POLYMORPHISM_SILENE_DECOUPLING_EVIDENCE_SYNTHESIS_20260925.md`

**Within Silene littorea**
- PAL positive-frequency range = 8–21%;
- PAL median positive frequency = 15.5%;
- WAL positive-frequency range = 0.05–0.86%;
- WAL median positive frequency = 0.21%;
- every positive PAL frequency exceeds the maximum positive WAL frequency.

**Across the Del Valle et al. literature table**
- PAL systems = 13;
- WAL systems = 13;
- median PAL lower bound = 5%;
- median numeric WAL upper bound = 0.1%;
- largest quantified WAL upper bound = 1.4%;
- 7/13 PAL lower bounds exceed that 1.4% maximum.

**Allowed claim**
Flower-restricted anthocyanin-loss phenotypes reach substantially higher reported natural frequencies than whole-plant anthocyanin-loss phenotypes in the source tables.

**Not demonstrated**
An unbiased cross-species effect size, equal mutation rates, or causal proof that tissue restriction alone raises natural frequency.

### C2. Direct bridge to the FCP image cohorts is not estimable

Machine-readable feasibility receipt:
- `results/polymorphism_pal_wal_h2_overlap_20260928/result.json`.
Taxonomic crosswalk audit:
- `docs/POLYMORPHISM_PAL_WAL_TAXONOMIC_ALIAS_AUDIT_20260928.md`.

A taxonomy-aware cross-reference preserved the 25 source strings but separately resolved current accepted-name matches. The natural-frequency bridge subset yielded:
- discovery D-eligible overlap = 0;
- reserve D-eligible overlap = 4;
- prospective selected overlap = 0;
- PAL natural overlap = 2;
- WAL natural overlap = 2.

The PAL overlaps are *Gymnadenia rhellicani* (D = 0.6515) and *Silene gallica* (D = 0.2344). The WAL overlaps are source *Delphinium nelsonii* resolved to accepted *D. nuttallianum* (D = 0.1659) and source *Silene dioca* resolved as an orthographic candidate to accepted *S. dioica* (D = 0.2055). Source *Mimulus guttatus* resolves to *Erythranthe guttata* (D = 0.0000) but its Table S1 frequency is greenhouse-only and remains excluded from the natural-frequency bridge. Source *Mimulis lewisii* resolves taxonomically to *Erythranthe lewisii* but has no FCP cohort overlap.

Descriptively, the two PAL D values have median 0.443 and the two natural-context WAL values median 0.186. These four values are not an inferential comparison.

**Decision:** no PAL-versus-WAL bridge test to H2 or D is estimable. After taxonomic resolution there are only two PAL and two natural-context WAL overlaps, all in the reserve cohort, with no discovery or prospective overlap. The PAL/WAL evidence therefore remains a bounded external mechanistic clue rather than a direct molecular validation of the recurrent white-axis result.

---

## D. PAL/WAL biochemical phenotype definition

**Status:** peer-reviewed source-derived empirical anchor; not newly recalculated from the public frequency tables.

Del Valle et al. (2019; doi:10.1186/s12870-019-2082-6) used HPLC-DAD-MS^n tissue profiling to distinguish:
- PAL: anthocyanins absent from petals but retained in photosynthetic tissues;
- WAL: anthocyanins absent from petals and photosynthetic tissues.

The current study uses this published biochemical definition to interpret the *S. littorea* frequency reanalysis. That within-species contrast is therefore not based on visual white-flower labels alone. The 13 PAL / 13 WAL cross-system entries come from the source paper's compiled Table S1 classifications and were not independently tissue-profiled by the current study.

**Allowed manuscript role**
Methods/Results phenotype definition and evidential anchor for the PAL/WAL frequency contrast.

**Not allowed**
Calling the HPLC phenotype contrast a newly generated molecular dataset of the current study.

---

## E. F3H / SlMyb1a evidence-status boundary

The quantitative F3h1/Myb1a values now appear in Results only as a structured extraction of the published *S. littorea* experiment. The underlying RNA-seq and HPLC measurements remain source-derived.

**Allowed manuscript role**
A source-derived single-species molecular anchor that is reported quantitatively in Results and linked to the PAL persistence reanalysis.

**Not allowed**
Calling the extraction a new transcriptome analysis, claiming Myb1a is the causal mutation, or treating F3h1/Myb1a as the molecular basis of the cross-species white/nonwhite axis.

---

## Evidence hierarchy used by the manuscript

1. **Study-derived tested result:** BIO5 within-cohort association, observer sensitivities and failed transport.
2. **Source-derived quantitative molecular anchor:** structured extraction of published *S. littorea* RNA-seq/sequence/HPLC results.
3. **Study-derived reanalysis:** PAL/WAL natural-frequency contrast.
4. **Source-derived empirical anchor:** HPLC-DAD-MS^n tissue phenotype definition for *S. littorea* PAL/WAL.
5. **Feasibility audit:** taxonomy-aware PAL/WAL × FCP overlap showing no estimable direct bridge.

This hierarchy permits BIO5, the *S. littorea* molecular anchor and PAL/WAL persistence to appear in Results while preserving the boundary between current-study analyses and source-derived molecular measurements.
