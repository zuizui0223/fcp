# FCP current measured results and New Phytologist hypothesis paper — 22 October 2026

**Current biological study:** within-species photographed flower-colour variation (FCP), not the separate between-species / seasonal `mois` programme.

## Current status: the 447,500 original photographs have been measured and analysed

The 17 October precompletion ledger is an **immutable interim evidence snapshot**, **not the current analysis status**. A later successful original-method image, statistical and archive workflow completed the source dataset and frozen estimators. It must not be inferred that a green workflow with skipped statistical jobs completed the analysis; the specific genuine science and publication runs are cited below.

| Actual cohort | Species | Original photograph measurement records |
|---|---:|---:|
| Discovery | 500 | 50,000 |
| Validation | 500 | 50,000 |
| Previously designated prospective sample | 499 | 49,900 |
| Fourth cohort | 498 | 49,800 |
| Newly measured fifth cohort | 2,478 | 247,800 |
| **All five** | **4,475** | **447,500** |

The completed source contains **227,107 four-state classifiable photographs**. The remaining **220,393** records are original photographic measurement outcomes but lack an eligible four-colour label; they **cannot** be assumed to be one specific hue or to represent absence of genuine natural-population variation. The 4,730-species opportunity frame includes **252** species with fewer than 100 eligible source photographs and **three** whole-taxon source-identity exclusions. Neither 4,730 nor the earlier 42,111-species opportunity frame is a prevalence denominator for naturally inherited colour polymorphism.

**Exact measured-photo CSV SHA256:** `8e7f10f0427927208485b2da6fa9c46bb681d6f3611e53b8dfdc606cfc8dabe6`.

Primary audited routes:

- [Original new-U100 image ROI workflow #38997169400](https://github.com/zuizui0223/fcp/actions/runs/38997169400), with all 256 original image-worker partitions.
- [Real five-cohort source join, original 50-km, H1/D and W run #39049598596](https://github.com/zuizui0223/fcp/actions/runs/39049598596), with all five required jobs successful.
- [Permanent publication and independent verification #39050465077](https://github.com/zuizui0223/fcp/actions/runs/39050465077).
- [Immutable measured-data and original-inference Release](https://github.com/zuizui0223/fcp/releases/tag/fcp-all4730-original-method-actual-colour-results-20261017).

## Completed original-method ecological results

At the original **post hoc** 50-km scale, the unchanged per-species eligibility rule requires at least 40 classifiable records and at least 30 within-radius photo pairs. **1,761 species** satisfy that rule, with species-equal range-wide minus local coarse-colour mismatch **+0.017614902**, upper P = **0.005** against 199 species-composition- and position-preserving vertex permutations. The original three cohort effects were reproduced unchanged; the fifth cohort separately had **1,020 eligible species**, +0.017598698, P = 0.005. This is **photographic geographic structure**; it does not establish genetic, within-breeding-population polymorphism, natural selection or historical gene flow.

The five-cohort historical-method observer-disjoint repeated-D check and frozen-target W inference also finished. The historical-method W primary eligibility tier has **1,223 species**, W = **0.518795** versus its structured-null median **0.455365**, upper P = **0.001** (999 permutations). This tier is **not** the all-4,475-species W estimator described below; the earlier 158-species prospective analysis is retained as historical test chronology rather than used as the current main-paper denominator.

A stricter exploratory same-year/month, independently observed **5-km** comparison **did not** detect additional local photo-colour clustering above its month-conditioned null: P = **0.300** (194 species with ≥5 eligible pairs) and P = **0.275** (14 species with ≥30 pairs). This null must not be silently dropped or confused with the supported 50-km result. Recurrent co-located photographs can also represent multicoloured flowers, organ/age differences or exposure artifacts.

## Historical positive measurement and diversity–geography evidence retained

Historical **H1** tested four-colour diversity D repeatability on 200 observer-disjoint partitions. The first-frozen validation median rho was **0.789103** and passed; a later stricter single split (rho **0.792728**) did **not** meet the required 0.80 threshold. The earlier **H2** prospective W test (158 species; W **0.517246** vs null median **0.457143**; P **0.001**) remains frozen validation of a selected phenotype direction, not the denominator for P1.

A distinct ecological discovery/validation result links greater photographed species-wide D with greater geographic colour organization after observed-span and technical-failure controls (discovery **369 species**, partial rho **+0.126637**, P **0.007**; validation **363 species**, **+0.099288**, P **0.025**; validation flower-minus-background **+0.116241**, P **0.010**). This relationship was **not rerun on all 4,475 species** and is not equivalent to P2's species-conditional 50-km mismatch.

The older **H3a** broad tree-wide D signal was unsupported in the 341-taxon reserve across three phylogenetic placements (P = 0.2716–0.4134). The older **H3b** photographed sampled-span association was not replicated in its 363-taxon reserve (rho **−0.002586**, P **0.958602**). These non-supports bound the earlier, smaller-sample analyses and **do not demonstrate phylogenetic or true-range-size irrelevance** in the full roster.

Full historical/modern hypothesis crosswalk, exact evidence paths and tested negative outcomes: [H1–H3 × P1–P4 source-locked synthesis](FCP_NP_H1_H3_P1_P4_ECOLOGICAL_SYNTHESIS_20261022.md).

## New Phytologist Full Paper: four full-source photographic hypotheses

**Ecological priority:** The main biological observations are nonrandom 50-km within-species colour geography (P2), a replicated but historically smaller D–spatial-organisation association, and weak heterogeneous BIO4/BIO5 photo-climate associations (P3). Geography depletion and signed BIO5 slope themselves had a weak positive cross-taxon relationship (n = 3,292; partial rank r = +0.04888, genus 95% interval [+0.01580,+0.08160]); this may involve shared photo-source artifacts and is not a causal or fitness selection signal. The recurring white–chromatic direction (P1) is complementary and by itself poorly predicts those ecological outcomes (P4). See the [source-locked synthesis](FCP_NP_H1_H3_P1_P4_ECOLOGICAL_SYNTHESIS_20261022.md).

The active Full Paper source is [`docs/POLYMORPHISM_MANUSCRIPT_NEW_PHYTOLOGIST.md`](POLYMORPHISM_MANUSCRIPT_NEW_PHYTOLOGIST.md) and its audited main-integration change is [PR #254](https://github.com/zuizui0223/fcp/pull/254). The original development draft [PR #252](https://github.com/zuizui0223/fcp/pull/252) is a historically stacked **source branch** and should not be merged in full. The integration copies exactly identified source blobs without importing the entire dependency history. This is a **working hypothesis paper, not a submitted or accepted article**. The common source registry contains 4,475 species; an unobserved locality contrast or white/nonwhite comparison is **NA**, not evidence for an actual zero effect. Its **post hoc** photographic hypotheses are:

| Prediction tested against the complete source roster | Actually estimable subset and observed result | Interpretation |
|---|---|---|
| P1: excess white–chromatic direction after matched coarse-composition opportunities | 4,475 species; W = 0.478872 vs matched-null median 0.432039; 199-world upper P = 0.005 | Repeatable **photographed phenotype-space** direction; not genetic recurrence or directional whitening |
| P2: local conspecific flowers are more similar than fixed-range-colour null predicts | 3,600 species with any ≤50-km classifiable local pair; mean depletion +0.016116, P = 0.005 | Geographic distribution of photographed variants, not inferred local morph maintenance |
| P3: white/nonwhite labels have nonzero location-climate association | 4,112 species; mean univariate Firth BIO5 +0.021888 log odds per °C, genus bootstrap 95% interval [+0.00846,+0.03489] | Small, source-dependent **association**; two earlier cohort means are negative and spatial transfer tests fail |
| P4: stronger white-axis alignment predicts stronger geographic/thermal separation | W vs spatial depletion: n = 3,600, partial rank +0.03375 [−0.00226,+0.06920]; W vs BIO5 slope: n = 4,112, −0.00927 [−0.03874,+0.02151] | No strong *observed* conditional rank relation; confidence intervals include zero, not proof of evolutionary independence |

These are exploratory results from **one photo provider** and must not be represented as preregistered trials or independent-provider replication. The source-verified result files and figures are included in the [clean main-integration PR #254](https://github.com/zuizui0223/fcp/pull/254), with the original source and development chronology preserved in PR #252. The paper's ecological positive finding is geographic structuring; its additional inference is that a recurring photographed colour axis by itself weakly predicts the species-specific strength of that structuring or the signed thermal association.

The older multivariate **1,574-species / 86,114-photo** model's +0.03233 log odds per °C and the expanded **4,112-species univariate** +0.021888 estimate refer to **different estimands and inclusion rules**. Their numerical difference alone is not a biological decline in thermal association. Similarly, the white-versus-chromatic contribution to absolute geographic mismatch cannot be interpreted as uniquely preferential geographic selection after accounting for baseline colour opportunities.

Neither the full photographic atlas nor external gene-expression case studies directly measure matched individual-genotype flower morphs, pigment investment costs, survival, mating or reproductive fitness. Adaptive selection, neutral population history, drift, gene flow and plasticity remain explanations to discriminate with independent evidence.

## Documentation and merge policy

- This file **supersedes status wording**, not source bytes or numerical historical frozen decisions, in [the 17 October precompletion snapshot](FCP_U100_REAL_199700_EMPIRICAL_EVIDENCE_AND_447500_PENDING_20261017.md). Preserve that snapshot unchanged for chronology and its existing frozen regression guards.
- The historical canonical paper and original frozen results remain source-locked in `docs/POLYMORPHISM_MANUSCRIPT.md` and their original provenance paths; the active manuscript is the 4,475-species P1–P4 Full Paper. **The new Full Paper has not been submitted to the journal.**
- PR #252 was too broad for direct main integration. Its paper-specific 45 source files were imported into a **single, clean, main-parent integration commit**, with dedicated old/new manuscript, source-scope and script/figure self-tests in PR #254. Neither numerical source measurements nor the original frozen inference was changed by this documentation or integration.
- The original finite-sample and source-photo classification constraints remain present, and the current paper does **not** estimate worldwide genetic polymorphism prevalence.
