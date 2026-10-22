# New Phytologist FCP: source-locked synthesis of historical H1–H3 and full-source P1–P4

**Date:** 22 October 2026. **Article status:** working manuscript in [clean main-integration PR #254](https://github.com/zuizui0223/fcp/pull/254), not submitted. **Evidence class:** source-measured photographic phenotype; historical cohort-specific inference and later exploratory tests are distinct. This document **does not introduce any new statistical fit or change a frozen inference**.

## Central biological question and supported answer

**Question:** Do flowering plants repeatedly exhibit similar *directions* of intraspecific colour variation, and are those differences reproducibly but heterogeneously *allocated* across geographic environments?

**Evidence-bounded answer:** A recurrent white–chromatic photographed colour direction and species-conditioned geographic clustering of flower-colour labels are both observed. A separate older, two-cohort high-depth comparison also found that species with greater photo-colour diversity D tend to have stronger within-species geographic organization. Temperature-related photo-label associations are small and heterogeneous. White-axis alignment W itself predicts little of the variation among species in the amount of geographical clustering or the signed BIO5 response. **This is a pattern of repeatable photographed phenotype geometry and geographic placement, not proof of inherited morph maintenance, adaptation, phenotypic accessibility mutations or a common pigment-cost mechanism.**

## Historical analysis H1–H3: measurement and alternative-explanation guards

| Historical test | Supported result or preserved negative outcome | Proper inferential role |
|---|---|---|
| **Old H1: observer-disjoint reliability of four-colour D** | In the historical reserve/validation cohort, **200** observer-disjoint partitions gave median rank rho **0.789103**, 5th percentile **0.765165**, median Lin CCC **0.854831**. First-frozen rule **passed**. A later stricter deterministic split had rho **0.792728** and **failed** its >=0.80 threshold. The new all-4,475 source-wide deterministic observer split gives rho **0.82259**. | **Measurement/reliability foundation**, not evidence that 227,107 photo labels identify stable genetic morphs. The strict failure must remain visible. |
| **Old H2: recurrent white–chromatic displacement direction** | The earlier **158-species prospective** primary test gave W **0.517246** versus category-preserving null median **0.457143**, P **0.001**. | Historically prospectively frozen support of the **white–chromatic target**. It is preserved in the archive but is not the denominator of the current P1. |
| **Old H3a: broad phylogenetic signal in D** | Fresh reserve cohort **341** species with frozen trees: scenario S1 K **0.0710**, P **0.2716**; S2 K **0.0601**, P **0.4134**; S3 K **0.0708**, P **0.2674**. All fail the prespecified support decision. | **Unsupported in this historical reserve**, not proof of no phylogenetic dependence in the later 4,475-species frame. Earlier genus-level structure and H3a reserve null may coexist. |
| **Old H3b: D versus photographed geographic span** | Discovery n = **369**, rho **+0.1799**, P **0.0009**; fresh reserve n = **363**, rho **−0.0026**, P **0.9586**, rank partial rho **+0.0055**, P **0.9162**. Prespecified replication **failed**. | A *sampled photographic span*, **not true species range size**. One positive discovery correlation does not establish a transferable ecological gradient. |

**Exact older source paths (on main and this integration branch):** `results/polymorphism_h1_observer_disjoint_reliability_20260913/result.json`; `results/polymorphism_h1_observer_disjoint_D_reliability_20260913/result.json`; `results/polymorphism_h2_third_cohort_prospective_white_axis_20260917/result.json`; `archive/fcp_submission_20260925/h3a_signal/result.json`; `archive/fcp_submission_20260925/h3b/result.json`. Historical conclusions are frozen and **not retested on all 4,475 species** by merely including the prior files as references.

## Earlier positive ecological connection distinct from P2

A separately frozen **reporting-only**, matched-source ecological association asked whether high photographed species-wide four-colour diversity **D** co-occurs with strong species-specific **geographic organization**. It is a **correlation between diversity and spatial organization**, not the P2 range-wide-minus-local mismatch statistic and not evidence that spatial organization *causes* greater diversity:

- Discovery **369 species**: partial rho **+0.1266367493**, matched spatial-null P **0.007** after sampled-span and clear technical-failure controls.
- Species-disjoint reserve **363 species**: partial rho **+0.0992877113**, P **0.025** under the same controls.
- Reserve flower-minus-background contrast: partial rho **+0.1162411363**, P **0.010**; the recorded ambiguity-endpoint sensitivity tests were also positive.

**Source:** `results/polymorphism_spatial_organization_clue_20260918/result.json`, itself explicitly a reporting receipt from previously frozen analyses. These are historical high-depth cohorts drawn from the **same photo provider**, not independent-provider or genome-level replication; the corresponding relationship is **not claimed to have been newly re-estimated in all 4,475 taxa**.

## Current Full Paper tests P1–P4: one 4,475-species source registry, different estimability

| Directly tested photographic prediction | Estimable species, observed result and inference |
|---|---|
| **P1, recurrent white–chromatic displacement beyond category construction** | **4,475 species**, W **0.478872** versus the matched 199-world source-cohort-and-colour-label-preserving null median **0.432039**; exploratory one-sided P **0.005**. This is a new, **post hoc** whole-source test, different in selection and null replication from old H2. |
| **P2, species-conditional geographical allocation** | **3,600 species** with at least one classified conspecific pair <=50 km; mean global-minus-local photograph mismatch **+0.016116**, within-species fixed-composition spatial null P **0.005**. The **historical stricter 1,761-species** gate yielded **+0.017614902**, also P **0.005**. This establishes photo-colour spatial structure, not local inherited-morph coexistence. |
| **P3, observed environment association** | **4,112 species** with photographed white/nonwhite contrast; equal-species univariate Firth BIO5 coefficient **+0.021888 log odds per °C** (genus bootstrap 95% **[+0.00846,+0.03489]**). BIO4 and BIO5 alone had positive simultaneous seven-axis location-contrast intervals. Positive **average**, not a shared per-species thermal rule: earlier discovery and validation mean slopes were negative, and spatial transfer/background controls weakened the claim. |
| **P4, does a recurrent phenotype direction imply a common geographical or climatic fate?** | W versus 50-km depletion: **3,600 species**, controlled partial rank r **+0.03375**, genus 95% interval **[−0.00226,+0.06920]**. W versus signed BIO5 slope: **4,112 species**, r **−0.00927 [−0.03874,+0.02151]**. Neither estimates evolutionary independence or proves absence of fitness effects. |

**Source of modern result counts and hashes:** `results/fcp_np_all4475_unified_species_ledger_receipt_20261021.json`, `results/fcp_all4475_white_axis_structured_null_actual_20261021.json`, `results/fcp_np_entire4475_real_photo_source_h1_h5_20261021.json`, `results/fcp_np_all4475_W_geography_climate_fate_20261021.json`, `results/fcp_np_all4475_observer_split_D_and_univariate_BIO5_Firth_20261021.json`, `results/fcp_np_all4475_seven_abiotic_axes_verified_20261021.json`. Each taxon remains in the **4,475-row** registry; P2/P3/P4 inability to estimate within a taxon is **NA, not zero**.

## Interpretation priority for the Full Paper

1. **Foundational and positive:** Old H1 shows limited but meaningful observer-disjoint reliability; old H2 supports historical prospective white-axis direction, now supplemented by P1 over the full source. They are **not two independent worldwide replications**.
2. **Primary ecological positive:** P2 establishes nonrandom geographic photo-colour allocation on the species-wide composition-conditioned null; the older D-versus-spatial-organization association is a **distinct, replicated positive correlate** and should not be lost.
3. **Conditional environmental result:** P3 identifies small mean BIO5/BIO4 photo-climate associations, with clear spatial-transfer, source/observer and missingness limits.
4. **Explanatory boundary:** P4 reveals weak W–fate association; old H3a/H3b provide **negative historical attempts**, **not** proof that shared ancestry or true range size are globally irrelevant.
5. **Unsupported causes:** Genetic drift, dispersal, morph-specific reproductive fitness, floral/whole-plant pigment regulation, selection and physiological cost–benefit remain possible, not quantified.

**Important label collision:** Historical H1/H2/H3a/H3b (above), current directly evaluated **P1–P4**, and the **H1–H5 proposed evolutionary mechanism diagram in Fig. 5** have different meanings. Do not use results from one numbering scheme as though they were direct tests of a different scheme. The Fig. 5 developmental/historical/selection/maintenance stages are **conceptual**, not newly measured causal effects.

**Source-classification and independence boundary:** Of **447,500** originally measured photographic records, only **227,107** have one of the four coarse colour labels; **220,393** do not. The reported D and species-level colour records are conditional on classification opportunity. No native-only/wild-only source filter or individually matched genotype/petal chemistry/fitness information establishes genetically inherited population polymorphism. Thus this paper addresses **photographic intraspecific phenotype organization**, with biological hypotheses and independent validation explicitly deferred.
