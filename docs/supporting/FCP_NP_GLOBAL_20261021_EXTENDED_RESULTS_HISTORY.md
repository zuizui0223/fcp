# New Phytologist global synthesis — extended chronological Results provenance

The following text preserves detailed original Results subsections compressed in the hypothesis-led main paper on 21 October 2026. It is archival Supporting Information from the pre-refocus PR #248 manuscript; no new independent experiment or primary claim is made.

## Original result 1

### H1: species-level polymorphism is reproducible across observer-disjoint photo sets

The first-frozen repeated-partition H1 test supported observer-disjoint reproducibility in the validation cohort. Across 200 partitions, the median number of paired species was 329. Median split Spearman rho was **0.7891**, with a 5th percentile of **0.7652** and a 95th percentile of 0.8109. Median Lin CCC was **0.8548**, and median Spearman-Brown projected reliability was **0.8821**.

The later deterministic stress test retained 363 validation species with zero observer leakage and yielded rho = **0.7927** (bootstrap 95% interval **0.7419–0.8324**) and CCC = **0.8474**. This split missed its deliberately stricter prespecified rho = 0.80 floor by 0.0073. We therefore retain the first-frozen H1 support while explicitly rejecting a claim of near-perfect or split-invariant reliability.

The later fresh-image transport check provided a stronger same-system replication of D. Among **136** species evaluable in both the earlier frozen analysis and the fresh execution, D showed Spearman rho = **0.9681**, Lin CCC = **0.9720**, and a calibration slope of **0.9691**. The median absolute change in D was **0.0148**, while mean signed fresh-minus-prior change was **+0.0032**. Thus the between-species ordering and scale of D transported strongly across a fresh photo set and a separate measurement execution within the same iNaturalist/FCP measurement system.

These results admit D as a reproducible high-depth species phenotype for the subsequent geometry analyses, but they do not estimate global polymorphism prevalence.

## Original result 2

### Finite-sample correction does not alter D-based conclusions

Finite-sample correction had negligible influence on the species-level phenotype. Raw and corrected D rankings were almost identical in discovery (Spearman rho = **0.999977**) and validation (**0.999978**); mean absolute changes were **0.00486** and **0.00469**, respectively. The primary D–spatial result was unchanged after correction: discovery partial rho changed from 0.126637 to **0.126256** with p = **0.007**, and validation from 0.099288 to **0.099561** with p = **0.025**. The validation flower-minus-background sensitivity likewise remained supported (corrected partial rho = **0.116155**, p = **0.010**).

The sampled-span conclusion was also unchanged (validation raw rho = -0.0026, p = 0.9586; corrected rho = **-0.0024**, p = **0.9629**), and corrected validation Blomberg-K effect sizes remained nearly identical to raw values across S1-S3. Thus the main D-based conclusions are not explained by finite-sample bias in the plug-in diversity estimator.

## Original result 3

### H2 discovery and audit: recurrent geometry localizes to white versus nonwhite

The original H2 analysis first established recurrent label-free geometry before naming its biological direction. At the primary 0.10 tier, discovery had 152 vector species with leading-axis concentration lambda1 = **0.541412**; validation had 129 vector species with lambda1 = **0.535045**. Validation mean squared projection onto the frozen discovery axis was **0.524157**, and the independently fitted discovery and validation leading axes had absolute alignment of approximately **0.986**. Discovery concentration and validation frozen-axis transport each exceeded both the isotropic reference and the subsequent coarse-state-preserving structured null (structured-null p = **0.001** in both cases). The strict 0.20 tier showed the same direction of support.

Audit then localized that recurrent component to white versus nonwhite. In the original discovery and validation cohorts, the fixed white-axis statistic exceeded the construction-preserving structured null at both admissibility tiers.

At the primary 0.10 tier, discovery contained 152 eligible species with W = **0.514625** (null median 0.430808, p = **0.001**) and validation contained 129 species with W = **0.514586** (null median 0.466546, p = **0.001**).

At the strict 0.20 tier, discovery contained 75 species with W = **0.542355** (null median 0.443626, p = **0.001**) and validation contained 65 species with W = **0.510517** (null median 0.469943, p = **0.008**).

When the white contrast was projected out, residual directional concentration no longer exceeded the same construction-preserving null. A non-white-only diagnostic likewise did not support a recurrent chromatic hue axis. Thus the original-cohort H2 signal was narrow: the repeatable component was concentrated along an achromatic–chromatic white-versus-nonwhite direction rather than a general shared hue direction.

Because q_white was identified after the original broad H2 geometry had been opened, these results establish target localization but not untouched prospective confirmation.

## Original result 4

### Prospective confirmation measurement gate

The prospective-confirmation run completed all **256 / 256** terminal partitions and produced **49,900** terminal rows from **499** species. All measurement IDs were unique; no species or rows were replaced; image pixels were not persisted.

Of the 49,900 rows, **25,788** were classifiable into the four biological states and 24,112 were nonclassifiable. Using the prespecified minimum of 40 classifiable rows per species, **377** species were measurement-evaluable, exceeding the required minimum of 250. The measurement-support gate therefore passed before H2 was opened.

## Original result 5

### What the prospective H2 contrast confirms, and post hoc validity diagnostics

The structured-null baseline was itself strongly white-axis aligned: the isotropic eight-dimensional expectation is 0.125, whereas the frozen primary structured-null median was **0.457** and observed W was **0.517**. The confirmatory result therefore concerns the **increment above a coarse-state-preserving construction baseline**, not the entire difference from isotropy. Of the 158 primary vector species, **137 (86.7%)** had white as one of their two leading coarse morphs; their mean per-species white-axis contribution was 0.562, compared with 0.224 among the remaining 21 species.

A post-confirmatory null diagnostic re-applied the continuous minor-cluster gate within each of 299 structured-null worlds, starting from the 185 species that passed the coarse-state gate. The null median was **0.447** (95% interval **0.429–0.470**; plus-one p = **0.0033**), lower than the frozen conditional-null median of 0.457. Thus the preregistered null was slightly more conservative with respect to this gate asymmetry.

The direct one-shot digital-highlight control completed all 49,900 frozen reacquisitions with zero acquisition failures and zero source-byte drift. Highlight metrics were available for 44,098 rows. The response-blind q95 rule set the high-clip threshold at near_clip_fraction > **0.4933**, identifying **2,205** rows. In the within-species conditional logistic model, 23,320 originally classifiable rows from 461 species contributed to inference. A one-SD increase in within-species near_clip_fraction increased the odds of frozen white classification by **1.444** (95% CI **1.389–1.502**). The entire interval lies above the predeclared negligible-coupling upper bound of 1.25, providing direct evidence that the measured white state is exposure-coupled.

Removing the response-blind high-clip set did not remove H2 support. The sensitivity retained **142** primary-tier vectors, with W = **0.503**, structured-null median = **0.454**, and upper-tail p = **0.001**. This corresponds to **89.9%** retention of the original 158 vectors, one vector below the prespecified >=90% requirement. The frozen executable checked this retention criterion before its coupling-CI branch and therefore returned the terminal machine state **INDETERMINATE**. The separately frozen prose protocol would classify the coupling interval itself as a FLAGGED condition because its complete 95% CI exceeds 1.25. Because both the prose rule and executable precedence were frozen before reacquisition, we do not recode the terminal machine state after observing the result. Instead, we report the two facts directly: near-clipping is substantially coupled to white classification, while excess H2 alignment remains supported after removal of the response-blind high-clip set.

## Original result 6

### The prospective BIO5 association is observer-sensitive and does not transport as a common rule

The post-confirmatory environmental-filter test retained **281** prospective-cohort species and 12,583 rows after the response-blind high-clip exclusion. BIO5 was the only one of the three frozen environmental predictions to pass its full gate. The median within-species white-minus-nonwhite BIO5 contrast was **+0.0690 SD**; 57.3% of species had a positive contrast. The species-level Wilcoxon p-value was 0.0118 and remained significant after Holm correction across BIO5, BIO14 and solar radiation (**Holm-adjusted p = 0.0354**). In the species-stratified conditional logistic model, a one-SD within-species increase in BIO5 was associated with OR = **1.073** for white classification (95% CI **1.029–1.119**, p = **0.000919**) after continuous near-clip adjustment. BIO14 and mean solar radiation did not pass their frozen gates (Holm-adjusted p = 0.692 for each).

Post hoc observer controls weakened that BIO5 signal. In the strict observer-paired design, **106 species** contributed 144 species-observer strata containing both colour states; the median BIO5 contrast was **0.000 SD** (Wilcoxon p = **0.485**) and the conditional estimate reversed direction without statistical support (OR = **0.787**, 95% CI **0.544–1.138**, p = **0.202**). In the broader observer-balanced design, **352 species** retained a positive median contrast of **+0.0541 SD**, but the Wilcoxon test was not significant (p = **0.0750**; sign test p = **0.0484**). These sensitivities were post hoc and do not overwrite the frozen within-cohort gate, but they show that the magnitude and support of the BIO5 association depend on observer conditioning.

The fixed species-disjoint BIO5 transport test did not reproduce that association across both original high-depth cohorts. In discovery, **271** eligible species had a median contrast of **-0.0068 SD** (Wilcoxon p = **0.743**), and the species-stratified estimate was OR = **1.033** (95% CI 0.990–1.078, p = **0.131**). In validation, **260** species had a directionally concordant median contrast of **+0.0662 SD**, but the species-level test missed the frozen criterion (p = **0.0541**); the row-level model was positive (OR = **1.048**, 95% CI 1.004–1.094, p = **0.0319**). Because both cohorts were required to pass at both inferential levels, the prespecified transport criterion was not met.

Thus the frozen primary prospective analysis provides a within-cohort association between white states and warmer maximum-temperature environments, but post hoc observer controls weaken that association and the fixed transport test does not support a common cross-cohort BIO5 rule (Fig. S9a).

## Original result 7

### Greater D is associated with stronger within-species geographic colour organization

The positive D–spatial association reproduced across the two species-disjoint high-depth cohorts. The raw association was rho = **0.0892133** (p = **0.034**) in discovery and rho = **0.1016008** (p = **0.025**) in validation.

After controlling for sampled geographic span and clear ROI/flip technical-failure rate, the geometry-preserving analysis remained positive in discovery (partial rho = **0.1266367**, p = **0.007**) and validation (partial rho = **0.0992877**, p = **0.025**). In validation, the matched flower-minus-background response was also positive (partial rho = **0.1162411**, p = **0.010**).

The validation result remained supported under exact uniform ambiguity-endpoint completions: primary D_min4 rho = **0.0970781** (p = **0.029**) and D_max4 rho = **0.1252858** (p = **0.008**); flower-minus-background p-values were **0.009** and **0.006**. Species with greater species-wide sampled colour-state diversity therefore tend to show stronger internal geographic organization, although this association neither partitions D into within- versus among-location components nor identifies the causal process that creates or maintains that organization.

## Original result 8

### Post-freeze fourth-cohort replication under the unchanged estimators

After the three-cohort results were frozen, **498 additional species and 49,800 genuinely measured photographs** were processed with the original image and statistical methods. At **50 km (post hoc)**, **190** new species met both unchanged classifiability and local-pair gates, with species-equal mean depletion **0.0172936** (199-permutation upper-tail **p = 0.005**). Across all four measured cohorts, **741** eligible species gave pooled depletion **0.0176372** (p = **0.005**) from **199,700** photographic measurements. The original discovery, validation and prospective-cohort effects reproduced without changes to their numerical values.

In the new fourth cohort, D reproducibility across 200 observer-disjoint splits had median Spearman rho **0.8056** (median 316 paired species; median concordance correlation **0.8441**). The fixed white-axis test gave W **0.5368** versus structured-null median **0.4363** (primary 0.10 gate, n = **134**, p = **0.001**); the strict 0.20 gate gave W **0.5473** versus **0.4581** (n = **67**, p = **0.001**). This is a *post-freeze, taxon-disjoint, same-iNaturalist-source extension*, not independent-source replication or a new untouched prospective H2 confirmation. W measures squared axis alignment, not a directional rise in white morphs. The completed fifth-cohort and full 447,500-photo results are reported below. The four-cohort results are preserved as their chronologically earlier analysis snapshot and are not re-described as the frozen prospective H2 test.

## Original result 9

### Post hoc robust spatial effect distribution under matched nulls

A further species-level audit asked whether the positive *cohort mean* was driven solely by unusually large effects. Of the **741** eligible species across four cohorts, **417** (56.3%) had positive individual depletion, **61** zero, and **263** negative. The species median was **0.005974**, below the species-equal mean **0.017637**. After excluding the highest-effect **149 species** (20%), the raw mean was **−0.006553**. This sign change by itself does not refute spatial sorting: when the **same top-20% exclusion was performed anew in each of the 199 original species-conditional permutation worlds**, the null mean was **−0.018278**, giving exploratory upper p = **0.005** (Holm p = **0.020** within four newly examined robust statistics). The median, fraction positive (observed **0.56275** vs null mean **0.45007**), and top-10% trimmed mean (**0.001911** vs null **−0.011139**) also exceeded their matched nulls (all pooled within-four Holm p = **0.020**). Across individual cohorts, these Holm p-values ranged **0.020–0.030** (Supporting Information S5e). This analysis remains post hoc and does not control all prior exploratory tests.

## Original result 10

### H3a: no detectable broad tree-wide conservation of D

Validation-cohort Blomberg-K tests were nonsignificant under all three frozen tree-placement scenarios:

- S1: K = **0.0710190**, p = **0.2716**;
- S2: K = **0.0601476**, p = **0.4134**;
- S3: K = **0.0707577**, p = **0.2674**.

Pagel's lambda was small in each scenario and its tests against lambda = 0 were also unsupported. Opportunity-adjusted K sensitivities remained nonsignificant. The frozen verdict was `H3A_PHYLOGENETIC_SIGNAL_NOT_SUPPORTED`.

This result closes the tested claim of broad tree-wide signal under the frozen validation-cohort design; it does not imply that phylogeny is irrelevant to flower-colour polymorphism at all evolutionary scales.

## Original result 11

### H3b: discovery span effect collapses in validation

Discovery showed a positive association between D and sampled photographic span (n = 369, rho = **0.1798786**, p = **0.00089996**). The species-disjoint validation cohort did not reproduce that effect (n = 363, rho = **-0.0025855**, p = **0.9586021**). The observer/classifiability-adjusted partial-rank result was likewise near zero (rho = 0.0055187, p = 0.9162042), and S1-S3 rank-PGLS sensitivities were unsupported.

The frozen verdict was `H3B_SAMPLED_SPAN_REPLICATION_NOT_SUPPORTED`.

