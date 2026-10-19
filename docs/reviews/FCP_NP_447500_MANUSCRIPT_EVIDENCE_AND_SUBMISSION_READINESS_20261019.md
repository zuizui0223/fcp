# New Phytologist manuscript update and evidence verification — 19 October 2026

**Target manuscript:** `docs/POLYMORPHISM_MANUSCRIPT_NEW_PHYTOLOGIST.md` on draft PR #248. **Status:** Updated, reproducibility-checked within the stated scope, *not ready to submit without shortening and source/review packaging*.

## What changed

1. **Summary and Methods:** The five disjoint measured-source cohorts now correctly total **4,475 species and 447,500 authentic measured photo records**, of which 227,107 have usable original four-state labels. Stage 5 is clearly later than the original prospective H2 test; the temporary October 17 “not completed” status was removed.
2. **Results:** The **1,761-species 50-km** species-equal mismatch decrease, the pooled matched permutation distribution and effect robustness are added. Pairwise mismatch partitioning now distinguishes **83.81% raw white–chromatic contribution** from its **79.57% global mismatch baseline**, with a conditional, baseline-normalized component comparison.
3. **Temperature:** The **1,574-species / 86,114-photo / 727-genus** post hoc joint within-species Firth plus two-stage normal-REML model is included. The positive white–BIO5 association is reported as **equal-species log-odds slope +0.03233 per 1°C [genus-bootstrap +0.01786, +0.04728]**, alongside the distinct **REML +0.01708 [modified-HKSJ +0.00897, +0.02518]**. BIO5 is WorldClim maximum temperature of warmest month, *not observed photographic-date temperature*. BIO14 and radiation means did not have intervals excluding zero.
4. **Discussion and Conclusion:** The research question is why photographed conspecific colour variants are **geographically allocated**. Contrasting candidate causes are (i) divergent selection, (ii) history/gene flow/drift, and (iii) plasticity/observer/image classification. The paper does **not** infer genetic polymorphism, allele-specific fitness or pigment production costs from photographic colour frequencies.
5. **Data availability:** Completed original measured-data release and new result folders are recorded. The original frozen 199,700-photo and prospective-test data lineage remains separately identified, and the new additions are explicitly labelled exploratory.

## Numerical claim verification

Values independently read and compared from authentic Git-tracked source CSVs:

| Claim | Observed from source files | Manuscript / guard |
|---|---:|---|
| Pooled mean depletion | 0.0176149021685411 | 0.0176149 |
| Positive fraction | 984 / 1761 = 0.5587734242 | 55.88% |
| Species median depletion | 0.0068712780577 | 0.006871 |
| Upper-20% removed observed / null | −0.0062466977 / −0.0182779327 | Both signs and selection repeated in null |
| White component | 0.0147632414013 | 0.0147632 |
| Chromatic component | 0.0028516607673 | 0.00285166 |
| Normalized white / chromatic values (636 conditioned taxa) | 0.05695979 / 0.06945287 | 0.05696 / 0.06945 |
| Equal-species BIO5 log-odds per °C | +0.0323286448565 | +0.03233 |
| Genus bootstrap interval for above | +0.0178556507 / +0.0472766744 | +0.01786 / +0.04728 |
| Two-stage REML BIO5 per °C | +0.01707632377 | +0.01708 |
| Original SD-unit two-stage REML | +0.0546730256 | +0.05467 |

All 11 quantitative checks passed; 17 expected numeric/text tokens appear in the updated manuscript; historical fixed W and frozen BIO5-transport results remain in place. No interim “447,500-photo analysis pending” statements remain.

## Source artifacts

- `results/fcp_real_447500_robust_matched_null_20261018/group_statistic_summary.csv` — blob `17ea6f0a33103c8a5f5df1142e9d8db7975a9cb1`
- `results/fcp_real_447500_white_axis_decomposition_20261018/component_summary.csv` — blob `eba2f31d4478b4035420b750ae743d9dce884f0a`
- `results/fcp_environment_native_units_20261019/primary_summary.csv` — blob `9022b7809d786f392289d9a3fb76d64410f4f831`
- `results/fcp_real_447500_environment_model_20261018/meta_summary.csv` — blob `0b7cf025d79c006ad3347db7132644603622ac4a`
- `results/fcp_white_axis_exchangeability_20261019/baseline_normalized_spatial.csv` — blob `2fc5e2750e18d9a869ef4e5663798abf239cdc36`

## Journal-fit checks and remaining issues

As of the 19 October check, the manuscript has a four-bullet **183-word Summary** and six alphabetical keywords; 5 main figures plus 1 main table (6 display items); Discussion is approximately **27.5%** of the manuscript's whitespace-delimited main-plus-table count. These match the relevant journal presentation ranges. The **main text plus embedded table is approximately 10,647 words**, exceeding New Phytologist's usual **6,500–7,500-word** Full Paper length by about 3,147 words over the upper target (longer submissions may be considered with a cover-letter justification). Official source: https://nph.onlinelibrary.wiley.com/hub/journal/14698137/about/author-guidelines

**Before submission:**
- Reduce the submission-oriented version to around 7,000–7,500 words, retaining the positive BIO5 effect and the full-data geographic mechanism question. Most reductions should come from repeated old-stage protocol and negative-test details moved to well-indexed SI, *without concealing inferential failures*.
- Materialize and number new Supporting Information Figs S12–S14 and ensure figure labels and legends actually correspond to generated artifacts.
- Freeze PR #248 merged artifacts to an immutable commit/release, and update the figure/data availability references in the submitted file accordingly.
- Verify original photo-white botanical interpretation where possible and maintain classification/exposure qualifications. A strict same-site/year-month comparison was not supported.
- Fill author/affiliation/correspondence/ORCID placeholders; prepare the required 3-question cover letter and line/page-numbered submission PDF.
- Do not claim climate *causes* colour morph differentiation, fitness advantage, repeated genetic white loss, universal pigment maintenance cost or local FCP from photo labels.

## Inference / revision rule

**Main ecological conclusion:** geographically structured within-species photographed flower-colour variation is observed across many plant species. **Secondary supported descriptive association:** white-labelled flowers are slightly more common at warmer BIO5-normal localities within eligible species on average. **Unresolved causes:** local divergent selection versus demographic history/gene flow/drift versus phenotype plasticity/photo formation. These alternatives are hypotheses, not demonstrated mechanisms.

No new raw photograph processing, genetic experiments or model fitting were conducted by the manuscript editing itself.
