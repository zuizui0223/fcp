# Journal of Biogeography Supporting Information index

This index maps the manuscript identifiers used in `docs/jbi_manuscript_editorial_revision.md` to committed, source-backed files. All tables were derived from the successful `PR climatic niche breadth` workflow run 57 for commit `5c7bd0f22aef54fee9b690e8697e9d8f580f82b2` (run ID `29972327794`; artifact `climatic-niche-breadth-analysis`).

| Identifier | Repository file | Rows | SHA-256 | Source artifact |
|---|---|---:|---|---|
| Table S1 | `docs/supporting/jbi_table_s1_spatial_scale_20_specifications.csv` | 20 | `555faa49c7dc8f41190fffe9e68467c8d3ef2118fe1037ef1b07b640fa0b28e3` | `climatic_niche_metric_by_threshold_table.csv` |
| Table S2 | `docs/supporting/jbi_table_s2_matched_control_models.csv` | 40 | `af520503bca8d73cce8ad8253d47b8d2239037e602afb6e75e5934d1f6b38ace` | `climatic_niche_conditional_logit.csv` |
| Table S3 | `docs/supporting/jbi_table_s3_source_stratified_robustness.csv` | 6 | `167097445eaaf7808ee6fe0fcaf9ca6afd1120e048b3922f4d7e112b01c12f6f` | `climatic_niche_spatial_scale_robustness.csv` |
| Table S4 | `docs/supporting/jbi_table_s4_leave_one_family_out.csv` | 108 | `23648a81736741f1b1b46ba4bf91c886b213ac609d4d18f378f3b5b7d60bc751` | `climatic_niche_spatial_scale_leave_one_family_out.csv` |
| Table S5 | `docs/supporting/jbi_table_s5_occurrence_cloud_alternatives.csv` | 38 | `e84084628750c1d78877a5183b7b753c78bebb0057c67de69f25f35819c1a936` | `spatial_fragmentation_scale_models.csv + environmental_turnover_robustness_models.csv` |
| Table S6 | `docs/supporting/jbi_table_s6_frozen_classification_manifest.csv` | 34 | `416949addd664d6e89230df00fc1e89adad261b51268f24f60cd42770559e217` | `baseline_unambiguous_frozen_manifest.csv` |
| Table S7 | `docs/supporting/jbi_table_s7_classification_correction_log.csv` | 0 | `e84da7bea829b1c8286a9e7dddf92a8de76d14c7959b1cffbe8c10a68caba` | `baseline_unambiguous_correction_log.csv` |

## Figure mapping

| Identifier | Repository file | Verified content |
|---|---|---|
| Figure 1 | `docs/figures/moisture_effect_strict_vs_enriched.svg` | Baseline-unambiguous OR 0.426 (95% CI 0.184–0.985) and all-evidence OR 0.563 (95% CI 0.292–1.085) |
| Figure 2 | `docs/figures/moisture_leave_one_family_out.svg` | Leave-one-family-out OR range 0.317–0.481 across 25 represented families |

## Interpretation boundaries

- Tables S1–S5 contain observational comparative analyses and do not establish climatic causation, local adaptation or morph-specific tolerance.
- Table S6 is the frozen model-eligible baseline manifest, not the complete set of 111 validated cases.
- Table S7 is empty because no post-freeze corrections were recorded in the source artifact. The pre-freeze classifier correction described in the manuscript occurred before this log was frozen.
- Permanent repository and GBIF derived-dataset DOIs remain `Not verified`.
