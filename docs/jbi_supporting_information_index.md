# Journal of Biogeography Supporting Information index

| ID | File | Rows | SHA-256 |
|---|---|---:|---|
| S1 | `docs/supporting/jbi_table_s1_spatial_scale_20_specifications.csv` | 20 | `555faa49c7dc8f41190fffe9e68467c8d3ef2118fe1037ef1b07b640fa0b28e3` |
| S2 | `docs/supporting/jbi_table_s2_matched_control_models.csv` | 40 | `af520503bca8d73cce8ad8253d47b8d2239037e602afb6e75e5934d1f6b38ace` |
| S3 | `docs/supporting/jbi_table_s3_source_stratified_robustness.csv` | 6 | `167097445eaaf7808ee6fe0fcaf9ca6afd1120e048b3922f4d7e112b01c12f6f` |
| S4 | `docs/supporting/jbi_table_s4_leave_one_family_out.csv` | 108 | `23648a81736741f1b1b46ba4bf91c886b213ac609d4d18f378f3b5b7d60bc751` |
| S5 | `docs/supporting/jbi_table_s5_occurrence_cloud_alternatives.csv` | 38 | `e84084628750c1d78877a5183b7b753c78bebb0057c67de69f25f35819c1a936` |
| S6 | `docs/supporting/jbi_table_s6_frozen_classification_manifest.csv` | 34 | `e3f51b381ffb10c2ad13a24245c9f3b9801dde7cad5400b4806d438a16a61517` |
| S7 | `docs/supporting/jbi_table_s7_classification_correction_log.csv` | 0 | `e5958dc538aed5da964d0e0ea27bdde7e1a7fa1e9aafdca0cb199ce93a2aed1d` |
| S8 | `docs/supporting/jbi_table_s8_baseline_model_dataset.csv` | 34 | `25c5e11a34f2324f7f96e9002fa6b23f724ec7b8d7163babf12a4fbe0f94cc8d` |
| S9 | `docs/supporting/jbi_table_s9_paginated_quality_filtered_model_dataset.csv` | 34 | `0a701a510ea6e8c24c2dd96984141f20b70ffc880702497df302e9238f1e67a8` |
| S10 | `docs/supporting/jbi_table_s10_gbif_taxon_resolution_audit.csv` | 34 | `f171afb4916e4770475087aeb5ed672248cfc47b5eb1fbb142985ce02154458b` |
| S11 | `docs/supporting/jbi_table_s11_paginated_robustness.csv` | 6 | `e96b7469eac8b8c72b8c184c284c63216a80ee1993346513c0a5c5a05b8d1c48` |
| S12 | `docs/supporting/jbi_table_s12_paginated_leave_one_family_out.csv` | 100 | `e0af7697674bbb47ef25b026c3053f1c15b3774a32e1be41c2613b5f1fa8d3cd` |
| S13 | `docs/supporting/jbi_table_s13_phylogenetic_sensitivity_summary.csv` | 2 | `a00fc52f80bbc23657debea20ef2d670ed7203d4f2897ce98493e69f4d1ce86c` |
| S14 | `docs/supporting/jbi_table_s14_opentree_name_resolution.csv` | 34 | `47c28ec5f6ac5f3f9780424fbc0a76ae986a86bd31e1b3f2b7268deb93b15f5c` |
| S15 | `docs/supporting/jbi_table_s15_phylogenetic_replicates.csv` | 200 | `f1702d60f0493c6e657ab87d1eb55c42da6092b011ee40251e2d2760221f7a15` |
| S16 | `docs/supporting/jbi_table_s16_dated_phylogeny_sensitivity_summary.csv` | 6 | `23c6392e5f7b8da9003174c2b5e7c350fc1aa91c1d48cfb9f83df67dfc237d84` |
| S17 | `docs/supporting/jbi_table_s17_vphylomaker2_species_placement.csv` | 34 | `cc1083e97ed0586fb97d749882286ffd4fd3795f753568378cf0963fd151042d` |

Tables S1–S7 derive from workflow run `29972327794`. Tables S8–S15 derive from successful Open Tree phylogenetic-sensitivity run `30067762848`, artifact `8586932030`, digest `sha256:a3ce368fa0dc42bcc26edfca7f09286a8bfe8b609d1b9e58fc75b6f096baf16f`. Tables S16–S17 derive from the fixed-seed dated-megaphylogeny run `30076757379`, artifact `8590190840`, digest `sha256:8f11f59a12758f67124647f719fcc79532651c0512f9e0c199a6afa80d178a68`.

Additional source-backed files:

- **Appendix S1:** `docs/jbi_literature_search_provenance.md`, reconstructing automated discovery and follow-up/enrichment activity from 16–19 July 2026 and distinguishing it from unrecorded human-review metadata; SHA-256 `c8d347ebc48452d05bf77fe611a09e597310a05e9a259876a69ca6cc5f9e53cd`.
- `docs/supporting/jbi_gbif_paginated_qc.json`: GBIF retrieval and filtering audit.
- `docs/supporting/jbi_gbif_paginated_metrics_qc.json`: climate extraction and PCA audit.
- `docs/supporting/jbi_opentree_induced_topology.tre`: induced Open Tree topology used for the topology-based sensitivity models.
- `docs/supporting/jbi_phylogenetic_sensitivity_manifest.json`: Open Tree model settings, seeds, package versions and artifact provenance.
- `docs/supporting/jbi_dated_phylogeny_s1.tre`, `jbi_dated_phylogeny_s2.tre` and `jbi_dated_phylogeny_s3.tre`: fixed-seed, time-scaled V.PhyloMaker2 trees; SHA-256 values are `d21d02d7ffa7c44fa5c9606b8c3e72427618dd96946cad8bc3ef5a76df74b5cc`, `b58a3d2b0f5c9f9290c04ee4aff0d4cea24d047e6f295854d94f94613766e086` and `44cd8e6d1b63150f75c94d0496afefec5e3bf4be8b1fdb1ed0d697ae7804a2f6`, respectively.
- `docs/supporting/jbi_dated_phylogeny_manifest.json`: fixed seed, backbone, placement audit, package commit and dated-model artifact provenance; SHA-256 `99e2843874ebb3f06fd3039c4a523c2e92334e15913e2c1322cc9fda05616d5a`.

The paginated sample remains capped. The Open Tree models use a synthetic topology with Grafen branch lengths, whereas the V.PhyloMaker2 models inherit time scaling from `GBOTB.extended.LCVP`; six species were inserted under placement assumptions. Permanent repository and GBIF DOIs remain `Not verified`.
