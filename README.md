# FCP — within-species flower-colour variation

This repository contains the current **New Phytologist** manuscript and reproducible analysis for a comparative study of within-species flower-colour variation.

The study asks four linked questions:

1. Can repeated community-science photographs recover a reproducible species-wide flower-colour distribution?
2. When a species is colour-variable across its range, does that variation coexist locally or is it partitioned among geographic localities?
3. Does geographic colour turnover contain only distance/history structure, or also a nonredundant environment-associated component?
4. Does continuous within-species colour displacement repeatedly align along a common achromatic–chromatic direction?

The active manuscript keeps measurement validation, geographic allocation of ITV, spatial-process decomposition and prospective phenotype-space confirmation as distinct inferential layers.

## Repository scope

This repository is the paper-facing FCP analysis surface. The generalized trait/interaction **spatiotemporal turnover** programme now lives in [`zuizui0223/turnover`](https://github.com/zuizui0223/turnover) and is intentionally excluded from the active FCP manuscript, workflows, results, tests, and development documents here. Earlier pre-split development remains recoverable from Git history.

## Start here

1. **Manuscript:** [`docs/POLYMORPHISM_MANUSCRIPT_NEW_PHYTOLOGIST.md`](docs/POLYMORPHISM_MANUSCRIPT_NEW_PHYTOLOGIST.md)
2. **Data lineage:** [`docs/POLYMORPHISM_DATA_LINEAGE_MAP_20260925.md`](docs/POLYMORPHISM_DATA_LINEAGE_MAP_20260925.md)
3. **Supporting Information:** [`docs/POLYMORPHISM_SUPPORTING_INFORMATION_20260918.md`](docs/POLYMORPHISM_SUPPORTING_INFORMATION_20260918.md)
4. **Figures:** [`docs/figures/polymorphism_20260918/`](docs/figures/polymorphism_20260918/)
5. **Reproducibility contract:** [`CURRENT_PAPER_REPRODUCIBILITY.md`](CURRENT_PAPER_REPRODUCIBILITY.md)
6. **Frozen claim ledger:** [`docs/POLYMORPHISM_CURRENT_CLAIM_LEDGER_20260918.md`](docs/POLYMORPHISM_CURRENT_CLAIM_LEDGER_20260918.md)
7. **Figure plan:** [`docs/POLYMORPHISM_FIGURE_PLAN_20260918.md`](docs/POLYMORPHISM_FIGURE_PLAN_20260918.md)
8. **Secondary mechanism evidence hierarchy:** [`docs/POLYMORPHISM_SECONDARY_MECHANISM_EVIDENCE_LEDGER_20260928.md`](docs/POLYMORPHISM_SECONDARY_MECHANISM_EVIDENCE_LEDGER_20260928.md)
9. **Frozen provenance release receipt:** [`archive/fcp_submission_20260925/NP_PROVENANCE_RELEASE_RECEIPT.md`](archive/fcp_submission_20260925/NP_PROVENANCE_RELEASE_RECEIPT.md)

## Current evidence surface

The manuscript uses high-depth, species-level repeated photographs with location-blind colour measurement and observer-disjoint validation. Its strongest current ecological pattern is **distributed polymorphism**: at a fixed 50-km scale, nearby conspecific observations are less colour-diverse than expected from each species' exact range-wide colour composition in discovery, validation and the third cohort. This persists after excluding same-observer pairs, removing white records and using continuous nine-colour distances. Continuous colour turnover contains a larger IBD-like component and a smaller BIO5-associated phenotypic IBE-like residual, but the latter fails stricter background/observer controls and is not treated as local adaptation. Separately, a species-disjoint prospective cohort confirms excess alignment with a pre-frozen achromatic–chromatic axis under a construction-preserving null. The prospective cohort remains within the same iNaturalist source and measurement system, so it is not an independent-source replication.

The 42,111-species metadata frame is an **opportunity frame**, not a denominator for global polymorphism prevalence. Claims are controlled by the frozen machine-readable results, protocols, claim ledger, and data-lineage map linked above.

## Reproduce the paper

The supported local entry point is:

```bash
python -m pip install -e .
python -m pip install pytest
python -m pytest \
  tests/test_polymorphism_manuscript_claims.py \
  tests/test_polymorphism_new_phytologist_submission.py \
  tests/test_make_polymorphism_manuscript_figures.py \
  tests/test_repository_layout.py -q
```

GitHub Actions additionally replay the frozen primary and secondary numerical analyses. Exact environments, immutable input hashes, release assets, and replay tolerances are documented in [`CURRENT_PAPER_REPRODUCIBILITY.md`](CURRENT_PAPER_REPRODUCIBILITY.md).

## Repository layout

| Path | Role |
|---|---|
| `docs/POLYMORPHISM_*` | Current manuscript, protocols, claims, lineage and supporting information |
| `results/polymorphism_*` | Frozen/current machine-readable study results |
| `scripts/analysis/` | Study analysis and figure code |
| `tests/test_polymorphism_*` | Claim, figure and submission regression guards |
| `archive/fcp_submission_20260925/` | Immutable provenance inputs, checksums and release receipts |
| `.github/workflows/` | Current manuscript validation and reproducibility workflows |
| `archive/workflows/` | Historical and superseded paper workflows retained for provenance but inactive in Actions |

## Reproducibility boundary

Repository organization must not alter frozen numerical results, biological decision rules, immutable source commits, or checksums. When prose and machine-readable artifacts differ, the authority order in [`CURRENT_PAPER_REPRODUCIBILITY.md`](CURRENT_PAPER_REPRODUCIBILITY.md) applies.
