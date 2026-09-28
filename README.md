# FCP — within-species flower-colour variation

This repository contains the current **New Phytologist** manuscript and reproducible analysis for a comparative study of within-species flower-colour variation.

The study asks three linked questions:

1. Can repeated community-science photographs recover a reproducible species-level flower-colour phenotype?
2. Does continuous within-species colour displacement repeatedly align along a common achromatic–chromatic axis across species?
3. Is the amount of within-species colour diversity associated with stronger geographic colour organization?

The active manuscript keeps measurement validation, prospective confirmation, geographic structure, and environmental follow-up analyses as distinct inferential layers.

## Start here

1. **Submission manuscript:** [`docs/POLYMORPHISM_MANUSCRIPT_NEW_PHYTOLOGIST.md`](docs/POLYMORPHISM_MANUSCRIPT_NEW_PHYTOLOGIST.md)
2. **Canonical manuscript:** [`docs/POLYMORPHISM_MANUSCRIPT.md`](docs/POLYMORPHISM_MANUSCRIPT.md)
3. **Reproducibility contract:** [`CURRENT_PAPER_REPRODUCIBILITY.md`](CURRENT_PAPER_REPRODUCIBILITY.md)
4. **Claim-to-input lineage:** [`docs/POLYMORPHISM_DATA_LINEAGE_MAP_20260925.md`](docs/POLYMORPHISM_DATA_LINEAGE_MAP_20260925.md)
5. **Frozen claim ledger:** [`docs/POLYMORPHISM_CURRENT_CLAIM_LEDGER_20260918.md`](docs/POLYMORPHISM_CURRENT_CLAIM_LEDGER_20260918.md)
6. **Supporting Information map:** [`docs/POLYMORPHISM_SUPPORTING_INFORMATION_20260918.md`](docs/POLYMORPHISM_SUPPORTING_INFORMATION_20260918.md)
7. **Figure plan:** [`docs/POLYMORPHISM_FIGURE_PLAN_20260918.md`](docs/POLYMORPHISM_FIGURE_PLAN_20260918.md)
8. **Canonical figures:** [`docs/figures/polymorphism_20260918/`](docs/figures/polymorphism_20260918/)
9. **Provenance receipt:** [`archive/fcp_submission_20260925/NP_PROVENANCE_RELEASE_RECEIPT.md`](archive/fcp_submission_20260925/NP_PROVENANCE_RELEASE_RECEIPT.md)

## Current evidence surface

The manuscript uses high-depth, species-level repeated photographs with location-blind colour measurement and observer-disjoint validation. A species-disjoint prospective cohort tests a pre-frozen white-versus-nonwhite colour-space axis under a construction-preserving null. Separate analyses quantify within-species geographic colour organization and bounded environmental associations.

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

## Reproducibility boundary

Repository organization must not alter frozen numerical results, biological decision rules, immutable source commits, or checksums. When prose and machine-readable artifacts differ, the authority order in [`CURRENT_PAPER_REPRODUCIBILITY.md`](CURRENT_PAPER_REPRODUCIBILITY.md) applies.
