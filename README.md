# FCP — within-species flower-colour variation

**Repository boundary — 14 October 2026:** This **fcp** repository is exclusively the current **within-species flower-colour variation** paper: image measurement validity, D, repeated achromatic–chromatic phenotype-space axis W and geographic distribution within sampled plant species. The **global between-species flower-colour × rainfall/vapour/soil study**, including its original macro scripts, 1,280 direction rows and 576 matched spatial coefficients, has been migrated to [**zuizui0223/mois**](https://github.com/zuizui0223/mois). Historical Git branches and FCP photo-frame provenance remain accessible; the superseded integrated one-paper PR #196 should not be merged into FCP main. [Split scope and all-three-cohort within-species axis note](docs/FCP_WITHIN_SPECIES_SCOPE_AFTER_MOIS_20261014.md).

**Seasonal photo-white × observation-month climate analyses now live in [`zuizui0223/mois`](https://github.com/zuizui0223/mois).** The October 2026 source-only date audit, observed-season matched locality white-rate QC, month-specific TerraClimate matching, and follow-up seasonally balanced photo-selection research were transferred to MOIS and passed its independent CI (see [MOIS seasonal plan](https://github.com/zuizui0223/mois/blob/main/docs/MOIS_SEASONAL_WHITE_MOISTURE_RESEARCH_20261016.md), [merged transfer PR #2](https://github.com/zuizui0223/mois/pull/2)). **FCP remains the upstream owner of original source-photo ID provenance, immutable iNaturalist photo/observer identity, blind ROI-v4/EfficientSAM image measurements, and the separate within-species D/W/distributed-polymorphism paper.** Historical FCP season code/Actions are retained for SHA provenance and must not be interpreted as a competing active MOIS seasonal manuscript. As of the transfer no empirical seasonal climate–white association or fully flowering-verified 3.43M images have been established.

This repository contains the current **New Phytologist** manuscript and reproducible analysis for a comparative study of within-species flower-colour variation.

## Repository split: fcp = within-species variation; mois = global interspecific moisture

**This repository (`fcp`) is the active research home for within-species photographic flower-colour variation.** Its New Phytologist paper tests observer-disjoint reliability of species-wide colour diversity `D`, geographically structured local colour differences, and a previously frozen white–chromatic colour-space direction `W` (158 untouched prospective species; 439 measured vectors across three cohorts are descriptive only). The 42,111-species discovery frame describes sampling opportunity, **not** genetically verified natural-population polymorphism prevalence.

**The separate [`mois` repository](https://github.com/zuizui0223/mois) is the active research home for *between-species* global floral colour × rainfall/atmospheric-moisture analysis.** Its source-linked manuscript, 1,280-row environmental coefficient ledger, 576-row matched spatial-support audit, original analysis code, uncertainty/missingness sensitivity records and global moisture figure were copied from their historical FCP branches. The transfer preserves original Git blobs and scientific negative results; it does not rerun the source-photo models. See [the FCP↔mois split receipt](docs/FCP_MOIS_RESEARCH_SCOPE_SPLIT_20261014.md).

**Publication boundary:** the integrated global-plus-ITV Draft [FCP PR #196](https://github.com/zuizui0223/fcp/pull/196) is historical/superseded, not the active FCP manuscript; the [original within-species New Phytologist manuscript](docs/POLYMORPHISM_MANUSCRIPT_NEW_PHYTOLOGIST.md) remains authoritative. Prior macro development PRs [#193](https://github.com/zuizui0223/fcp/pull/193) and [#194](https://github.com/zuizui0223/fcp/pull/194) remain source-provenance records and are **not** a second active FCP paper. Neither repository has direct pigment chemistry, genotype-frequency or morph-specific reproductive-fitness evidence.

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
