# FCP — within-species flower-colour variation

This repository contains the current **New Phytologist** manuscript and reproducible analysis for a comparative study of within-species flower-colour variation.

## Single New Phytologist integrated submission candidate — 14 October 2026

**One article, not companion papers:** [**Unified Full Paper**](docs/POLYMORPHISM_MANUSCRIPT_NEW_PHYTOLOGIST_UNIFIED_20261014.md) — *Global flower-colour diversity reveals recurrent within-species variation and moisture-associated ecological patterns*.

**Paper sequence — global photographic colour variation → within-species variation structure → testable pigmentation trade-off:** Worldwide original photographs include 42,111 named taxa, 18,457 classifiable single-photo flower colours and 797 differing colour pairs among 3,196 observer-disjoint, both-classifiable conspecific cross-cell pairs (original full pair opportunity 13,416). The 20,694-species congeneric climate panel finds broad moisture-associated photo-white gradients, **but an original 21,448-photo/6,143-species matched comparison does not recover a common positive signed rainfall slope within species**; the picture-quality and spatial cautions remain. **Separate** high-depth sampling of 149,900 photographs from 1,499 species establishes reproducible colour distributions, existing white–chromatic colour-space directions in **439** evaluable species (152 discovery +129 validation +158 prospective), including the **158-species frozen prospective confirmation**, and replicated geographic colour allocation. The 439-species arithmetic summary is descriptive, with no new pooled permutation P-value. The global photo-white fraction (46.95%) is not a genetic 50:50 equilibrium; source tests reject geographic colour-component compensation and do not support middle-rainfall local colour mixing. Pigment physiology, cost and relative fitness remain unmeasured. These datasets are **not nested**. Source-preserving ecological Figures 1–6, the world atlas/cross-cell receipt and supplementary evidence are linked below.

**What the photographs measure:** ROI-v4 flower localization and EfficientSAM segmentation yield a fixed nine-coordinate visible-colour palette, then four coarse photo-hue states when classifiable. This is **not calibrated reflectance/UV, pigment concentration or a genetic local morph census**. Read the [pixel-to-species measurement and identical-photo global/within bridge](docs/POLYMORPHISM_CROSS_SCALE_COLOUR_MEASUREMENT_AND_INTERPRETATION_20261014.md), [original matched-photo coefficients and source-only verification](results/np_unified_20261014/cross_scale_matched_photo_moisture_bridge.json), and [SI S21](docs/POLYMORPHISM_SUPPORTING_INFORMATION_UNIFIED_20261014.md). Original model thresholds, class assignments and prospective H2 remain frozen.

**Version interpretation:** `docs/POLYMORPHISM_MANUSCRIPT_NEW_PHYTOLOGIST.md` and the older frozen claim ledger below remain unchanged as the *historical high-depth scientific baseline*, **not a second planned submission**. The new integrated manuscript is the only journal-facing editorial candidate on this development branch. Main branch and historical provenance receipt are unchanged. Submission still requires author approval/metadata, remaining supplementary visual assets and one new unified permanent source archive.


The study asks three linked questions:

1. How widely do photographed flower colours vary across named plant species and geographic records worldwide, and how is the observed composition associated with moisture?
2. Which repeatable colour-space direction and geographic organization describe **within-species** flower-colour variation under high-density, observer-disjoint photographs?
3. What can be inferred—and what remains untested—about moisture-dependent colour-pigmentation benefits, costs and genuine genetic morph maintenance?

The active manuscript keeps measurement validation, geographic allocation of ITV, spatial-process decomposition and prospective phenotype-space confirmation as distinct inferential layers.

## Repository scope

This repository is the paper-facing FCP analysis surface. The generalized trait/interaction **spatiotemporal turnover** programme now lives in [`zuizui0223/turnover`](https://github.com/zuizui0223/turnover) and is intentionally excluded from the active FCP manuscript, workflows, results, tests, and development documents here. Earlier pre-split development remains recoverable from Git history.

## Start here

The [ecology-first research and balance audit](docs/POLYMORPHISM_ECOLOGY_WHITE_AXIS_TRADEOFF_RESEARCH_AUDIT_20261014.md) distinguishes (a) observational photo-white moisture associations, (b) prospectively repeated achromatic–chromatic colour geometry, and (c) **unsupported global genetic/pigment fitness equilibrium**. Original classified photo-white is **21,574/45,953 = 46.95%**, but 54,590 photo hues are unclassified; within-species latitude null P=.214 and geographic cancellation null P=.997. Source receipts and arithmetic-only source checks cover both global photo colours and all 439 measured H2 axes under `results/np_unified_20261014/`, `scripts/analysis/reconcile_np_global_to_within_axis_20261014.py`, and `tests/test_np_global_to_within_axis_20261014.py`. The older 2,000 frozen additional taxon identities remain **metadata-only, not high-depth photographed species**.


**Single integrated Full Paper (current journal-facing candidate):** [unified manuscript](docs/POLYMORPHISM_MANUSCRIPT_NEW_PHYTOLOGIST_UNIFIED_20261014.md); [unified Supporting Information](docs/POLYMORPHISM_SUPPORTING_INFORMATION_UNIFIED_20261014.md); [one-paper source ledger](docs/POLYMORPHISM_UNIFIED_EVIDENCE_AND_CLAIM_LEDGER_20261014.md); [Figure 1 water-associated ecology](docs/figures/polymorphism_20261014/figure1_moisture_association.svg).

The numbered list below preserves the **frozen historical high-depth documentation**, not another intended submission.

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
