# FCP Validation: frozen metadata-to-image-quality review gate — 2026-10-09

## Status and scope

The original source-only 10-km species-cohort feasibility gate is **HOLD** and remains unchanged. This note concerns a *retrospective supplementary quality-review opportunity*, not the frozen New Phytologist H1/H2 outcome, a genetic flower-colour polymorphism prevalence or an independent 2,000+730-species future test.

Source-verified GitHub Actions:

- [Gap-one first-page metadata, run 37798611129](https://github.com/zuizui0223/fcp/actions/runs/37798611129): **15/22** Validation species with a potentially usable unused observation/photo ID under fixed site-month-year, observer, licence, positional-accuracy and prior-photo-exclusion rules. Archive 11558978651, ZIP SHA256 `8ac04d7408b714d7b1cd9945414ffed8c322d3ae035e93881952c92934536db7`.
- [Gap-two independent-observer year cells, run 37801452055](https://github.com/zuizui0223/fcp/actions/runs/37801452055): **10/60** Validation species meet **both** required missing observer-photo slots; 50 first-page-insufficient, zero unresolved query errors/pagination. Exactly 120 one-page taxon-year requests, synthetic tests passed, full API responses checkpointed before aggregation. Archive 11561291901, ZIP SHA256 `8ca425141ff669d81702b94e4f30f5bb0c1ca5cb4e05cee8fd56ea2314e7c3ea`.
- Permanent checked summary receipts: `results/fcp_validation_gap1_live_metadata_20261008/terminal_pilot_receipt.json` and `results/fcp_validation_gap2_live_metadata_20261009/terminal_pilot_receipt.json`.

The two immutable GitHub Actions ZIPs were decoded **without opening any photograph pixels** and reconciled against the frozen `source_only_observer_photo_query_queue.csv` from those artifacts. The retrospective **25 distinct candidate species** have **115 distinct public photo metadata IDs** (39 attached to gap-one taxa, 76 attached to gap-two taxa); the full response-ID records remain in those source archives and a locally generated quality-review JSON ledger with SHA256 `5e3a13c2439354cdb97cc589861cbac13b6c0c50e39fa5e932d996a54ea5ad41`. None of the 115 source images has yet been reviewed for actual flower visibility, organ/identity, exposure or colour measurability.

## Gate: what a photograph must establish

1. Preserve the **exact observation/photo ID** in the archived metadata; never add a colour-selected replacement or treat absent photos as monomorphic.
2. Validate the species and flower organ directly from the source image, rather than relying on the observation's iNaturalist taxon assignment alone.
3. Audit focus, field of view, camera exposure/clipping and usable floral colour, recording missing or ambiguous as such. White visible-photo labels have known exposure dependence.
4. Retain the exact *original* target year, calendar month, site-anchor photo and observer-disjoint requirements. For gap-two taxa, **every** deficient year cell must have sufficient independent, successfully inspected photographic evidence.
5. Assign one of `ACCEPT_QUALIFIED_FLOWER_IMAGE`, `REJECT_NOT_FLOWER_OR_WRONG_ORGAN`, `REJECT_SPECIES_MISMATCH`, `REJECT_PHOTOMETRIC_UNUSABLE`, `REJECT_DATE_SITE_OBSERVER_OR_LICENSE`, `UNREVIEWED`; do not silently drop failures or claim genotypes/fitness.
6. Only after such QA can the frozen, separate photo-colour measurement pipeline be run on accepted photos. That is a *new retrospective supplement*, not untouched species-disjoint prospective confirmation.

## Sharp feasibility threshold

Before adding any new photo pixels: 10 original species meet the conservative 10-km source-photo requirement. If all 15 gap-one and all 10 gap-two metadata-candidate species later meet their full image-quality criteria, the Validation metadata-only **upper bound** is **35**. The required pre-defined sample-size gate is **30**, so >=20/25 candidate species (80%) must become fully qualified. There can be **at most 5 failed candidate species**; one unusable gap-two year cell makes that entire gap-two species unqualified. Without new image measurement, the realised validated sample remains 10, not 35.

**No year-specific climate, biochemical pigment, neutral demographic process, gene flow or selective fitness effect is identified.** Do not relax the 10-km locality definition or contaminate the 2,000+730 prospective reserve to produce a positive result.
