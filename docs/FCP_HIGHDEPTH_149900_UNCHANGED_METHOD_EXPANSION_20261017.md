# FCP: enlarge the 149,900-photo original study WITHOUT changing the methods

Date: 2026-10-17. Active project: fcp (within-species photographic flower-colour structure, not mois seasonal climate study).

## Only scientific objective

Increase the number of truly newly measured original photographic records and species, holding the actual inferential machinery of the accepted 149,900-photo analysis fixed. No new geographic weighting, null model, 0.25-degree test, new colour classifier, environmental response, or sample selection based on observed colours. This document supersedes the wrong-target PR #205 (closed, not merged).

## Original fixed high-depth study

- Measured input: 149,900 photographs of 1,499 sampled taxa; 500 discovery + 500 validation + 499 prospective species, each originally allotted 100 photos.
- Original candidate frame U0 42,111 species; U100 4,730 species (metadata photo capacity after original observer cap).
- Legacy high-depth = 1,000; historical P500 selected = 500; historical third-cohort selected = 500 (499 ultimately photographed; do NOT reuse the excluded 500th selected taxon).
- Therefore, original outcome-blind 3,230 pre-third candidates minus all 500 previously selected third-cohort taxa = 2,730 species identities still eligible for deterministic outcome-blind selection. This is a photo-opportunity frame, not proven open-flower/flower-colour variation.
- New wave4 first draw: 500 distinct taxon identities, target 100 original photos each, at most +50,000 photo metadata records if every species meets the exact historical gate. Total would then be 199,900 photos after blind ROI classification (not now).
- Freeze selected 500 species from those 2,730 by SHA256 of FCP_HIGHDEPTH_WAVE4_20261017_FIXED_SELECTION_V1|inat_taxon_id|species, ascending hash then taxon ID. The source 3,230 CSV must exactly reproduce SHA 7fc0337074b55a91c0b50773a8d5c3ef82e877074c8b3cacc21f9af58e0ba77e. Original P100/P500 SHA checks and 2026-09-16 500-taxon selection manifest SHA prevent re-use. No phenotype, photo-pixel, climate or geographical span enters species selection.

## Unmodified photo metadata acquisition

The new expansion script invokes the same existing freeze_h9_metadata() function and iNaturalist client from immutable historical commit ca69930986e39e3ea9b2d2f97ab247045ba9db0b. Parameters match the prior 499-species cohort:

- research quality, iNaturalist flowering observation annotation term_id=12, term_value_id=13, species taxon rank, open/unobscured/geo records, GPS uncertainty up to 5 km, same five licence values;
- exactly one randomly ordered API page of 200 per selected species; request interval >=1.05 s, timeout 45 s, zero retries;
- at most 2 photographs per observer, same geographic maximin selection, fixed 100 photo IDs per species;
- exclude all previous ID ledgers (the original 178,462-ID union and 49,999 third-cohort candidate IDs = expected 228,461 distinct original observation/photo IDs);
- retain every failed/missing photo as technical non-availability, never biological monomorphism, without species or photo replacement;
- if request errors >5% or <300 selected species have a complete 100-photo metadata denominator, freeze the original not-evaluable outcome and do not authorize image opening.

The stage 1 workflow uploads a checksum-audited selection manifest, species-level acquisition audit, candidate metadata, the exact 100/species authorized original photo metadata and source receipt. It does not claim images were measured or W/D/50-km statistics recomputed yet.

## Stage 2 (after authenticated stage 1 metadata pass): same original image and statistics machinery

1. Freeze stage 1 successful run ID, authorized metadata SHA and photo-ID set; require no original photo-ID overlap in the real 149,900-photo image study.
2. Recover unchanged detector ROI-v4 + EfficientSAM (original weight/model SHA checks, original nine-colour palette, original four-state mapping, no predictor metadata in blind image worker). Retain 100 frozen rows per authorized species and exactly one terminal classification or technical failure per original photo ID. No photo replacement or new classifier.
3. Re-run existing H1 D calculation and observer-disjoint repeatability with exactly the frozen code. Re-run the original 50-km geographic local depletion test (same 199 within-species complete-vertex permutations, same 40-classifiable/30-local-pair support gates) on the new species, then present results separately and combined with the original cohorts. The output can legitimately change because the sample size changes, but the formulae and gates must not.
4. Re-run frozen H2 W (original achromatic–chromatic axis, 0.10/0.20 vector gates and 999 construction-preserving null worlds) on the additional independently selected species. New batch is a same-source replication of the fixed axis, not a retrofit to the 2026-09-17 prospective test.
5. Archive exact joined photo input SHA and new cohort-level D/local/W statistics. Compare original 149,900 (preserved) vs added 100×N metadata/ROI terminal vs expanded dataset. Report actual classifiable and evaluable numbers, not only nominal downloaded rows.

## Machine boundaries and status

Source: scripts/acquisition/expand_existing_fcp_highdepth_20261017.py.
Workflow: .github/workflows/fcp-highdepth-wave4-original-method-metadata-20261017.yml.
Regression tests: tests/test_expand_existing_fcp_highdepth_20261017.py.

All older scientific claims remain unchanged. New data/biological effects must not be claimed before the post-pixel terminal artifact and inferential replay genuinely exist. Outcomes from mois are not used as FCP species-level estimates.
