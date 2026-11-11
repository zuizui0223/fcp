# FCP — Independent NHM botanical flower-colour descriptions do not verify the east rainfall white–blue/purple signal

**Date:** 2026-10-13  /  **Status:** verified original-source GitHub Actions **38405377887 SUCCESS**; this is observational **independently derived colour phenotype**, not independent source climate or an out-of-sample genetic/fitness study. Corresponding FCP Draft PR #191 is outside main. No FCP New Phytologist H1/H2 or historical PR #164/#188–190 effects are changed.

## Central scientific decision

The fixed FCP photographic result (east 90–180° white-minus-blue/purple rainfall contrast) was previously positive in full-source, species-equal 4-colour images under a selected family of 24 contrasts but weakened under genus-equal weighting or high image-quality restriction. **The independent botanical-description measurement did not give statistically conclusive positive evidence for that east rainfall hue direction under any predefined genus weighting / four-sector family test.** A positive point estimate alone is **not confirmed support**. This does not establish an opposite sign, biological null effect or refute a different four-hue photographic estimand.

## Exact independently sourced phenotype and FCP environmental opportunity

- **Independent colour**: Araujo, Hartley & Brummit (2026), Natural History Museum dataset, DOI `10.5519/p3dm31kc`, CC BY 4.0. GitHub raw `NaturalHistoryMuseum/plant-traits` commit `4843c22bb2fc91018d9828dcddebf067c83b1097`, file `output/raw/ocr-traits/combined.csv`, Git blob SHA-1 `e1cfd2d5c518c2061ec91d1a9e8b64fc7ca69389`, byte length 17,279,412. **59,369 original named plant records**; 6,005 strict white, 3,580 strict blue/purple/violet, 10,579 other-only colour, **39,205 HOLD/unknown/mixed/conflicting descriptions**. No unclassified botanical name is imputed.
- **Original FCP geographical opportunity, not its measured photo hue**: 100,543 original photo IDs from 42,111 named species, SHA256 `2d6097ba1b0434632e13bbe79dd78180608a287ff5dabebc770b5c0a2cfb52a0`; **100,255** climate-and-geolocation complete original site records from 42,055 species. This includes records classified as *unknown* by the original FCP image palette.
- Exact matching unambiguous scientific binomials yields **2,476 botanical colour-labeled species in 1,144 genera** with **8,120 FCP site records**; within those: **1,538 white species, 938 blue/purple species**.
- Fixed main east-vs-rest group requires >=3 botanically measured species on BOTH sides of a common genus: **58 matched genera**, 670 distinct taxa / **2,355 FCP photo-site environmental opportunities**; 323 source species photographed East and 473 elsewhere (overlap across geography, so do not sum to 670). East recorded 186 independently described white and 137 blue/purple species; outside 265 white and 208 blue/purple.

These source labels are not confirmed individual genets or within-species morph frequencies. The 2,476 exact FCP matched species are neither a random sample of all 42,111 photo opportunity species nor a new independent species occurrence source. Mixed-colour species (biologically important) are conservatively excluded from the strict binary contrast.

## Pre-fixed main results (min 3 named species per side per matched genus)

4 sectors: `W180–W90`, `W90–E0`, `E0–E90`, **`E90–E180`** (selected east sector already known from FCP parent). Unit: **percentage-point contrast per one SD of BIO12 within all FCP photo-opportunity climates (approximately 763 mm)**. This is the slope difference between botanical WHITE and BLUE/PURPLE probabilities in the binary pure-colour subset. It is **not identical in denominator or sampled species to earlier photographic four-hue PR #189 contrasts**.

| Sector versus other sectors | Species-equal contrast, pp | Species-equal bootstrap95 | Genus-equal contrast, pp | Genus-equal bootstrap95 |
| --- | ---: | --- | ---: | --- |
| W180–W90 | -16.19 | [-44.85,+8.30] | -18.29 | [-52.24,+10.61] |
| W90–E0 | -27.53 | [-49.98,+12.21] | -16.18 | [-41.22,+14.08] |
| E0–E90 | +10.46 | [-13.04,+33.70] | +0.58 | [-16.27,+16.85] |
| **E90–E180** | **+2.07** | **[-13.87,+19.07]** | **+3.56** | **[-15.34,+22.60]** |

The 999-replicate paired-genus bootstrap and 4-sector simultaneous maxT intervals retain zero in **all four regions**. Primary east 4-region maxT ~95% bounds:
- species-equal: **[-19.58,+23.72] pp per SD**;
- genus-equal: **[-20.59,+27.70] pp per SD**.

This is **NO CONVERGENT INDEPENDENT COLOUR-SOURCE SUPPORT** at the fixed inferential gate. Since intervals are wide, this is an **inconclusive cross-measurement test**, not evidence of physiological absence or exact equality.

## Frozen support sensitivities (none upgrade the result)

- Primary >=3 different taxa per geographic side per genus, **>=2 original FCP photograph opportunities per focal taxon**: 40 matched East genera /352 taxa, east species-equal +3.39 pp CI [-13.81,+24.77], genus-equal +12.20 pp CI [-5.92,+32.92]. Global 4-sector maxT HOLD because the W180–W90 cohort had only 24 qualifying genera (pre-fixed minimum 30).
- Relaxed >=2 taxa per side per genus (exploratory): East 112 genera /931 taxa, species-equal +8.85 pp CI [-5.01,+25.20], genus-equal +13.03 pp CI [-2.38,+30.09]. 4-sector maxT also spans 0 for both.
- Relaxed >=2 per side +>=2 FCP site observations: East 91 genera /561 taxa, species-equal +6.27 pp CI [-11.91,+28.09], genus-equal +13.49 pp CI [-9.91,+36.59]. 4-sector maxT both span 0.

None of these nested taxon opportunity subsets is an independent biological confirmation. No subset was promoted in place of the frozen main gate.

## Why this changes the scientific argument

The **original photo phenotype** yields a positive, taxonomically weighting-sensitive regional photo-colour association. The independently text-described **botanical phenotype** does *not* currently show a comparable statistically identified regional direction. The following explanations cannot be separated here:

- true association confined to certain species, taxa, and/or habitats not covered by descriptive trait extraction;
- differential image exposure/flower-organ interpretation vs conservative prose colour coding;
- within-species local colour differences excluded by species-wide botanical descriptions;
- the different denominator (four-colour photo probabilities versus strict two-colour textual species probabilities);
- spatial, lineage and composition confounding; geographic habitat occupancy and original citizen-science sampling biases;
- low trait overlap and imprecise genus-block uncertainty.

Critically, the evidence **does not currently justify** saying annual precipitation causes global white-flower adaptation or that blue/purple pigment chemistry has been suppressed by humid climate. Nor does the new failure prove no plant evolutionary environmental filtering. A stronger phenotype-independent conclusion is that the fixed photo-based east white–blue geographical contrast is **not cross-measurement robustly established**.

## Next genuinely independent science

A prospective test needs both **botanically/physically verified organ colour** and **independent species occurrence or sampled population climate data**, ideally across multiple lineages with adequately sampled *white and blue/purple* states. Sampling should pre-fix region, hue definitions, allowed multi-colour species handling, taxonomic matching, photo-quality/missingness rules, the number of distinct genera, and the 4-sector multiplicity family. New colours for the same NHM taxa would not independently replicate the same source.

## Independent original-data verifier — SUCCESS

[Final GitHub Actions #38405377887](https://github.com/zuizui0223/fcp/actions/runs/38405377887) re-downloaded the exact NHM original (17,279,412 bytes and pinned Git blob SHA-1), the FCP site-source archive (pinned SHA-256), and compared the full frozen `result.json` with the source replay. A **second script independently read both sources and recalculated 30 non-HOLD region×support×weight coefficients by direct species-row weighted least-squares**; 2 additional conditions were correctly marked predeclared HOLD. Maximum coefficient difference **8.33e-16 probability units**; synthetic tests 6/6 passed. The bootstrap 95% CI draws were **not independently resimulated** by this verification script; source-replay and frozen numerical output were compared by the GitHub workflow.

Artifact **11818501773**, GitHub Actions archive SHA256 **12e1889f7e3d8518e19f7d5228de4db68e1e877e91b5c1a65e60513b2119befc**. No new independent biological sample was generated.

Source: [NHM data portal](https://data.nhm.ac.uk/dataset/angiosperm-functional-traits-extracted-from-botanical-descriptions); [original GitHub data](https://github.com/NaturalHistoryMuseum/plant-traits); [FCP original geographic hue PR #189](https://github.com/zuizui0223/fcp/pull/189), [PR #190 genus leverage](https://github.com/zuizui0223/fcp/pull/190). All statistics are conditional observational estimates and source-verified in the new Actions result JSON.
