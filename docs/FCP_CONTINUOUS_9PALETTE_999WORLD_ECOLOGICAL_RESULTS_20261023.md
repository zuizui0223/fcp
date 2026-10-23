# FCP: authenticated nine-palette continuous geographic structure (23 October 2026)

**Status: SOURCE-VERIFIED, EXPLORATORY, NOT A GENETIC/MORPH ADAPTATION RESULT.** The 999-world run and all existing repository guards succeeded: [GitHub Actions 39921400989](https://github.com/zuizui0223/fcp/actions/runs/39921400989). Exact machine-read source/results receipt: [999-world JSON](../results/fcp_all4475_continuous_9palette_geography_999worlds_20261023.json). Full 999-world summary and two taxon-only CSVs reside in this Actions run's artifact 12486224900; raw coordinates/photo/observer IDs are deliberately not published.

## What was actually tested

A fixed 4,475-taxon register contains 447,500 original photos; 227,107 pass the frozen four-state source photo classification and have usable nine-hue flower ROI composition. Each photo's original nine-dimensional palette fractions were compared by base-2 Jensen–Shannon divergence. A species-equal statistic measured how much nearer photographs within 5, 25, 50, 100 and 250 km were in hue-composition space than photos across that species' sampled geographic range.

Three 999-world species-conditional nulls permuted complete nine-hue vectors while retaining (1) species sampling positions, (2) additionally month of photo observation, or (3) additionally the photo's frozen coarse four-colour label. Original fixed four-class 50-km P2 was reproduced *exactly* before these new post hoc analyses: all-eligible 3,600 taxa +0.016115883798; strict 1,761 taxa +0.0176149021685. Analysis scripts and six synthetic falsification tests are part of PR #261.

## New source-supported ecological results

Under the strict >=40 classifiable images and >=30 local pair definition, 1,761 species (832 genera) permit a 50-km continuous analysis:

| Nine-palette JSD contrast | Source-verified 50-km excess | Genus bootstrap 95% CI |
|---|---:|---:|
| Month-conditioned distribution | **+0.0114075** | [+0.0095764,+0.0133275] |
| Month AND frozen 4-state label conditioned | **+0.0051116** | [+0.0042504,+0.0059493] |
| One photo per observer, month conditioned (1,655 species) | **+0.0112671** | [+0.0092318,+0.0131267] |
| One photo per observer, month + label conditioned | **+0.0053571** | [+0.0044341,+0.0062717] |

Each corresponding one-sided Monte-Carlo p is **0.001** (minimum attainable for 999 worlds), five-radius Holm is **0.005**, and Holm over all **15 radius × null-family** comparisons is **0.015**. Every tested radius from 5 to 250 km had a positive mean residual under month+label conditioning and the same 15-test Holm adjustment 0.015. Positive aggregate effects **do not** mean that every species is structured or that a photographed flower has an inherited discrete colour morph.

With the species set fixed at the 372 taxa having >=30 even 5-km photo pairs, the month+coarse corrected mean excess declines from approximately +0.0054 at 5 km to +0.0018 at 250 km: *paired same-species* difference +0.003644, genus 95% [+0.001587,+0.005546]. The short-radius caveat is especially severe because original allowed reported coordinate accuracy up to 5 km. The earlier, genuinely stricter same-year/month <=5-km distinct-observer null **remains unsupported**, and must not be overwritten by this different new test.

## Interpretation

The strongest evidence-bounded result is **hierarchical spatial organisation of photographed floral colour**: geographically concentrated variation is present both across the original coarse colour categories and within those category labels, in measured continuous nine-palette composition. Restricting the null to category × calendar month retains spatial evidence, so the broad-scale signal is not *solely* a consequence of those coarse bins or a simple calendar-month composition shift.

This is interesting for floral intraspecific phenotypic variation because assigning every species "2/3/4 morphs" is not necessary to find its geographic structure. It **does not prove** a second genetic hue gradient, local adaptation, a clinal genotype frequency change, an extra independent morph or pigment conservation cost/benefit. The 9-channel feature is a photographed ROI palette mixture, not calibrated petal reflectance. Residual locality-correlated lighting, background contamination, petal stage, observer/camera and cultivation status remain realistic alternatives. Nearly half of the original photo set failed the coarse classification gate; classification missingness could bias species and locality results.

## Next evidence discriminator

1. Blind independent original-petal organ/stage review, sampling across hue scores, geography and classification failures.
2. Flower-minus-background spatial negative control and within-year/month cross-observer matching for the continuous nine-palette residual.
3. On independently validated continuous vectors, compare unimodal gradients against stable two- to four-mode *species-specific* distributions via observer-disjoint cross-validation; do not force K=2 or equate photograph bin count with inherited morph number.
4. Only for confirmed natural multi-phenotype taxa, relate geography to environmental gradients, gene flow and molecular/fitness mechanisms.

This is an addition to the active FCP intraspecific paper's evidence; it does not modify the original published-equivalent historic hypotheses, and the current New Phytologist working manuscript remains unsubmitted.
