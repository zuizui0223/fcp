# Spatial versus temporal allocation of field-counted *Linanthus parryae* colour frequencies

**2026-10-13; post-publication, source-exposed exploratory reanalysis.** No new raw biological observations, no claim of prospective confirmation. This is deliberately a different response from Shannon/Gini–Simpson diversity in PR #172 and not a cross-species estimate.

## Ecological question

Across fixed census sections repeatedly visited over seven calendar years, did the geographical placement of the blue/white floral morph-frequency cline remain stable despite changes in annual abundance? Does side of the ravine predict an unseen year's four section-level white fractions more accurately than a species-wide constant? Does knowing the north/south transect beyond side materially improve transfer?

Original source: Schemske & Bierzychudek (2007), Evolution 61:2528–2543, DOI 10.1111/j.1558-5646.2007.00219.x, Appendix. Published field records are individuals counted as white or blue, pooled over many 2m×2m quadrats into four **large** fixed census sections: north/south transect × east/west of one ravine. There is no section-level true replicate spatial block, no per-plant genealogy, no climate timeseries, and no independent-species replication. Counts were transcribed in PR #172, whose source identity is SHA256 `fa408c7ad6e521a8410252d903b024a057b1ee4f7e27406927e77fcc22bf4ac7`.

## Frozen support gate

Main 4-section × year cell exists if counted blue+white >=20 in **each** of four sections that year, selected without using flower-colour values. Expected evaluable years 1992,1993,1995,1996,1998; 1994 has sparse/empty sections and is HOLD, 1997 has zero plants and is UNDEFINED. The full raw ledger is retained; no missing frequency is replaced by zero. Sensitivity minimum >=100 for all four sections, yielding 1992,1993,1995,1998.

## Main estimands (descriptive)

Let p_{s,t} = white / (blue+white) in fixed section s of census year t (four sections, equal weight; five main years, equal weight), \bar p grand average, \bar p_s section average over years, \bar p_t year average over sections. Define

- V_total = mean_{s,t} (p_{s,t}-\bar p)^2
- V_space = mean_s (\bar p_s-\bar p)^2
- V_time = mean_t (\bar p_t-\bar p)^2
- V_residual = mean_{s,t} (p_{s,t}-\bar p_s-\bar p_t+\bar p)^2.

Check exact identity V_total = V_space + V_time + V_residual. V_space/V_total describes **variance of observed section proportions** on this particular support; it is NOT an ecological variance explained by environment, a census-state Shannon fraction, F_ST, or a confidence interval. One individual is not a replicate for the cross-year effect.

## Leave-one-year-out (LOYO) geographic transfer

For each main year t, train on the other complete years only (no data from the held-out year for any model):

- GLOBAL: equal-section-year grand white fraction over the training years.
- SIDE: per-side equal section-year mean over training years, ignoring north/south.
- SECTION: per-section mean over training years, including north/south.

Score each test section's squared deviation of observed **frequency** from the prediction. First average equally over four sections, then equally over held-out years. This comparison is a descriptive difference of fraction-MSE (not plant-count Brier likelihood) with no pseudo-independent binomial p-values. Report full held-out ledger and each year's gain. Purely descriptive relative MSE reduction = (GLOBAL-SIDE)/GLOBAL; global >0 required.

## Additional audit

Show density ratio max/min across eligible years; both sides' equal-transect white frequencies per year; annual east-minus-west contrast; all four fixed sections' blue/white nonzero supports; analogous sensitivity minN=100; 1994 and 1997 unaltered. Do not reclassify 1997 as monochromatic or a demographic extinction. Preserve the source's 1992 south/west percentage discrepancy (printed 1.2% vs 26/1383=1.87997%): perform calculations from integer counts only.

## Scientific ceiling

Even a large stable contrast does not identify whether geographic distribution was shaped by localized selection, dispersal limitation, founder events or planting history, and it does not rule out colour-frequency variation at the unsupplied 2m² scale. A transect side is a real landscape split but not an independent biogeographical realm. The source's reciprocal transplant and neutral allozymes provide independent contextual evidence **from that publication**, not newly conducted validation. Do not upgrade this study into a universal temporal law or a mechanistic climate/hypervolume test.
