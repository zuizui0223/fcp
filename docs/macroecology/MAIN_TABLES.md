# Main tables — macroecology companion manuscript v0.1

## Table 1. Source/analysis populations and independence of units

| Evidence unit | Data source | Original numerator / usable N | Unit of inference | Independence and caveat |
|---|---|---:|---|---|
| Outcome-blind discovery opportunity | FCP iNaturalist original record union | 100,543 photographs / 42,111 named taxa | Photographic opportunity | Not a global species prevalence denominator or random biological sample |
| Classified photograph colour | FCP four-hue photo system | 45,953 photographs; 54,590 unclassified | Sampled visible photographic hue | Same system across photographic analyses; missing class not biological white |
| Climate26-complete photographic records | FCP+WorldClim/other climate | 45,836 photographs / 25,028 taxa | Photo-site environment + sampled hue | Source study pool, not eligible genus-regression size |
| Soil36-complete photographic records | FCP+climate+SoilGrids | 31,925 photographs / 20,786 taxa | Photo-site environment + sampled hue | Soil completeness alters the species set |
| Genus-conditioned climate direction (global, all photo quality) | FCP macro analysis #193 | 20,694 eligible named species / 1,934 genera | Among-species genus-centred slope | All species in qualified genera; broad coordinates controlled |
| Joint climate26 predictive test | FCP joint abiotic prediction #166 | 17,183 held-out species | Genus-equal out-of-geocell Brier | Same original images, new held-out folds; not a new biological source |
| Joint soil36 predictive test | FCP joint abiotic prediction #166 | 13,961 held-out species | Genus-equal out-of-geocell Brier | Different species and environmental completeness |
| Matched global B/C spatial comparison | FCP #194 | 9,370 species / 845 genera / 14,268 photograph records | Same species under two spatial adjustments | Fixes species/colour but changes spatial conditioning |
| Matched eastern B/C spatial comparison | FCP #194 | 2,886 species / 280 genera / 4,149 photo records | Same species under two spatial adjustments | East source 90–180 E chosen after prior outcome inspection |
| Botanical description diagnostic | NHM botanical texts + FCP original site climate | 2,476 exact matched white/blue species / 1,144 genera / 8,120 photo sites | **Different species-level binary colour endpoint** | Independent textual phenotype, **same iNaturalist geographic sampling source** |
| NHM eastern paired-genus cohort | NHM botanical texts + FCP photo sites | 670 species / 58 genera / 2,355 photo sites | East versus elsewhere two-colour response | Inference does not transfer exactly to four-hue photographic proportions |

## Table 2. Annual precipitation and white photographic probability: effect of spatial and quality conditions

Effect: **percentage-point change in photographed white-category probability per global source population SD in annual precipitation (BIO12)**, genus-equal unless otherwise stated. A simultaneous interval is corrected **only for that row's within-cohort environmental-anchor × 4-colour family**, not all historically explored regions/models. Paired B–C rows are paired estimates on an exactly identical species set.

| Population / comparison | Species | Genera | BIO12 → white (pp/SD) | 95% uncertainty | CI definition |
|---|---:|---:|---:|---|---|
| World, all classifiable, broad genus/geography | 20,694 | 1,934 | **+2.255** | [+0.319,+4.190] | 24-hue×anchor within-scope maxT |
| World, strict source photographic quality | 5,848 | 727 | +1.229 | [−2.046,+4.504] | 24-hue×anchor within-scope maxT |
| World, genus×coarse-cell conditioning (not same sample as above) | 9,370 | 845 | +0.808 | [−3.556,+5.172] | 24-hue×anchor within-scope maxT |
| East 90–180 E, all classifiable, broad geography | 5,774 | 689 | +4.716 | [+1.565,+7.866] | 24-hue×anchor within-scope maxT; eastern sector selected posthoc |
| East 90–180 E, strict source photographic quality | 1,459 | 216 | +2.300 | [−3.700,+8.301] | 24-hue×anchor within-scope maxT |
| East 90–180 E, genus×coarse-cell conditioning | 2,886 | 280 | −2.076 | [−8.794,+4.643] | 24-hue×anchor within-scope maxT |
| World, **same 9,370 species**, broad geography (B) | 9,370 | 845 | +1.326 | [−0.709,+3.575] | Stage-specific nominal 95% genus-bootstrap |
| World, **same 9,370 species**, genus×coarse-cell (C) | 9,370 | 845 | +0.808 | [−1.884,+3.754] | Stage-specific nominal 95% genus-bootstrap |
| **World paired B−C** | 9,370 | 845 | +0.518 | [−1.253,+2.289] | 24-hue×anchor simultaneous maxT on paired contrasts |
| East, **same 2,886 species**, broad geography (B) | 2,886 | 280 | +1.102 | [−2.601,+4.749] | Stage-specific nominal 95% genus-bootstrap |
| East, **same 2,886 species**, genus×coarse-cell (C) | 2,886 | 280 | −2.076 | [−6.463,+2.132] | Stage-specific nominal 95% genus-bootstrap |
| **East paired B−C** | 2,886 | 280 | +3.178 | [−0.258,+6.613] | 24-hue×anchor simultaneous maxT on paired contrasts |

**Do not append the NHM binary +3.56 pp estimate as another row of Table 2**: it has a different species phenotype, probability denominator, environmental source-SD and sampled taxon composition. It belongs to Table S3 / Figure 4 with its own units and population.

## Table S3 (working) — independent botanical-description diagnostic

| Region/comparison | Genus sample | Effect | 95% CI | Scientific status |
|---|---:|---:|---|---|
| East 90–180 E vs others, species-equal binary white−blue | 58 genera, 670 species | +2.07 pp per all-FCP-site BIO12 SD (~763mm) | [−13.87,+19.07] | Source-conditional, not conclusive |
| East 90–180 E vs others, genus-equal binary white−blue | Same | +3.56 pp per all-FCP-site BIO12 SD (~763mm) | [−15.34,+22.60] | Source-conditional, not conclusive |
| East species-equal after four-sector maxT | Same | +2.07 | [−19.58,+23.72] | Does not exclude zero |
| East genus-equal after four-sector maxT | Same | +3.56 | [−20.59,+27.70] | Does not exclude zero |

NHM phenotype sources classify species, not individual flowers. The climate locations remain the original FCP photographic observation sites, not independent GBIF preserved-specimen field localities.
