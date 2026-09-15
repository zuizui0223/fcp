# Measurement references for the H1–H3 manuscript

Verified 2026-09-15. Literature-only note: no FCP photographs were opened, no analyses rerun, and no validation outcome or admission rule changed. These references provide context, not external certification of FCP measurements.

## Laitly et al. (2021)

Laitly, A., Callaghan, C. T., Delhey, K., & Cornwell, W. K. (2021). Is color data from citizen science photographs reliable for biodiversity research? *Ecology and Evolution, 11*(9), 4071–4083. [DOI: 10.1002/ece3.7307](https://doi.org/10.1002/ece3.7307).

**Access:** Primary article full-text methods and results were available through [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC8093748/); a subsequent request encountered a browser challenge. This is not an abstract-only assessment.

**Evidence:** The study compared 9,441 photographs of 537 Australian bird species with museum spectrometry, and compared measurement methods for two plant species. Bird species means agreed better with spectral measurements than individual photographs. Its approximate 12–14-photo plateau concerns interspecific mean estimates, not a universal sample size for within-species diversity. Plant measurements used manually filtered photographs and sampled petal points. Reduced photographic control increased variation; larger samples only partly mitigated error, including saturation problems. Ordinary photographs did not recover UV information. Relevant sections: 2.1, 2.2, and 3.1. [Full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC8093748/).

**FCP interpretation:** Supports testing aggregation and repeatability, not assuming accuracy. It does not validate FCP's automatic ROI, diversity D, or white-axis statistic.

## Luong et al. (2023)

Luong, Y., Gasca-Herrera, A., Misiewicz, T. M., & Carter, B. E. (2023). A pipeline for the rapid collection of color data from photographs. *Applications in Plant Sciences, 11*(5), e11546. [DOI: 10.1002/aps3.11546](https://doi.org/10.1002/aps3.11546).

**Access:** [Publisher full text](https://bsapubs.onlinelibrary.wiley.com/doi/10.1002/aps3.11546) inspected, including Methods, image-correction Results, and Discussion; PMC initially returned a browser challenge.

**Evidence:** Analysed 4,886 cleaned Erysimum occurrences using manually selected pixels, avoiding shaded, overexposed, or withered petals. The principal colour range was yellow–orange–red, represented by HSV hue. A separate 439-flower ColorChecker collection supported comparisons with corrected images; regional patterns were similar, while grayscale panels tended to have greater error. The authors explicitly discuss background similarity and petal shading as problems for automated extraction. Relevant sections: Data acquisition; Impacts of image color correction; Discussion. [Full text](https://bsapubs.onlinelibrary.wiley.com/doi/10.1002/aps3.11546).

**FCP interpretation:** Supports feasibility in a restricted, manually screened system. It does not establish accuracy for automatic multi-species ROIs, white/nonwhite contrasts, or pollinator perception. Its calibration comparison cannot be transferred as FCP's own validation.

## Troscianko and Stevens (2015)

Troscianko, J., & Stevens, M. (2015). Image calibration and analysis toolbox – a free software suite for objectively measuring reflectance, colour and pattern. *Methods in Ecology and Evolution, 6*(11), 1320–1331. [DOI: 10.1111/2041-210X.12439](https://doi.org/10.1111/2041-210X.12439).

**Access:** [Publisher full text](https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.12439) inspected, particularly Equipment checklist and Taking and processing photos.

**Evidence:** The workflow extracts linear measurements from RAW files and normalizes using grey standards. Visual-system-specific conversion additionally requires camera spectral sensitivities; UV-sensitive systems require suitable UV imaging. These are acquisition/calibration requirements, not properties automatically conferred by converting ordinary RGB values to another colour space. [Full text](https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.12439).

**FCP interpretation:** Use image-derived colour descriptors, not calibrated reflectance or pollinator-perceived colour, unless those additional measurements exist. This methodological distinction does not imply that every uncalibrated-image association is invalid.

## Manuscript claim boundary (synthesis, not an external validation result)

H1 repeatability and tissue/colour accuracy are different questions. H2's white/nonwhite structure needs target-domain checks for exposure and ROI contamination. H3 associations or non-support should remain conditional on the image-derived response and sampled observations; these papers supply no causal evidence for FCP's ecological mechanisms. Neither increasing sample size nor converting colour space alone establishes those missing links. Do not use any of these papers to relax frozen thresholds, reinterpret a failed measurement check as a pass, or authorize P500 opening.
