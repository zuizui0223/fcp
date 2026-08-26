# JBI FCP mechanism hypotheses — v1

## Purpose

This registry separates **hypothesis generation from the FCP literature itself** from **independent comparative tests**. Mechanism text found in the same papers used to document C/S organization is never treated as independent confirmation of C/S causation. Keyword hits in the candidate queue are navigation aids only.

## Current synthesis

The working ecological model is no longer “larger climatic niche breadth causes one spatial state.” The sharper contrast is:

- **C (local coexistence):** mechanisms that keep alternative morphs in the same population by temporal fluctuation, opposing agents, frequency dependence, reproductive assurance, or persistence/storage.
- **S (spatial segregation):** mechanisms that make the direction of selection or demographic history differ among places, including pollinator geographic mosaics, edaphic/local adaptation, restricted gene flow, drift/founder history, or reproductive character displacement.
- **C+S:** both layers can operate at different spatial/temporal scales; it is not an intermediate or coding nuisance.

This is a hypothesis architecture, not a claim that C and S each have a single cause.

## Evidence boundary already established

1. Five static climatic breadth/heterogeneity metrics do not robustly explain the v6 C/S states.
2. Mean BIO4 (temperature seasonality) is a weak C-side candidate, but it loses robustness after mean-temperature/latitude controls and is **not** interannual temporal variability.
3. Range size, GBIF fragmentation and static spatial climate-turnover proxies do not robustly explain S.
4. Self-compatibility vs self-incompatibility is not explanatory. Autonomous reproductive assurance remains a weak C-side candidate and must be tested separately from compatibility.
5. CLO-PLA clonality coverage is only 4/18 pure-state species, so the clonality hypothesis is **not tested**, rather than negative. A separate global growth-form lookup also has state-imbalanced missing woodiness and is not used inferentially.
6. The same-literature screen finds direct spatial local-adaptation evidence in 3/7 S-only and 0/11 C-only species (two-sided Fisher p=0.0429). This is useful for prioritizing an independent soil/pollinator test, but it cannot be treated as confirmatory because the evidence source is coupled to the literature defining the spatial state.

## Literature anchors

- Sapir Y, Gallagher MK, Senden E. 2021. *What Maintains Flower Colour Variation within Populations?* Trends in Ecology & Evolution 36:507–519. doi:10.1016/j.tree.2021.01.011. The review explicitly covers multiple pollinators, opposing/fluctuating selection, heterozygote advantage, frequency dependence and neutral alternatives.
- Trunschke J et al. 2021. *Flower Color Evolution and the Evidence of Pollinator-Mediated Selection.* Frontiers in Plant Science 12:617851. doi:10.3389/fpls.2021.617851. It emphasizes that pollinator-mediated selection is context dependent and direct evidence is scarcer than simple preference observations.
- Sapir et al. 2017 Mediterranean FCP review (PMID 28430395). It highlights pollinator and non-pollinator selection, gene flow, drift, autonomous selfing and clonal reproduction, and explicitly treats within/among-population colour organization as a continuum of geographic outcomes.
- v6 focal examples used only for mechanism mapping include Disa porrecta (doi:10.1093/aob/mcaf074), Eruca sativa (doi:10.1093/jisesa/iez038), Leptosiphon parviflorus (doi:10.1002/ajb2.70018), Gentiana leucomelaena (doi:10.1002/ece3.2899), Silene littorea (doi:10.1111/plb.13209), and the Raphanus sativus pollination/herbivory series.

## Ordered tests

The next independent tests are ordered by biological specificity and attainable coverage:

1. **Temporal environmental heterogeneity → C:** frozen dynamic TerraClimate test at fixed occupied locations. Protocol is complete; real remote extraction is still pending and no result is claimed.
2. **Edaphic/local adaptation → S:** attach soil chemistry/texture heterogeneity to the same fixed occurrence geometry, then ask whether geographic environmental partitioning predicts S rather than total niche breadth.
3. **Pollinator geographic turnover → S vs multi-guild temporal turnover → C:** use independent interaction/occurrence data rather than counting pollinator words in FCP papers.
4. **Autonomous reproductive assurance → C:** expand quantitative autofertility coverage; do not substitute SC/SI.
5. **Gene-flow/connectivity → S:** population-genetic or dispersal-informed connectivity is required; GBIF fragmentation is retained only as a failed proxy.
6. **Persistence/storage mechanisms → C:** seed-bank persistence and clonality once a global dataset with balanced coverage is available.

## Reproducibility

`scripts/literature/build_v6_mechanism_candidate_queue.py` generates the 110-source navigation queue by joining the v6 core membership sources to the 12,064-record screen. Its regex flags are deliberately named `*_candidate`; manual verification or an independent external dataset is required before inference.
