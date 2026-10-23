# FCP: raw-photo organ and phenology falsification of apparent colour modes

**2026-10-23 / Source-locked visual pilot, NOT blinded human botanical validation.**

## Where the study stood

The original-source 447,500-photo FCP analysis (4,475 species, 227,107 eligible four-colour photo labels) found species-conditional local nine-hue similarity at 5–250 km. The later independently fitted phenotype geometry (Hellinger nine-hue PC1, BIC Gaussian mixtures against unimodal normal/skew-normal competitors) identified **27 recurrent photographic density-mode candidates** among the 1,761 strict photo/spatially eligible taxa, 43 unstable and 1,691 not resolved. All 27 had two detected *density peaks* in the tested 1D projection, even though four fitted K=3 components. The month-matched P2 outcome had an exploratory mode-candidate versus undetected association +0.032633 after continuous diversity, photo count and near-pair adjustment; 499-genus bootstrap interval [−0.000376,+0.066839] overlapped zero.

27 candidate and 27 matched control taxon-photo packs were then frozen in PR #263. Of 594 selected licensed, authentic iNaturalist photo IDs, 589 original image thumbnails were successfully recovered; **5 could not be recovered and remain missing**. Two reviewer packages are anonymous with respect to source colour, case/control and geographical location; the PI mapping must be kept separate and is not a confidential GitHub vault. **No independent human petal annotations have been completed.**

### Important new distinction

The photo-PC1 peak does not say what makes the colour differ. Author-side inspection of the *actual authenticated source images*, only after the review-packet source mapping was opened, suggests at least four different sources of multimodality:

1. different open **corolla petal colours**;
2. different *colour pattern on petals*, such as maroon blotches on white corollas;
3. **bud/anthesis/fruit-stage** colour changes within a taxon;
4. differences among **petals vs sepals/fruit/inflorescence or background**, sometimes with white petals at all photographed extremes.

The author-side image observations below are *AI-assisted qualitative pilot annotations*, not separate human reviewer measurements. No frequencies of these four sources are estimated across the 27 by this incomplete qualitative set.

## Visual pilot cases (actual recovered original source photos)

| Taxon | Image-supported observation from 3 low, 3 central and 3 high PC1 samples | Preliminary decision |
|---|---|---|
| Viola lutea | Open yellow corollas cluster among low-PC1 images and open violet/purple corollas among high-PC1 images; bicolour flowers appear in central samples | **Priority**: visible petal-hue contrast, no genetic inference |
| Leptosiphon parviflorus | Open yellow flowers are seen among low-PC1 images, mostly white among high and white/lavender in central; one low-PC1 photo includes unopened structures | **Priority**: photo petal hue contrast with stage caveat |
| Delphinium hesperium | Pale pinkish and violet flowers vs deep violet/purple flowers are visible; other yellow flowers appear in the scene in some low-PC1 images | **Priority**: verify the target petal/flower ROI and source identification |
| Xyris caroliniana | White and yellow small open flowers occur in different photo-PC1 groups | **Priority for taxonomic check**: visible open-flower contrast, botanical identity not independently verified |
| Cistus ladanifer | Corollas are broadly white; maroon basal petal blotches differ in presence/prominence across source photos | **Pattern phenotype**, not equivalent to a yellow-violet corolla hue polymorphism |
| Argemone albiflora | Open corolla petals appear white across inspected low, central and high images; yellow reproductive centres/background change | **Red flag**: photo modes do not visibly demonstrate distinct petal hue |
| Ochna serrulata | Low-PC1 images show bright red structures around black fruit; high-PC1 images show open yellow flowers | **Strong organ/stage red flag**: different botanical organs, not red-versus-yellow open corolla morphs |
| Oxalis tenuifolia | High-PC1 images show elongated red-pink mostly unopened flower buds; middle/low include pale/white open corollas | **Strong floral stage red flag**: buds versus open flowers |
| Anemone caroliniana | Most inspected flowers appear white despite photo-PC1 grouping, with an occasional bluish flower | **Manual taxon/flower review priority**: single-hue modes remain unverified |

For *Ochna serrulata*, the red sepals that persist after the petals fall and black mature fruits, versus earlier yellow petals, are independently documented by National Parks Singapore (https://www.nparks.gov.sg/florafaunaweb/flora/6/6/6680) and World Flora Online (https://www.worldfloraonline.org/taxon/wfo-0001087136). For *Cistus ladanifer*, true white petal morphs with or without maroon basal blotches are documented by Cistus & Halimium (https://cistuspage.org.uk/the-plants/cistus-ladanifer/) and Trees and Shrubs Online (https://www.treesandshrubsonline.org/articles/cistus/cistus-ladanifer/).

**Crucial inferential ceiling:** Even the visually clear yellow/purple *Viola* or yellow/white *Leptosiphon* contrasts could reflect taxonomic misidentification, environmental plasticity or cultivated introductions; no genotype-based or within-breeding-population morph frequency was measured.

## New source-exact post hoc numerical sensitivity

The three deliberately flagged taxa for the explicitly *hypothetical* exclusion stress test are **Argemone albiflora, Ochna serrulata, Oxalis tenuifolia**. Flags were chosen by author-side visual inspection **after** seeing low/central/high model scores; **there is no independent-blind provenance or pre-registration** and the excluded trio must not be treated as confirmed false positive ground truth. The remaining 24 taxa are not presumed to have valid petal morphs. The original 27 and 1,691 undetected controls remain intact.

Recalculate the mode-vs-month-adjusted continuous 50-km JSD spatial difference against all source controls, controlling for species-wide JSD, log source photo count, and log(1+50-km pair count), plus a second model adjusting for first PC variance and crude category richness. Use genus-cluster sufficient-statistic bootstrap (999 fixed draws), but report all intervals as exploratory. A CI that excludes zero **after visually selecting taxa to exclude** does not produce a newly confirmed effect.

CI uses the exact immutable source artifact CSV SHA:
- Mode+geometry 1,761-taxonomic records: f799b8b0c0f6a250e5ff0c1ad15bf3150854f3099a6a9c03ad5a26711ba1e7df
- Continuous nine-palette 1,761-species geography: 6dfb1b8ef5b9bf60d8903b97f27257c4b2ddb26e45231b9870032a8b1b5cf528

Code: scripts/analysis/fcp27_petal_organ_posthoc_sensitivity_20261023.py; output: source_conditional_organ_stage_sensitivity.json in source-linked Actions artifact.

## Next *actual* biological confirmation gate (now implemented, pending people)

The independent reviewer CSVs remain blank. Two independent human annotators must mark all original 594 photo IDs (including 'not_available' where images are missing), without opening the source PI mapping. A version-2 form adds:
- target_organ: corolla_petal / petaloid_sepal / sepal_bract / stamen / fruit / bud / vegetative / mixed / uncertain / not_available
- petal_pattern: uniform / basal_blotch / veined / bicolour_sector / irregular / unknown / not_available.

The original existing form must additionally specify: target_flower_visible, petal_region_unambiguous, observed_dominant_colour, flower_stage, lighting_white_balance, multiple_flower_or_organ_confusion, cultivation_or_taxon_concern, confidence. Both reviewers first freeze independent CSVs before PI joins on review_id. The adjudication code refuses to output species-level estimates on unfilled, incomplete or invalid data and writes 'human_review_incomplete_not_estimable' instead.

After two reviews, source strata (low/central/high PC1) are unveiled to the adjudication program; acceptable petal photos must be fresh and open, flower and petal clear, photographic hue assessable and reviewer labels concordant. Separately evaluate low vs high **whole-petal hue** and **basal/sector pigment pattern**, require at least two concordant images per photo-PC1 tail, report both candidate and matched-control groups and leave unstable taxa unknown rather than monochrome.

This procedure yields *validated photographic petal hue/pattern candidates* only. Genetic colour polymorphism, local maintenance, expression pathways, fitness and evolutionary pigment cost-benefit require other evidence. The earlier strict same-year/month <=5km different-observer null (P 0.300/0.275) remains unsupported.

Files:
- scripts/analysis/fcp27_two_reviewer_petal_adjudication_20261023.py
- .github/workflows/fcp27-organ-stage-falsification-20261023.yml
- Source 589-photo anonymous reviewer sets: run 39931967531, reviewer artifact 12492480632, PI mapping artifact 12492660357 (never distribute PI mapping to blind annotators)
