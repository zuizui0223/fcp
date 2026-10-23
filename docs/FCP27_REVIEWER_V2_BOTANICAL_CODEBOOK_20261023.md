# FCP27 independent photograph reviewer codebook v2

This file describes the **two independently blinded human photographic reviewers**. Their observations are not inferred by the assistant or the original model. Original image source: 589 successfully recovered photographs and 5 explicitly missing source images, 54 anonymous taxa total.

## Distribution

Give each reviewer the anonymous source-photo package and that reviewer's reviewerN_botanical_v2_blank.csv only. The PI source mapping, original taxon names, photo-PC1 group, machine hue, and candidate/control arm must remain hidden from both reviewers until their first-pass responses are frozen. Photo thumbnails are individually anonymous and bundle_id groups an unknown species. Do not edit review_id or bundle_id. Respect image CC licences and photographer credit before wider distribution.

Fill **all 594 rows**, including explicitly unavailable images. Use these exact lowercase machine-readable spellings, not translations.

| Column | Allowed values |
|---|---|
| target_flower_visible | yes / no / uncertain / not_available |
| petal_region_unambiguous | yes / no / uncertain / not_available |
| target_organ | corolla_petal / petaloid_sepal / sepal_bract / stamen / fruit / bud / vegetative / mixed / uncertain / not_available |
| observed_dominant_colour | white / cream / yellow / orange / red / pink / magenta / purple / blue / green / mixed / other / unassessable |
| flower_stage | open_fresh / bud / senescent / fruit_postfloral / mixed / unknown / not_available |
| petal_pattern | uniform / basal_blotch / veined / bicolour_sector / irregular / unknown / not_available |
| lighting_white_balance | no / acceptable / yes / uncertain / not_available |
| multiple_flower_or_organ_confusion | yes / no / uncertain / not_available |
| cultivation_or_taxon_concern | yes / no / uncertain / not_available |
| review_confidence | high / medium / low / not_available |
| reviewer_notes | Free text explaining uncertainty, flower organ, bud, sepal, fruit, camera condition or potential misidentification |

Definitions: target_organ identifies the organ **responsible for the colour**, not merely whether a flower appears somewhere. Petaloid sepals may be scored as a target perianth when botanically appropriate, but reddish post-floral sepals, stamens, bracts, fruits and unopened buds must be identified separately. Petal pattern means basal blotch, venation or sector colour within an open flower, not a distinct entire-petal hue. A lighting_white_balance response of no or acceptable denotes no material apparent photographic colour problem. Mixed organ evidence means a single interpretable floral hue cannot be isolated. A photograph alone usually cannot establish wild provenance; use uncertain for the cultivation flag unless there is direct visual evidence.

For the five unavailable source images, set target_flower_visible, petal_region_unambiguous, target_organ, flower_stage, petal_pattern, lighting_white_balance, multiple_flower_or_organ_confusion, cultivation_or_taxon_concern, and review_confidence to not_available. Set observed_dominant_colour to unassessable. Enter "image unavailable" in reviewer_notes; do not replace a missing original photo.

## Independence and analysis

Each reviewer must annotate without consulting the PI mapping or the other reviewer's answers. After both files are complete, the source mapping is joined by review_id. Source SHA256 and both reviewer CSV SHA256 values are retained in the adjudication result.

The FCP27 adjudication script requires two concordantly rated, fresh open flowers with a clear comparable petal or petaloid perianth target and acceptable lighting **within each source low- and high-PC1 group** to record an *observational photographic corolla hue candidate*. Petal pigment patterns are analyzed separately. If validity or concordance is insufficient the response is not_estimable, not monomorphic. Even a convincing photograph-supported difference does not by itself demonstrate inherited genotype morphs, reproductive fitness or local population polymorphism.

The code prevents blank or incomplete files from being counted as zero variation: without two complete human review sheets the only allowed result is human_review_incomplete_not_estimable. No independent human annotations have yet been completed in this project.
