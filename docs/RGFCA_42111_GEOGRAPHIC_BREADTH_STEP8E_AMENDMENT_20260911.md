# RGFCA Step 8E amendment — separate species breadth from geographic breadth

Date: 2026-09-11 JST
Status: frozen after the metadata-only 42,111-species transport gate and **before any Step-8E image pixel or flower-colour outcome is opened**.

## Reason for this pre-outcome amendment

The fixed one-anchor-per-species frame successfully retains 42,111 species but its selected anchor locations occupy 124 of the 128 equal-area cells represented anywhere in the V1+V2 discovery metadata. A one-photo-per-species frame is appropriate for species-equal global colour-state composition, but it is not an adequate sole sampling frame for a geographic colour surface because a species observed in several cells contributes only its arbitrarily hash-selected anchor cell.

No Step-8E colour or image outcome has been opened. This amendment therefore changes the sampling representation **before outcome access**, not in response to a colour result.

## Two distinct breadth layers

### A. Species breadth

Retain the already frozen one-anchor-per-species frame unchanged:

- denominator: 42,111 species;
- one deterministic V1/V2 anchor per species;
- global summaries weight species equally;
- this layer estimates species-equal composition of the selected observed floral-colour states;
- it is not used by itself to map global geography.

### B. Geographic breadth

Construct an additional deterministic frame containing exactly one frozen discovery observation/photo anchor for every unique `(inat_taxon_id, equal_area_cell_id)` link present in the V1+V2 discovery metadata.

The current metadata-only census contains 85,337 unique taxon–cell links over 128 occupied discovery cells. For each taxon–cell pair, select the row with minimum SHA-256 of UTF-8:

`20260911|cell_id|inat_taxon_id|observation_id|photo_id`

after deduplicating exact observation/photo IDs across V1/V2. Flower colour, image pixels, climate, literature, and later biological outcomes cannot enter this choice.

No selected taxon–cell anchor may be substituted after Step-8E pixels open.

## Geographic estimator

For each occupied cell:

1. measure the fixed taxon–cell anchor for each represented species;
2. give each successfully measured species at most one contribution in that cell;
3. retain unresolved/acquisition/measurement failures in the cell denominator;
4. report cell-wise palette/four-morph composition with species-equal weights among measured anchors and explicit missingness;
5. do not weight a cell more strongly merely because it has more observations per species.

The global species-composition estimator continues to use the 42,111-species layer, not the 85,337 taxon–cell rows. The two layers answer different questions and must not be pooled into one frequency estimate.

## Missingness and transport

The geographic frame is frozen from already stored V1/V2 observation/photo IDs. Current URL resolution may fail for some frozen anchors. Such failures are retained and not replaced. Current image URLs, licences and attribution may be re-resolved metadata-only before pixels.

## Measurement model boundary

Both layers use the same already frozen RGFCA ROI-v4 / fixed-palette implementation. The existing independent Monarda localization-agreement non-support remains unchanged. Geographic scale and larger sample size do not convert the model output into verified focal-species petal truth.

## Claim boundary

- Species breadth: global composition of one fixed observed colour state per discovered species.
- Geographic breadth: global spatial distribution of fixed observed colour states across species–cell links.
- Neither layer alone estimates within-species polymorphism, modal species colour, C*, S*, adaptation, or an unbiased census of all angiosperms.

The prior 42,111-species protocol remains in force except where this amendment explicitly separates the geographic estimator from the species-equal estimator.
