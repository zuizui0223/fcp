# FCP Validation: photo-quality review queue, 2026-10-09

## Current status

The historical FCP 10-km source-only repeated-site/year/month test remains **HOLD**.
Two **retrospective public metadata-only** screens identified 15 gap-one and 10
gap-two Validation species that may have sufficient distinct new observer
photographs at the frozen site/month/year windows. These *are not* botanically
validated flowering photos and do not demonstrate local genetic morphs.

* Gap-one pinned Actions source: [37798611129](https://github.com/zuizui0223/fcp/actions/runs/37798611129) (15 of 22 species).
* Gap-two pinned Actions source: [37801452055](https://github.com/zuizui0223/fcp/actions/runs/37801452055) (10 of 60 species satisfying every required source year).
* Original source-only Validation anchored species: 10. Even under ideal independent photo-quality review, the ceiling is **35** species. Reaching the >=30 benchmark requires 20 of the 25 candidate species to fully qualify; 6 unsuccessful species would make that threshold unreachable.

## New fixed photo review queue

The machine-buildable source is `scripts/analysis/build_fcp_validation_photo_quality_queue_20261009.py`. It reads exactly the original public metadata from those two previously completed Actions artifacts; verifies the SHA256 of each input and the 263-taxon source queue; keeps every valid metadata photo for the 25 positively screened species; and checks source observer, taxon, year, month, site anchor, pairwise maximum anchor radius, position accuracy and licence.

Expected output: **25 species, 115 distinct candidate photo IDs and observation IDs** (39 from gap-one, 76 from gap-two). Each photo is `UNREVIEWED`: its flower presence, organ, focal taxon, imaging/exposure acceptability, colour measurability and individual identity are all deliberately unknown (`null`). No photograph is downloaded or inspected.

A gap-two species can pass only after **each** originally deficient year contains the required number of independently photographed, independently reviewed usable flowers. Passing one year cannot compensate for failure in another year. No taxon, year, locality or photo-ID may be replaced after inspecting colour results. The fixed untouched 2,000+730 candidate species are entirely outside this retrospective queue.

## Inference boundary

These metadata candidates are not confirmed flowers, genotypes, biochemical pigments, reproductive outcomes or site-specific climate anomalies. The maximum 35-species tally is **not** a qualified species count; confirmed new flower images currently remain zero. The original 10-km climate-related HOLD remains unchanged.

The next independent step, if explicitly conducted, is an image-level *quality* review of precisely the frozen 115 candidate IDs, preserving failures and exposure/organ ambiguity. A quality-positive image remains a photographic observation, not proven genetic polymorphism or adaptation.

## Reproducibility

`.github/workflows/fcp-validation-photo-quality-queue-20261009.yml` downloads only the two checked historical **GitHub Actions artifacts**, never contacts the iNaturalist API or opens image pixels. Input SHA checks, tampering regression tests, and output-denominator gates run before the queue is uploaded and frozen. Results are stored under `results/fcp_validation_photo_quality_queue_20261009/` only if validation passes.


## Repeated-year fragility diagnostic (metadata only)

The fixed 25 candidate species divide into **10 with zero spare photographs in one or more required source-year windows** (five gap-one, five gap-two) and **15 with at least one extra public-metadata photograph in every required year** (10 gap-one, five gap-two). This is a source-opportunity property, not a prediction of quality or genetic polymorphism. A species with a zero-spare year becomes unqualified under this fixed queue if any necessary photograph fails quality review; a spare-photo species may still fail if all alternates are unusable.

The extra file `candidate_species_quality_review.csv` reports the original month, source site-photo anchor, missing year cells, required observer slots, raw available metadata photo counts and minimum spare count. Review all 115 frozen candidate-photo IDs; do not modify the selected species set based on flower colour outcomes.

For the 30-species target, at least 20 of the 25 species must pass full year-specific photo quality checks. Therefore even if all 15 currently spare-supported species pass, at least five of the 10 no-spare species must also pass. This is a worst-case planning observation, not a measured quality pass probability.
