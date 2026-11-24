# FCP 2026-11-24: Source-authenticated geographic open-petal holdout review packet

**Status:** heldout source photo IDs were sampled and forms created; **no humans have scored these photos**. Numerical FCP results still concern photographed flower-region colours, not independently verified hereditary open-petal morphs.

## Why a second review packet is necessary

The original 128 near-pair pilot was selected from a temperature-tail frame and oversampled original colour differences. It cannot validate the 3,600-species 50-km geography effect, nor can 66 readable one-AI-reader pairs establish population-wide classification accuracy. The new source-locked sample is a pre-analysis *blind scoring design* for a geographically relevant external organ check.

## Source and sampling

- Original [authenticated release](https://github.com/zuizui0223/fcp/releases/tag/fcp-all4730-original-method-actual-colour-results-20261017): archive SHA256 `99f7ee68c3fd009e8ab60fd0cfcdf1fb407e3a4e885cb7be49747c`; original measured 447,500 photo CSV SHA256 `8e7f10f0427927208485b2da6fa9c46bb681d6f3611e53b8dfdc606cfc8dabe6`; 4,475 taxa, 227,107 source photo-colour classifiable.
- One *near* source-colour pair (<=50 km) and one *far* pair (>=250 km) within a species. Both original photographers must differ, both geolocation uncertainties <=5 km, and photo records must be from the same calendar month (possibly different years). Species must have eligible opportunities in both geographic distance bands.
- Original same/different photo classes define **four ascertainment strata** for paired near/far source matches. There were 3,150 eligible species for near-same/far-same; 2,859 for near-same/far-different; 2,306 for near-different/far-same; 2,271 for near-different/far-different. **These are overlapping eligibility sets, not four separate species cohorts.**
- Precisely **24 source taxon draws per stratum = 96 taxon-stratum draws**, each contributing two distance bands, hence **192 near/far pairs**. Sampling uses deterministic seed 20261124, species uniform within each source stratum and a source pair uniform within that species and source stratum. Source-cell counts and conditional pair inclusion probabilities are retained in the investigator's key. The selected samples should never be treated as a naturally representative global photo-label distribution.
- An independent photographic label-quality panel contains 64 source-classifiable and 64 source-unclassifiable photos (one per sampled taxon and classifiability stratum), for **128 single-photo assessments**. This assesses whether apparent open petals were present among source classification failures as well as successes; source status is concealed from reviewers.
- **501 distinct photo IDs** across pair and singleton draws. We excluded the 250 distinct original pilot photo IDs and four original published example photo IDs; overlap with those 254 originals is **zero**. Selected pairs may share an image or species across source strata, so analyses must preserve pair/species identity and never count every raw photo as independent.

## Frozen artifacts and reproducibility

- [Verified source-run and quality gates](https://github.com/zuizui0223/fcp/actions/runs/45378374319) (original measured source SHA256 check, source criteria, counts, blind forms, and separated unblinding artifacts).
- **Reviewer-facing blank forms only:** [four reviewer CSVs](https://github.com/zuizui0223/fcp/actions/runs/45378374319/artifacts/14911677807) — reviewer A and B, each with geographic-pair and classified/unclassified single-photo queues, ordered differently, with only original photo-page links, comparable open-petal checks, and empty botanical scores.
- **Investigator-only:** [separate decode key and source opportunity ledger](https://github.com/zuizui0223/fcp/actions/runs/45378374319/artifacts/14912142270) — source photo labels, near/far status and distances, species, observational opportunities and selection probabilities. **Do not share this artifact with reviewers before their independent scores are frozen.**
- `scripts/analysis/fcp_geo_petal_heldout_sampling_20261124.py` regenerates the packet from the byte-identical original release. Result receipt: `results/fcp_geo_open_petal_heldout_sampling_verified_20261124.json`. Two-reader score validator: `scripts/analysis/fcp_geo_petal_two_human_review_20261124.py`; it returns `NO_HUMAN_PETAL_RESULT` when reviewer sheets are blank and refuses to compute apparent petal effects before independent reviewers (and a third adjudicator for disagreement) have submitted actual observations.

## Reviewer instructions and interpretation

Score **homologous fully open petals** only. Explicitly indicate whether each image visibly contains a resolvable open corolla and whether both images show comparable developmental stages. Exclude uncertain organs or stages (buds, fruits, sepals, stamens/anthers, pollen, inflorescence background); never guess a colour when the image is unresolvable. Assign the four predeclared source-compatible colour groups **only** for readable open petals. Reviewers must not see source original colour, distance category, cohort, prior AI visual categories, or selection stratum. The original photo page may still display species or metadata, so this is *label/distance-blind*, not completely metadata-blind.

After actual independent scoring, report source-stratum-specific open-petal far-minus-near mismatch contrasts, the number of taxon-stratum comparisons with **both** geographic pairs readable, unresolved-outcome bounds, and source-classifiable vs source-failed photo readability. Keep full A and B scores before adjudication. If both reviewers disagree on petal status/colour, require a third botanical reviewer. Do not pool the source-enriched strata into global 3,600-species biological rates without a validated poststratification or species-weighted refit.

The source geographic comparison (<=50 km vs >=250 km) is intentionally narrower than the main FCP range-wide-minus-50-km photographic estimator. A successful visual review could support inference for this **dual-opportunity source frame**; it cannot by itself establish genetics, pigment synthesis costs or pollinator-mediated selection.

## Decision

**Ready:** independent reviewer packets, original source opportunity counts, inference safeguards and reproducible sampling.

**Still blocked by external data:** two genuine independently scored botanical review files (and adjudication where they disagree). Synthetic CI scores are never scientific evidence. Until valid human labels exist, do not elevate photographed colour to verified flower petal colour in the main manuscript.
