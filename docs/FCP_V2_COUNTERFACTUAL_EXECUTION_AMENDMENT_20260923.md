# FCP v2 counterfactual execution amendment

Date frozen: 2026-09-23 JST  
Status: **PRE-PIXEL EXECUTION AMENDMENT; NO FRESH IMAGE PIXELS OPENED**

This amendment refines compute allocation and pass ordering only. It does not change the FCP v2 biological questions, species queues, metadata gate, exposure levels, neutral background, ROI-jitter magnitudes, q_white, or any downstream inferential threshold.

## Why an amendment is needed

Running every expensive detector/segmenter counterfactual on all 40,000 target images would spend most compute on repeated segmentation rather than improving the measurement-validity estimand. The measurement question can be separated into:

1. all-row technical diagnostics and fixed-mask photometric sensitivity; and
2. expensive full-pipeline perturbation sensitivity on a large response-blind within-species subset.

The split is frozen here before fresh v2 image pixels are opened.

## Two-pass source-byte design

### Pass T — technical only

For all authorized v2 rows:

1. reacquire the exact frozen image URL;
2. record source-image SHA256;
3. run the frozen ROI-v4 detector/segmenter;
4. compute exposure, background and ROI technical diagnostics;
5. compute fixed-mask photometric technical diagnostics at EV {-1,-0.5,0,+0.5,+1};
6. for the heavy-counterfactual subset only, run full-pipeline exposure, neutral-background and prompt-jitter technical perturbations;
7. persist no morph, palette fractions, D, q_white projection, W, spatial colour outcome or H3 variable.

The complete technical table is sealed and read back before biological Pass B starts.

### Pass B — biological measurement

Reacquire the same frozen URLs only after Pass T is sealed.

For every successfully reacquired row:

- require the Pass-B image SHA256 to equal the Pass-T source SHA256;
- if bytes drift, mark the row source-drifted and do not substitute another image;
- run the frozen biological flower-colour measurement.

Pass B may then calculate original palette/morph outputs and the already frozen biological counterfactual outputs.

This avoids persisting a many-gigabyte image archive while preserving exact source-byte identity across the outcome firewall.

## Compute allocation

### All 40,000 terminal images

Run:

- base ROI-v4 measurement;
- flower/background exposure metrics;
- flower/background luminance and Lab summaries;
- mask area and horizontal-flip diagnostics;
- fixed-mask EV {-1,-0.5,0,+0.5,+1} transformations.

In Pass T, fixed-mask EV transforms produce only technical pixel summaries. Morph/palette classification from these transforms is deferred to Pass B.

### Heavy counterfactual subset

Target: **20 images per terminal species**, therefore 8,000 images if the 400-species metadata gate passes.

Selection is response-blind and frozen before pixel opening:

`heavy_hash = SHA256("FCP_V2_HEAVY_COUNTERFACTUAL_20260923" | panel | inat_taxon_id | photo_id)`

Within each terminal species, sort by `heavy_hash` and take the first 20 photo IDs.

The subset is fixed from metadata identity only. Exposure, background, ROI quality, colour, morph and biological outcomes are not used.

On these images run:

- full-pipeline EV -1.0, -0.5, +0.5, +1.0 relative to the unperturbed pass;
- full-pipeline fixed mid-grey background neutralization;
- the seven-prompt deterministic ROI-jitter set.

This yields 8,000 same-flower technical perturbation units while keeping the all-row base diagnostics at 40,000.

## Estimand allocation

### MV1 image-level invariance

Primary all-row component:

- fixed-mask exposure sensitivity over all terminal images.

Primary heavy-subset component:

- full-pipeline exposure sensitivity;
- background-neutralization sensitivity;
- ROI-prompt sensitivity.

Equal-species summaries are primary for the heavy subset.

### MV2 species-level D invariance

Primary D comparisons use all 100 rows/species for:

- unperturbed measurement;
- fixed-mask exposure counterfactuals.

Full-pipeline/background/ROI perturbations are reported as 20-row/species sensitivity summaries and are not substituted for the all-row D estimand.

### MV3 H2 geometry invariance

q_white remains fixed.

All-row W sensitivity is computed for:

- unperturbed measurement;
- fixed-mask exposure counterfactuals.

Heavy-subset W-like alignment summaries for full-pipeline/background/ROI perturbations are explicitly labelled subset sensitivities, not replacements for all-row H2.

### MV4 spatial measurement sensitivity

The original all-row v2 spatial statistic uses the terminal 100-row/species measurement.

Heavy-subset counterfactual spatial summaries are secondary because the 20-row design has different finite-sample geometry.

## No-rescue rule

After any Pass-T pixel is opened:

- heavy subset membership cannot change;
- 20 images/species cannot be increased/decreased;
- EV levels cannot change;
- neutral background cannot change;
- ROI-jitter set cannot change;
- a failed Pass-T technical row cannot be replaced;
- Pass-B source drift cannot be repaired by a replacement photo.

## Current-paper boundary

None of these v2 pixels or results are required for the current New Phytologist submission. v2 remains a separate prospective measurement-validity programme.
