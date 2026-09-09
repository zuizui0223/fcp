# RGFCA global flower–background signal recurrence exploration

Date frozen: 2026-09-09 (JST)

## Status and purpose

This is an **exploratory, hypothesis-free signal-discovery lane**. It does not prespecify island whitening, a common flower-colour boundary, a pollinator syndrome, or the direction of any colour change. The previously opened pilot census/PCA is retained only as feasibility evidence; recurrence quantities below were defined before they were computed.

The measurement source is the frozen 500-species RGFCA reserve artifact. No image reacquisition or remeasurement is authorized by this protocol.

## Measurement boundary

The analysis uses photo-derived flower and matched-background 12-anchor palette counts. Because target-domain flower-region localization remains imperfectly validated, these quantities are reported as **photo-derived flower–background signal structure**, not calibrated animal-perceived colour or demonstrated pollinator attraction.

## Eligible photos

A row is eligible when all of the following hold:

1. `roi_status == automated_colour_state_admitted`;
2. `global_classifiable == true`;
3. flower 12-anchor palette mass > 0;
4. matched-background 12-anchor palette mass > 0;
5. coordinates are present;
6. positional accuracy <= 5000 m.

No island/mainland label, climate variable, latitude band, biogeographic realm, pollinator variable, or outcome colour is used to define eligibility.

## Signal vector

For each eligible photo, define the 12-dimensional matched signal vector

`d = flower_palette_fraction - background_palette_fraction`

in the fixed order:

`white, yellow, orange, red, pink, magenta, purple, blue, bronze, green, brown, black`.

Absolute flower colour is retained as a descriptive companion only. The recurrence analysis uses the matched signal vector.

## Local-state construction

Each realization uses EPSG:6933 equal-area coordinates and a 300 km square grid. The grid origin is independently shifted in x and y by Uniform(0, 300 km).

A species × cell local state is retained when it contains at least 3 eligible photos and at least 2 unique observers.

Within each local state:

1. compute each observer's mean signal vector from that observer's photos;
2. bootstrap observers with replacement, drawing the observed number of unique observers;
3. for each drawn observer, bootstrap one of that observer's photos;
4. average the drawn photo vectors equally across drawn observers.

Thus prolific observers do not dominate the local-state estimate.

A species enters a realization only when it has at least two retained local states whose centroids are separated by at least 500 km. Island/mainland membership is irrelevant to this gate.

## Species-equal signal decomposition

For each admitted species, subtract the species mean local-state vector from every local-state vector. If a species has `k` admitted local states, each of those rows receives weight `1/k`, so every species contributes total weight 1.

Compute the weighted covariance matrix of the 12-dimensional within-species state deviations and its eigen-decomposition.

The first five eigenvectors are retained as exploratory signal modes. Eigenvector sign is arbitrary and is aligned only for recurrence comparison; no biological direction is assigned in advance.

## Recurrence estimands

Use 200 realizations, master seed `20260909`.

For each mode rank 1–5 report:

1. distribution of explained variance ratio;
2. absolute cosine alignment to the full-data reference eigenvector of the same rank after optimal one-to-one matching among the first five modes;
3. recurrence probability that alignment is >= 0.80;
4. recurrence probability that alignment is >= 0.90;
5. frequency with which each palette coordinate is among the two largest absolute loadings of that matched mode.

Because adjacent PCs may swap order, mode matching is performed within the first-five subspace by maximizing total absolute cosine similarity over all 5! permutations. Signs are then aligned to the matched reference vector.

The reference decomposition uses the same 300 km grid with zero shift and observer-equal local-state means without bootstrap. It is a coordinate system for comparing realizations, not an independent validation set.

## Robustness summaries

Across realizations also report:

- number of admitted species;
- number of admitted local states;
- palette-coordinate species-equal mean absolute within-species deviation;
- first-five subspace similarity using principal angles.

No threshold is defined for declaring a biological discovery. These are recurrence/stability summaries for exploratory mode discovery.

## Ecological overlay hard stop

Until the recurrence outputs above are frozen, do **not** use island/mainland status, island area, isolation, climate, latitude, realm, pollinator guild, reproductive system, or other ecological covariates to choose, rotate, relabel, or retain signal modes.

After recurrence is frozen, ecological variables may be overlaid as a separate interpretation stage. Any association found there is exploratory unless independently replicated.

## Claim boundary

This lane may support statements about recurrent **photo-derived flower–background signal structure** across geography within species. It does not by itself establish:

- pollinator-perceived contrast;
- attraction or visitation effects;
- adaptive evolution;
- an island syndrome;
- a universal geographic boundary;
- causal effects of climate, geography, or pollinators.
