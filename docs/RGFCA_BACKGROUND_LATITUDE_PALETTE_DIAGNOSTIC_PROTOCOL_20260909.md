# RGFCA background-latitude palette diagnostic — post-hoc descriptive protocol

Date frozen: 2026-09-09 JST

This diagnostic is opened only after the ecological overlay showed no calibrated flower-only latitude association for M1–M3, while matched-background latitude slopes were strongest for M1 and M2. It is explicitly post-hoc and cannot upgrade those exploratory findings to confirmatory evidence.

## Fixed question

Which of the 12 retained palette fractions contribute to the observed within-species latitude structure in the flower (`F`), matched background (`B`), and flower-minus-background (`D = F - B`) components?

## Fixed data and realizations

Reuse the exact 200 recurrence realizations, grid shifts, observer/photo bootstrap, admitted species/local states, and absolute-latitude predictor from the ecological overlay frozen at `106ffd9986a5a2316a0909df14fb26f0dee28296`.

No island, climate, candidate-cell, or other ecological variable enters this diagnostic.

## Estimator

For every palette colour (`white`, `yellow`, `orange`, `red`, `pink`, `magenta`, `purple`, `blue`, `bronze`, `green`, `brown`, `black`) and each component F/B/D:

- compute species-centred local-state palette fraction;
- use absolute latitude per 10 degrees, centred within species;
- use the same species-equal state weight `1/k_s`;
- fit the same through-origin slope `beta = sum(w*x*y) / sum(w*x^2)`.

Reuse the ecological-overlay latitude null: within each admitted species and realization, permute latitude across retained states before fitting the same slope. Use null seed `20260912` for this diagnostic.

Report across 200 realizations the observed median and 2.5–97.5% quantiles, fraction positive, null median and quantiles, and the same two-sided p-like calibration. Report all 12 colours; do not select only colours that look strong.

## Claim boundary

This is a descriptive decomposition of an already outcome-opened background signal. It does not establish ecological causation, adaptation, pollinator perception, or a universal latitude rule. It is intended to identify whether the matched mode-level pattern is driven by green/white, brown/yellow, or other palette structure and to flag potential photographic/background nuisance mechanisms.
