# Polymorphism measurement missingness — Step 8 analysis result

**Reserve genus survives span + pure technical-failure control: `True`.**
**Reserve genus D_unbiased survives span + pure technical-failure control: `True`.**
**Reserve primary D-spatial survives span + pure technical-failure control: `True`.**
**Reserve flower-minus-background D-spatial survives span + pure technical-failure control: `True`.**

## What `technical failure` means

Only `not_evaluable_roi_or_flip_gate` is counted as technical failure. Ambiguous palette and no-biological-palette rows are kept separate.

## Measurement-process diagnostics

- discovery rho(D, technical failure) = **0.233427**
- discovery rho(D, ambiguous palette) = **0.727545**
- reserve rho(D, technical failure) = **0.217455**
- reserve rho(D, ambiguous palette) = **0.704373**

## Genus clustering after sampled span + technical failure

- discovery D: gain **0.158703**, p_lower **0.00914954**
- discovery D_unbiased: gain **0.156661**, p_lower **0.0110994**
- reserve D: gain **0.121211**, p_lower **0.0448478**
- reserve D_unbiased: gain **0.121764**, p_lower **0.0444978**

## Geometry-preserving D-spatial tests after sampled span + technical failure

- discovery primary: partial rho **0.126637**, p_upper **0.007**
- reserve primary: partial rho **0.099288**, p_upper **0.025**
- reserve flower-minus-background: partial rho **0.116241**, p_upper **0.01**

## Interpretation boundary

These tests separate clearly operational ROI/flip failures from ambiguity-related missingness. A surviving signal is not explained by sampled geographic span plus the rate of clear technical failures alone. Ambiguous palette rows remain unresolved measurement/biological ambiguity and are not treated as hidden morphs.
