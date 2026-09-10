# Polymorphism genus clustering — Step 7b adjustment decomposition result

**Decision: `SPAN_NOT_SUFFICIENT_BUT_SIGNAL_ENTANGLED_WITH_CLASSIFICATION_YIELD`.**
**Reserve raw D_unbiased replication: `True`.**
**Sampled span alone insufficient: `True`.**
**Reserve genus clustering survives n_classifiable-only adjustment: `False`.**

## Discovery

- raw D_unbiased: gain **0.204285**, p_lower **0.00164992**
- span-only D: gain **0.190610**, p_lower **0.00259987**
- span-only D_unbiased: gain **0.189455**, p_lower **0.00289986**
- n_classifiable-only D: gain **0.165441**, p_lower **0.00499975**
- n_classifiable-only D_unbiased: gain **0.164771**, p_lower **0.00484976**

## Species-disjoint reserve

- raw D_unbiased: gain **0.153907**, p_lower **0.020749**
- span-only D: gain **0.141642**, p_lower **0.0229989**
- span-only D_unbiased: gain **0.143566**, p_lower **0.0234988**
- n_classifiable-only D: gain **0.032849**, p_lower **0.313234**
- n_classifiable-only D_unbiased: gain **0.035106**, p_lower **0.301235**

## Fixed-denominator diagnostic

- raw photo denominator per species: **100 in discovery and 100 in reserve**
- discovery rho(D, n_classifiable) = **-0.537139**
- reserve rho(D, n_classifiable) = **-0.555901**
- discovery rho(D, sampled span) = **0.181031**
- reserve rho(D, sampled span) = **-0.002567**

## Shared 23 repeated genera — D_unbiased

- Spearman rho = **0.521739**, p_two_sided = **0.0109495**

## Interpretation boundary

A span-only surviving result rules out sampled geographic span as a sufficient explanation of the reserve genus pattern. Failure after n_classifiable adjustment is treated separately because n_classifiable is a post-acquisition classification-success variable under a fixed 100-photo denominator and is strongly entangled with mixed/uncertain colour measurement. Non-random classification missingness therefore remains a measurement limitation rather than a solved confound.
