# Sisyrinchium: paired-quadrat fruit set and two frequency axes (2026-11-26)

## Source
Takahashi, Takakura & Kawata (2016), PLOS ONE, DOI 10.1371/journal.pone.0164381.
Original CC0 data DOI 10.5061/dryad.j8q34. Verified public Zenodo mirror record 5003124, original data.zip 127502 bytes, MD5 008ec196aa5acb769dbfd36e835e32ba, SHA256 4859b803af66db1f8f18796fe08ce2581363258444de715f1d4851e9f2f988c5. Recovered from run 38327303795 artifact 11781395491, verified again by run 45778737837. Dryad API returned HTTP401; the original ZIP was not recreated or replaced.

## New analysis
The raw fitness.txt contains 601 unique records, with the author's population, year, line, quadrat and census fields in the same table. The separately archived frequency.txt has 512 quadrats; only 282 fitness records join by literal census-table keys, so no partial cross-table join was used. Original fitness-table census arithmetic was checked: niwa_t=v+w; ok=o+k; ok_p=ok/(ok+niwa_t) when the denominator is positive. Fruit bearing rate was reconstructed as fruits/(fruits+peduncles without fruits), excluding open flowers from the denominator, and checked against the author field. Filter original species code sp=n and color v/w, finite outcome and positive census/fruit-opportunity denominators. This leaves 570 records in 86 quadrats. Compare only quadrats containing observed fitness records for both morphs: 282 records, 32 quadrats, five sites (banpaku 6, ibusuki 5, meidai 7, nature2 11, ootaka 3). All census fields are constant within the analyzed quadrat.

Fit equal-weighted quadrat OLS:
delta_pp = 100*(mean white fruit-set rate - mean purple fruit-set rate) ~ 0 + population + 10*white conspecific frequency + 10*author ok_p + log(total count).
This contrasts the two morphs within quadrats and separates within-species color frequency and other-species relative abundance; population fixed effects and density are included. Shore distance was omitted due to substantial missingness (only 14 paired quadrats with nonmissing distance). This is a retrospective exploratory analysis, not a preregistered experimental test. HC3 uncertainty treats quadrats as units, not the individual plants, but does not correct possible residual spatial dependence.

## Observed results (independently calculated directly from run 45778737837 log tables)
| Predictor | Coefficient, percentage points | HC3 SE |
|---|---:|---:|
| White conspecific frequency per +10 percentage points | -0.8871916709145705 | 1.5164685295238218 |
| Author other-species frequency per +10 percentage points | -1.0676761597212447 | 2.0299706813305414 |
| Natural log of total plant count | 12.19954166515451 | 14.239136919855978 |

Residual df 24. Neither frequency slope is distinguishable from zero with this exploratory model. Both negative directions are consistent with rare-white advantage and other-species interference, respectively, but do not establish either or their independence. The observational sampling time/census may differ from the collected fitness individuals. A negative white-frequency slope does not alone establish stable polymorphism, invasion criteria, genetic mechanism or lifetime fitness. This result does not overturn the original study's analysis, which differs in estimand.

Reproduction script: scripts/analysis/fcp_sisyrinchium_paired_fitness_20261126.py.
CI: first run 45779001343 failed importing an incompatible statsmodels/pandas dependency; direct NumPy/SciPy HC3 replacement in ce3cf24dac0cc359f043a2b5e9d73d09b54da6fc retains the same model and asserts the two observed coefficients. Reverification run 45779101934.
