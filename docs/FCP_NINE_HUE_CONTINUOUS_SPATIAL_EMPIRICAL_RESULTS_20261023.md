# Original-photo empirical result: continuous nine-hue flower-colour geography

**Authenticated original-source execution complete: 23 October 2026.**
Workflow: https://github.com/zuizui0223/fcp/actions/runs/39918785955
Machine receipt: results/fcp_all4475_continuous_spatial_verified_20261023.json
Analysis PR: https://github.com/zuizui0223/fcp/pull/260

## Decision 1: geographic hue organisation is not an artifact of four predefined classes

The complete source has exactly 447,500 original photographed records from 4,475 taxa; 227,107 are original four-state-classifiable. The new test uses their original nine-element continuous palette compositions, base-2 Jensen-Shannon divergence, fixed observed coordinates and 199 permutation worlds per null condition. It retains 1,761 species under the frozen >=40-source-photo and >=30 50-km local-pair gate. The exact original four-state strict 50-km mean +0.0176149021685411 is independently reproduced before the new results are admitted.

**Calendar-month-preserving continuous palette-vector permutation:**

| Spatial radius | Evaluable taxa | Observed minus matched-month null, equal-species | Upper permutation P | Exploratory Holm across five radii |
|---|---:|---:|---:|---:|
| 5 km | 372 | +0.013337806 | 0.005 | 0.025 |
| 25 km | 1,232 | +0.012525764 | 0.005 | 0.025 |
| 50 km | 1,761 | +0.011456425 | 0.005 | 0.025 |
| 100 km | 1,761 | +0.009116586 | 0.005 | 0.025 |
| 250 km | 1,761 | +0.006334573 | 0.005 | 0.025 |

The original unconditioned source-vertex permutation gives mean positive continuous excess +0.013684417 at 50 km. Restricting to one photo per observer and recomputing the month-matched null retains 1,655 species at 50 km, with +0.011233602, P=0.005, Holm P=0.025. A photo-colour pattern confined to a month was correctly rejected in the synthetic test.

**Ecological inference:** on this source's photo-colour scale, conspecific flower-palette resemblance is geographically nonrandom, even without four hard category boundaries and after calendar-month composition matching. For the common 1,761 species at 50/100/250 km, the positive local excess decreases with expanding radius. This is not yet proof of isolated populations or independent spatial hue clusters.

A crucial negative result remains untouched: the prior <=5-km **same-year/same-month plus independent-observer coarse-colour test** was unsupported (P=0.300 with >=5 eligible pairs, P=0.275 with >=30 pairs). The new 5-km calendar-month test has substantially different selection and conditioning, and must not be presented as reversal of that stricter null.

## Decision 2: number of photographed colour labels does not add a proven spatial rule

With a minimum of five geographically valid photos supporting each original coarse label, and the strict 50-km opportunity gate:

| Recorded labels supported by >=5 photos each | Species | Mean 50-km month-adjusted nine-hue geographic excess |
|---|---:|---:|
| One colour | 948 | +0.006119440 |
| Two colours | 674 | +0.015638228 |
| Three colours | 136 | +0.028391194 |
| Four colours | 3 | -0.009277880 |

These are label support counts, not the number of inherited morphs or GMM modes. Raw positive two-versus-three difference: +0.012752966.

After adjustment for the species-wide continuous nine-hue divergence, log classified photo count and log(1+50-km local photo-pair count, the **three-minus-two adjusted contrast** is **-0.004520111**, 499-replicate genus cluster bootstrap interval **[-0.015969320, +0.006441172]**, 447 genera in 810 admitted taxa. Thus greater *recorded category richness* does not independently predict stronger geographic hue sorting in this exploratory linear adjustment.

One-photo-per-observer sensitivity: 604 two-colour and 120 three-colour taxa; adjusted difference **-0.004548704**, genus bootstrap interval **[-0.017780217,+0.007498755]**, 401 genera. Same qualitative conclusion. Four-colour strict support has only **three taxa**, so it is not adequate for a four-morph geographic population claim or even stable cross-species group effect.

Interpretation: it is useful to retain distinct estimands for (i) continuous species-wide palette diversity, (ii) the direction of hue change including white-chromatic contrast, and (iii) the spatial organization of palette differences. The total number of coarse assigned labels cannot stand in for all three.

## Technical and inference boundaries

- This is an expressly **post hoc, single-provider** analysis on selected and classified photos; 220,393 source images lack an original four-state label and are not assigned colours here.
- Source hue vectors reflect camera/image segmentation and nine-colour palette shares, not spectrally calibrated reflectance or species-wide CIELAB. Illumination, developmental stage, flower-organ composition, cultivation, regional observer bias and taxonomic confusion may still generate some observed similarity.
- Calendar-month label/vector preservation is **not** a joint location-by-year-by-month independent-observer null. A positive month-conditioned result cannot imply genetic within-population morph coexistence.
- A continuous geographic autocorrelation signal does not reveal whether the true floral distributions are unimodal gradients, discrete modes, mosaics or transition zones. K=2 k-means from historical H2 was an axis-estimation device, not evidence of two natural flower-colour morphs. Discrete morph numbers require independent measurement-validation and stability tests.
- The 199-world plus-one permutation P=0.005 is the lowest attainable resolution, not an estimate of a finer P value. Holm adjustment is within the exploratory radius family; neither genetic inheritance nor adaptive selection was tested.
- Genus-cluster bootstrap addresses shared-genus dependence partially, not a full dated phylogeny or MNAR measurement missingness. Conditional model adjustments cannot establish causal effects of the number of photographic categories.

## Repository evidence

1. GitHub Actions source execution: https://github.com/zuizui0223/fcp/actions/runs/39918785955
2. Machine-readable checked receipt: results/fcp_all4475_continuous_spatial_verified_20261023.json
3. Original source photo archive SHA256: 99f7ee68c3fd009e8ab60fd0cfcdf1fb407e3a8de07edd38341b27889310bc25
4. Authentic compressed photographed CSV SHA256: 8e7f10f0427927208485b2da6fa9c46bb681d6f3611e53b8dfdc606cfc8dabe6
5. Workflow artifact with two per-taxon tables and full summary: artifact 12485608012, ZIP SHA256 79515908d7e46f0af4c14f490c1848fee4bed68ea580966fcc42141a78b8d24e.

**Bottom line:** a continuous geographically organised photographed phenotype is supported. A higher discrete-colour-type count *independently causing or predicting* stronger spatial differentiation is not supported. True discrete vs continuous morph biology remains an open measurement and genetic question.
