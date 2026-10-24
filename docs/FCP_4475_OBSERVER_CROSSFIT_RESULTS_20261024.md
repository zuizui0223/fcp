# FCP original photo-colour × thermal gradient: independent-observer crossfit results

**24 October 2026. Completed original-source, observer-disjoint photographic cross-fitting and independent audit, stacked on PR #270.**

## Scientifically relevant conclusion

The original 447,500-photo FCP atlas (4,475 named taxa; 227,107 classifiable nine-palette photo ROI measurements) showed a recurrent species-equal positive thermal-context/PHOTO-white colour association and an exploratory stronger association in taxa with both white and nonwhite original photo colour labels.

That species-level **modifier claim required new scrutiny**: the original white/nonwhite baseline group and temperature/white slope had used the same source photos. The new study assigned original observers—not photo IDs independently—to two outcome-blind folds inside each named species, used the white photo baseline only from one fold and recalculated temperature/white partial-slope-like coefficients from the *other* disjoint observer fold, then reversed the roles. Four variants were analysed: forward/reverse observer folds, each with all eligible heldout source photos or exactly one photograph per heldout observer.

**Main new result:** The positive *across-species mean* source PHOTO–thermal white association is supported by heldout-only observer data in both fold directions. However, the *white/nonwhite photographed coexistence modifier* is not robustly separated from zero in the observer-disjoint, source opportunity-adjusted models. Moreover, *individual named-species thermal slope estimates* have almost no crossfold agreement at the available half-sample size. The correct inference is therefore a recurring average SOURCE PHOTO climate-hue association, **not a validated species-specific thermal phenotype selection gradient, genotype morph frequency rule or pigment cost-benefit explanation**.

## Source verification and observer separation

Immutable original source: 447,500 photograph IDs / 4,475 species, 227,107 originally classifiable nine-colour flower ROI photos. Compressed original source CSV SHA256: 8e7f10f0427927208485b2da6fa9c46bb681d6f3611e53b8dfdc606cfc8dabe6. Original source WorldClim 2.1 BIO archive SHA256 00513224583665ec0f2f955a4ec252730c4deb2004cce9e793492a3f26df4dcf; monthly SRAD source archive SHA256 c72ee7f4f9a0eb4b5f6cd7a003eddc05bda22a2e5666d968ed4e817fd36b9026.

Original run: https://github.com/zuizui0223/fcp/actions/runs/40080001765 (**PASS**). Full taxon-level original source artifact: 12554693263, ZIP sha256 50d5d15a2d145b0db7a176f945435701251768ba614f659d9dac496011b33177. Source CSV and original result JSON independently SHA-verified in a second run: https://github.com/zuizui0223/fcp/actions/runs/40080546906 (**PASS**).

Photo fold allocation was based on observer IDs from ALL original photographs within each species (before photographic colour/classification or climate filtering), with deterministic SHA256 observer ranking. No photographer's photographs in a given species appear on both baseline and heldout sides; missing observer IDs are excluded from this test. Baseline photos require >=10 classified images and >=8 observers. Heldout requires >=18 valid complete climate/geo source images, >=10 disjoint observers, >=3 original WorldClim 10-arcminute raster cells, non-collinear shared BIO1/BIO5 temperature contrast conditioned on precipitation, climate seasonality, solar, broad XYZ geography and photo month. Source month-matched photographic nulls have 199 real-photo hue permutations; whole taxonomic genera are bootstrapped 999 times. **All 4,475 source taxa remain explicitly in each of the four panels' ledgers** even when half-sample photo/climate support is not identifiable.

## Result 1: source photographer-independent species-wide temperature–white direction remains positive

| Heldout observer PHOTO sampling | Species with eligible source response | Mean observed temperature–white photo slope | 999-genus marginal 95% | Fraction eligible species with positive slope |
|---|---:|---:|---:|---:|
| Baseline fold A, test fold B | 2,403 | **+0.012879** | **[+0.007331,+0.018743]** | 54.8% |
| Baseline fold B, test fold A | 2,430 | **+0.005747** | **[+0.000118,+0.011283]** | 51.6% |
| A -> B, one eligible original photo per heldout observer | 2,324 | **+0.011011** | **[+0.005585,+0.016732]** | 54.1% |
| B -> A, one per heldout observer | 2,364 | **+0.006159** | **[+0.000006,+0.011839]** | 51.7% |

These four observational CI results are *separate exploratory marginal bootstrap intervals*, not joint familywise controlled new pre-registered confirmation. The reverse one-photo-per-observer interval barely excludes zero, and classification/camera/platform/geography can remain confounded. At a population mechanism level these data are still photographs, NOT genetics, pigment chemistry or field physiology. The paired fold directions come from the same source pool; they are not two fully independent biological studies.

## Result 2: the supposed stronger photo-colour coexistence effect is uncertain in independent-observer comparisons

The main baseline PHOTO-group contrast was originally defined as >=3 baseline photographs labelled white AND >=3 labelled nonwhite, versus no baseline photograph labelled white. All baseline photographs are from different observers than those supplying each target species' thermal-white photo response. A second design raises the threshold to >=5 of each class to match the earlier source study's support definition.

The adjusted source-group model controls source baseline/heldout photo sample counts, georeferenced climate grid cells, heldout BIO5 photographic q90–q10 climate breadth, residual identifiable warmth-climate variance and original sampling cohort. It uses whole-genus 999-resample marginal CIs plus a 199-world heldout photo-month matched null.

| Observer-disjoint white/nonwhite photo-group contrast | Eligible taxa | Adjusted coexistence-minus-no-white coefficient | Genus marginal95 | Photo-month-null p (two-sided) |
|---|---:|---:|---:|---:|
| A baseline -> B response | 2,403 | +0.010504 | [−0.004635,+0.026140] | 0.180 |
| B baseline -> A response | 2,430 | +0.005819 | [−0.009495,+0.021990] | 0.410 |
| A -> B with one PHOTO per heldout observer | 2,324 | +0.006301 | [−0.009971,+0.022733] | 0.455 |
| B -> A with one PHOTO per heldout observer | 2,364 | +0.006086 | [−0.009772,+0.023271] | 0.445 |

All four group point estimates retain a positive sign, **but none is robustly separated from zero** under the tested source conditions. Requiring >=5 baseline photographs for both photo labels does not repair the non-support; adjusted genus CI includes zero and the monthly photographic null also does not reject for all four photo panels. Do not convert a strong group mean in one fold into evidence that the group-to-group difference or a genetic polymorphism maintenance mechanism has been established.

The previously reported within-source same-photo-group adjusted coefficient (+0.008141 on 2,837 taxa, marginal genus95 [+0.000615,+0.015923], observer thinning interval crossing zero) remains an **exploratory photographic association** and is neither deleted nor silently promoted to independent biology. Fold-based sample sizes and colour grouping thresholds differ from earlier source analyses.

## Result 3: baseline original-photo class states are moderately repeatable, but per-species climate coefficients are unstable

On 4,246 species with sufficient photographed baseline observations in **both independent observer photo halves**, four-class photo white-state agreement between source folds was **69.05%**, chance agreement from marginal group abundances **33.47%**, yielding Cohen's **kappa=0.535**. It measures stability of the *observed-photo label composition*, not of biological petal-colour genotype and not of baseline source instrument calibration.

On 1,919 species with separately estimable heldout thermal/white slopes in both directions, their source PHOTO partial slope correlation was **Pearson r = +0.0181**, 999-whole-genus CI **[−0.0374,+0.0767]**. After one original photo per heldout observer, 1,843 source species remained, crossfold **r = −0.0061**, CI [−0.0641,+0.0575].

**Interpretation:** Species-specific source PHOTO/climate coefficients are poorly repeatable after halving the data. That could reflect sampling error and reduced within-species climate breadth, measurement noise, genuinely heterogeneous/unstable source photo conditions, or actual biological variability. It does not demonstrate a universal absence of response or rule out meaningful ecological thermal colour relationships.

A separate exploratory baseline continuous photo whiteness quartile comparison gave a similar intermediate-white photo peak in the forward heldout split, but reverse-fold quartile means differed; this pattern is not taken as a validated evolutionary optimum.

## Results and manuscript framing

> Across many named flowering-plant taxa, a small species-mean source photographic white-fraction response to shared habitat warmth remained positive when the photo-based temperature association was recalculated using observers independent of those defining baseline white states. The proposed amplification of this response in taxa with co-observed white and coloured PHOTO classes did not survive as a well-identified effect in observer-disjoint, photo-opportunity-adjusted cross-fitting. Although photo-state category assignments were moderately concordant across observer halves, individual species' half-sample temperature–photo-hue coefficients were not reproducible. These results support a shared average photographic colour–climate geometry and delimit, rather than establish, its biological morph-maintenance explanation.

All findings are bounded by source original photo-colour classification missingness (220,393/447,500 records), flower-age/organ/camera/background/cultivation uncertainties, historical long-term source WorldClim normals and no individual genotype/fecundity data. The earlier independently unsupported flower-minus-background BIO5 source validation (rho −0.00328, p=0.82) remains a material caveat. The source's strict <=5km same-year/month observer-independent photo-geographic null remains a separate unsupported sensitivity and is not made an arbitrary universal pass/fail gate.

## Reproducibility

- Original genuine source observer-disjoint PHOTO-fold estimator: scripts/analysis/fcp4475_crossfitted_observer_independent_thermal_white_20261024.py.
- Source-release authenticating CI: .github/workflows/fcp-4475-crossfit-observer-thermal-white-20261024.yml; passing full genuine source run https://github.com/zuizui0223/fcp/actions/runs/40080001765.
- Independent source-photo taxon output fidelity, fold-photo group stability and 999 genus-bootstrap slope-correlation audit: scripts/analysis/audit_fcp4475_observer_crossfit_photo_slope_20261024.py, passing source artifact re-read https://github.com/zuizui0223/fcp/actions/runs/40080546906.
- Immutable original derived source artifact ID 12554693263, SHA256 50d5d15a2d145b0db7a176f945435701251768ba614f659d9dac496011b33177; complete four-panel 17,900-row taxon CSV and full machine source JSON without original photo IDs, raw coordinates or observer IDs.
- Permanent full source JSON: results/fcp4475_observer_disjoint_thermal_white_verified_20261024.json.
- Independent audit JSON: results/fcp4475_crossfit_photo_slope_reliability_verified_20261024.json.

**Neither source photos nor parent PR #270 source results were changed.** The crossfit adds genuine conditional PHOTO evidence; it cannot replace independent field-verified fresh-petal phenotyping and genotype/photo matching if the scientific target is the maintenance of genetic flower-colour polymorphism.
