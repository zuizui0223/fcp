# PAL/WAL × FCP taxonomic alias audit — 2026-09-28

## Role

This audit resolves taxonomic naming differences only for the PAL/WAL × FCP overlap-feasibility analysis.

It does **not** alter the source registry, the Del Valle et al. (2019) PAL/WAL frequency summaries, or any frozen FCP species identity. Source strings remain preserved exactly as parsed from Supplementary Table S1. Resolved names are stored as an additional crosswalk solely to avoid false non-overlap caused by synonyms or obvious orthographic variants.

Machine-readable result:

- `results/polymorphism_pal_wal_h2_overlap_20260928/result.json`

## Resolutions relevant to the FCP cohorts

| Source string | Resolved FCP name | Resolution | Authority / note | FCP overlap |
|---|---|---|---|---|
| *Gymnadenia rhellicani* | *Gymnadenia rhellicani* | exact | source and FCP strings identical | reserve; PAL; D = 0.6515 |
| *Silene gallica* | *Silene gallica* | exact | source and FCP strings identical | reserve; PAL; D = 0.2344 |
| *Delphinium nelsonii* | *Delphinium nuttallianum* | accepted synonym | ITIS lists *D. nelsonii* as a synonym of accepted *D. nuttallianum*; POWO accepts *D. nuttallianum* | reserve; WAL; D = 0.1659 |
| *Mimulus guttatus* | *Erythranthe guttata* | accepted synonym | POWO lists *Mimulus guttatus* as a homotypic synonym of accepted *Erythranthe guttata* | reserve; WAL; D = 0.0000; source frequency is greenhouse-only |
| *Silene dioca* | *Silene dioica* | orthographic candidate resolved to accepted name | POWO accepts *Silene dioica*; the source string is retained unchanged and the one-character correction is used only for overlap auditing | reserve; WAL; D = 0.2055 |
| *Mimulis lewisii* | *Erythranthe lewisii* | source orthographic candidate + accepted synonym | POWO lists *Mimulus lewisii* as a homotypic synonym of accepted *Erythranthe lewisii*; the source string `Mimulis lewisii` is retained unchanged | no FCP overlap |

Authority pages used at audit time:

- POWO, *Silene dioica*: https://powo.science.kew.org/taxon/157239-1
- POWO, *Mimulus guttatus*: https://powo.science.kew.org/taxon/urn:lsid:ipni.org:names:162201-2
- POWO, *Erythranthe guttata*: https://powo.science.kew.org/taxon/urn:lsid:ipni.org:names:77120232-1
- POWO, *Mimulus lewisii*: https://powo.science.kew.org/taxon/317398-2
- POWO, *Erythranthe lewisii*: https://powo.science.kew.org/taxon/urn:lsid:ipni.org:names:77120063-1
- ITIS, *Delphinium nelsonii*: TSN 512147; accepted name *Delphinium nuttallianum*
- POWO, *Delphinium nuttallianum*: https://powo.science.kew.org/taxon/urn:lsid:ipni.org:names:77438-2

## Natural-frequency bridge subset

The Del Valle registry contains 26 rows representing 25 unique source species strings.

After taxonomic resolution, the natural-frequency bridge subset contains only four D-eligible FCP species, all in the reserve cohort:

- PAL: *Gymnadenia rhellicani* — D = 0.6515
- PAL: *Silene gallica* — D = 0.2344
- WAL: source *Delphinium nelsonii* → *D. nuttallianum* — D = 0.1659
- WAL: source *Silene dioca* → *S. dioica* — D = 0.2055

Descriptive medians are PAL = 0.443 and WAL = 0.186.

No inferential PAL-versus-WAL comparison is permitted: n = 2 versus 2, every overlap comes from the same reserve cohort, and there is no discovery or prospective overlap.

The *Mimulus guttatus* → *Erythranthe guttata* overlap is tracked separately because the source frequency entry is explicitly greenhouse-only and was already excluded from the natural-frequency PAL/WAL panel.

## Claim ceiling

This audit supports only:

> Taxonomic normalization increases the observed PAL/WAL–FCP overlap from two exact PAL matches to a four-species natural-frequency subset (two PAL, two WAL), but the overlap remains far too small and cohort-restricted for an inferential bridge test.

It does not support:

- a PAL-versus-WAL effect on D;
- a molecular validation of H2;
- a claim that the recurrent white axis is caused by anthocyanin loss;
- using the greenhouse-only *Mimulus guttatus* entry as natural persistence evidence;
- modifying the original source-table taxonomy or frequency analysis.
