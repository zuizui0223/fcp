#!/usr/bin/env python3
"""Apply the frozen acquisition-provenance wording to paper v0.1, fail closed."""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PATH = ROOT / "paper/polymorphism_v0_1/MANUSCRIPT.md"

REPLACEMENTS = [
    (
        "**Version:** v0.1-post-Step-9, 2026-09-10  ",
        "**Version:** v0.1-post-Step-9+acquisition-provenance, 2026-09-11  ",
    ),
    (
        "The analysis was built from georeferenced community photographs assembled under a fixed-photo design. Each candidate species entered image measurement with exactly 100 raw photographs. The fixed denominator was chosen to prevent common, heavily photographed species from dominating species-level comparisons simply because more images were available. Species inclusion in the final analytical frames was conditional on the upstream source, acquisition and measurement gates; the resulting species sets are therefore not random samples of all angiosperms and are not used to estimate angiosperm-wide prevalence.\n\nThe discovery polymorphism frame contained 369 species with at least 40 observations admitted to the four-state colour representation. A separate reserve campaign supplied 363 eligible species and was species-disjoint from discovery. The reserve was used to test transfer of the already defined polymorphism and spatial patterns rather than to search for a new threshold, colour axis or predictor set.",
        "The analysis was built from georeferenced community photographs assembled under a fixed-photo design. Both analytical tranches trace to the same frozen iNaturalist v1 observation candidate pool. The pre-pixel acquisition contract required `quality_grade=research`, photographs and georeferences, species-rank taxa, a frozen flowering-annotation gate (term 12, value 13), positional accuracy no worse than 5 km, unobscured coordinates and one of five allowed photo licences (`cc0`, `cc-by`, `cc-by-sa`, `cc-by-nc` or `cc-by-nc-sa`). Candidate-page and candidate-species selection did not use flower colour, and candidate image pixels remained unopened during acquisition. Within each species, an observer contributed at most two photographs and the final 100 photographs were chosen by deterministic geographic maximin selection.\n\nThe frozen candidate pool contained exactly 1,000 species with 100 raw photographs per species. The original fixed hash-ranked 500-taxon measurement budget defined the discovery tranche; the reserve was the entire species-disjoint complement of the other 500 taxa from the same frozen candidate pool. The reserve constructor replays the original discovery budget against the historical Git object and verifies the candidate-pool SHA-256 before taking that complement. Thus discovery and reserve differ in species membership but inherit the same frozen acquisition frame. Each candidate species therefore entered image measurement with exactly 100 raw photographs. The fixed denominator was chosen to prevent common, heavily photographed species from dominating species-level comparisons simply because more images were available. Species inclusion in the final analytical frames was conditional on the upstream source, acquisition and measurement gates; the resulting species sets are not random samples of all angiosperms and are not used to estimate angiosperm-wide prevalence.\n\nThe frozen acquisition query did **not** apply a native-range restriction or explicit `captive=false`/`wild=true` parameter, and the local metadata parser did not add such a filter. We therefore use Research Grade only as the iNaturalist quality-grade criterion and do not interpret it as proof of native-range-only or exclusively wild-population sampling. Geographic estimands in this paper describe organization among the observed community-photograph records.\n\nThe discovery polymorphism frame contained 369 species with at least 40 observations admitted to the four-state colour representation. A separate reserve campaign supplied 363 eligible species and was species-disjoint from discovery. The reserve was used to test transfer of the already defined polymorphism and spatial patterns rather than to search for a new threshold, colour axis or predictor set.",
    ),
    (
        "Four limitations define the present claim. First, the discovery and reserve frames are fixed-photo analytical samples, not random samples of angiosperms; the observed fractions of polymorphic species are therefore descriptive of admitted frames. Second, `D` summarizes four broad visible colour states and does not represent continuous reflectance spectra, pigment chemistry or pollinator-perceived colour. Third, ambiguous-palette observations remain biologically unresolved. `D_min4` and `D_max4` are exact only under the stated assumption that ambiguity can be completed into the four existing states, and applying one endpoint uniformly across all species is not an adversarial optimization over arbitrary species-specific latent allocations. Fourth, geographic organization is observational and does not identify selection, migration, drift or environmental causation.",
        "Five limitations define the present claim. First, the discovery and reserve frames are fixed-photo analytical samples, not random samples of angiosperms; the observed fractions of polymorphic species are therefore descriptive of admitted frames. Second, `D` summarizes four broad visible colour states and does not represent continuous reflectance spectra, pigment chemistry or pollinator-perceived colour. Third, ambiguous-palette observations remain biologically unresolved. `D_min4` and `D_max4` are exact only under the stated assumption that ambiguity can be completed into the four existing states, and applying one endpoint uniformly across all species is not an adversarial optimization over arbitrary species-specific latent allocations. Fourth, geographic organization is observational and does not identify selection, migration, drift or environmental causation. Fifth, although the common frozen acquisition frame required iNaturalist Research Grade, georeferencing, flowering annotation and the other stated metadata gates, it did not impose a native-range filter or explicit `captive=false`/`wild=true` parameter. Horticultural, introduced or otherwise non-native-range records may therefore remain, so the spatial result is an estimand of organization in the observed community-photograph sample rather than a native-range-only cline.",
    ),
    (
        "All analytical contracts, frozen input ledgers, result receipts, null summaries, synchronized paper-number ledgers, figure-data tables and figure-generation code are versioned in the project repository. The paper-facing numerical authority is `paper/polymorphism_v0_1/paper_numbers.json`; main figures are generated from the synchronized `paper/polymorphism_v0_1/figure_data/` tables. Repository/archive DOI and final public-data statement will be inserted at submission freeze.",
        "All analytical contracts, frozen input ledgers, result receipts, null summaries, synchronized paper-number ledgers, figure-data tables and figure-generation code are versioned in the project repository. Acquisition provenance is independently reconstructed by `scripts/analysis/audit_polymorphism_acquisition_provenance_20260911.py` and frozen in `results/polymorphism_acquisition_provenance_20260911/`; that audit opens neither candidate image pixels nor flower-colour outcomes. The paper-facing numerical authority is `paper/polymorphism_v0_1/paper_numbers.json`; main figures are generated from the synchronized `paper/polymorphism_v0_1/figure_data/` tables. Repository/archive DOI and final public-data statement will be inserted at submission freeze.",
    ),
]


def main() -> int:
    text = PATH.read_text(encoding="utf-8")
    original = text
    for old, new in REPLACEMENTS:
        count = text.count(old)
        if count != 1:
            raise RuntimeError(f"fail-closed manuscript patch: expected exactly one target, found {count}: {old[:80]!r}")
        text = text.replace(old, new, 1)
    if text == original:
        raise RuntimeError("manuscript patch made no change")
    PATH.write_text(text, encoding="utf-8", newline="\n")
    print(f"patched {PATH.relative_to(ROOT)} with {len(REPLACEMENTS)} frozen provenance replacements")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
