#!/usr/bin/env python3
"""Idempotently insert canonical figure callouts/captions into JBI paper docs."""
from pathlib import Path

MANUSCRIPT = Path("docs/jbi_manuscript.md")
SI = Path("docs/jbi_supporting_information_index.md")


def update_manuscript() -> None:
    text = MANUSCRIPT.read_text(encoding="utf-8")

    evidence_marker = (
        "The retained evidence path contained 664 candidate species from 140 families and a resolved "
        "111-species review queue. After the binary spatial-state and climatic-eligibility requirements "
        "were applied, the final comparison contained 34 species from 25 families: 20 within-population "
        "and 14 geographically structured cases. Because discovery effort was literature-dependent and "
        "the candidates were not a random sample of angiosperms, these counts should not be interpreted "
        "as prevalence estimates.\n"
    )
    evidence_add = (
        "\nThe geographic context of the 34 focal species is shown in Figure 1, and species-specific "
        "occurrence maps are provided as Supporting Figure S1. Those maps use the broader exact GBIF "
        "occurrence subset retained for citation and distribution-context auditing; they are not presented "
        "as the exact occurrence sample that generated the checksum-locked climatic metrics, nor are the "
        "records labelled by flower-colour morph.\n"
    )
    if evidence_add.strip() not in text:
        if evidence_marker not in text:
            raise RuntimeError("Evidence-base insertion marker not found")
        text = text.replace(evidence_marker, evidence_marker + evidence_add, 1)

    result_marker = (
        "After Holm adjustment across five metrics, the moisture-breadth clustered Wald p-value was 0.184 "
        "and the permutation p-value was 0.212. No metric retained conventional statistical support after "
        "this correction. Moisture breadth is therefore the strongest observed association within a broader "
        "directional pattern, not a uniquely established climatic driver.\n"
    )
    result_add = (
        "\nThe common direction of the five effect estimates is summarized in Figure 2. Figure 3 shows the "
        "34 species themselves across the five standardized climatic metrics, so the overlap, extreme "
        "observations and modest sample size remain visible rather than being represented only by model "
        "coefficients.\n"
    )
    if result_add.strip() not in text:
        if result_marker not in text:
            raise RuntimeError("Primary-result insertion marker not found")
        text = text.replace(result_marker, result_marker + result_add, 1)

    if "## Figure captions" not in text:
        text += """

## Figure captions

**Figure 1. Geographic context of the 34 focal species.** The world map shows the broader exact GBIF occurrence subset retained for the same 34 focal taxa, with symbols distinguishing species classified as within-population coexistence or geographic differentiation. The species strip reports the number of occupied climate cells used in the checksum-locked primary climatic analysis. The mapped occurrence archive is used for geographic context and auditability; it is not represented as the exact occurrence sample that created the frozen climatic metrics and is not morph-labelled.

**Figure 2. Five climatic-niche metrics show the same direction of association with spatial organization.** Points are odds ratios from the production binomial models (`among ~ metric_z + effort_z`) and horizontal intervals are 95% Wald confidence intervals from family-clustered sandwich standard errors. Odds ratios below one indicate lower odds of geographically structured rather than within-population flower-colour variation as occupied climatic breadth increases. The figure is selected as the central result because it displays all five prespecified climatic summaries symmetrically rather than privileging the strongest unadjusted p-value.

**Figure 3. Species-level climatic metrics underlying the comparative models.** Each point represents one of the 34 focal species (20 within-population, 14 geographically structured). Climatic metrics are standardized for visualization only so that all five dimensions can be displayed on a common scale; horizontal bars mark category medians. This is a visualization of the frozen observations, not an additional inferential analysis.

**Supporting Figure S1. Geographic occurrence context for each of the 34 focal species.** Species-specific maps show the broader exact GBIF occurrence subset retained for citation and distribution-context quality control, with panel labels indicating the spatial-organization class and occurrence count. These maps provide an auditable view of sampled geographic ranges but do not represent morph-specific ranges or the exact primary occurrence sample used to construct the frozen climatic metrics.
"""

    MANUSCRIPT.write_text(text, encoding="utf-8")


def update_si() -> None:
    text = SI.read_text(encoding="utf-8")
    if "## 6. Canonical figures" not in text:
        marker = "\n## Submission boundary\n"
        if marker not in text:
            raise RuntimeError("Supporting-index insertion marker not found")
        section = """

## 6. Canonical figures

- `docs/figures/figure1_geographic_context.png` / `.pdf`
  - main-text geographic context for the 34 focal species; broader exact GBIF occurrence subset is shown only as context/QC.
- `docs/figures/figure2_five_metric_forest.png` / `.pdf`
  - central five-metric family-clustered effect-size figure generated directly from the frozen 34-species dataset.
- `docs/figures/figure3_raw_species_metrics.png` / `.pdf`
  - all 34 species shown across standardized versions of the five frozen climatic metrics; visualization only.
- `docs/figures/figureS1_34_species_distribution_context.png` / `.pdf`
  - species-specific geographic occurrence maps using the broader exact GBIF subset as supporting distribution context/QC.
- `scripts/make_paper_figures.py`
  - canonical figure generator.
- `docs/FIGURE_PLAN.md`
  - explains why these figures were selected from the pipeline and why robustness diagnostics remain supporting material.

Figure 1 and Supporting Figure S1 must retain the distinction between the broader exact GBIF citation subset and the exact occurrence sampling path that produced the frozen primary climatic summaries.
"""
        text = text.replace(marker, section + marker, 1)
    SI.write_text(text, encoding="utf-8")


if __name__ == "__main__":
    update_manuscript()
    update_si()
    print("Updated manuscript figure callouts/captions and Supporting Information index")
