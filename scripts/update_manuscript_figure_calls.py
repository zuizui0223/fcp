#!/usr/bin/env python3
"""Idempotently insert canonical figure callouts/captions into JBI paper docs."""
from pathlib import Path

MANUSCRIPT = Path("docs/jbi_manuscript.md")
SI = Path("docs/jbi_supporting_information_index.md")


def insert_once(text: str, marker: str, addition: str, label: str) -> str:
    if addition.strip() in text:
        return text
    if marker not in text:
        raise RuntimeError(f"{label} insertion marker not found")
    return text.replace(marker, marker + addition, 1)


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
    text = insert_once(text, evidence_marker, evidence_add, "Evidence-base")

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
    text = insert_once(text, result_marker, result_add, "Primary-result")

    loo_marker = "Thus, no single represented family generated the shared negative direction.\n"
    loo_add = (
        "\nFigure 4 displays all 125 leave-one-family-out refits across the five climatic metrics. "
        "The figure makes the directional stability visible while retaining the important distinction "
        "that family deletion diagnoses concentration in individual families and is not a phylogenetic correction.\n"
    )
    text = insert_once(text, loo_marker, loo_add, "Leave-one-family-out")

    cr2_marker = (
        "The corresponding CR2 p-values were 0.6245 for temperature breadth, 0.3962 for climatic "
        "heterogeneity, 0.4550 for PCA dispersion and 0.1174 for PCA hull area.\n"
    )
    cr2_add = (
        "\nFigure 5 compares the primary family-clustered estimates with CR2/Satterthwaite, Open Tree/Grafen "
        "and dated-phylogeny treatments for all five metrics. Point estimates remain below one across treatments, "
        "whereas confidence intervals broaden under the phylogenetic and finite-cluster analyses.\n"
    )
    text = insert_once(text, cr2_marker, cr2_add, "Inference-sensitivity")

    power_marker = (
        "This diagnostic indicates that directional recovery under an effect of the observed magnitude can be "
        "substantially more stable than conventional significance in a 34-species design.\n"
    )
    power_add = (
        "\nSupporting Figure S2 expands this design diagnostic across the prespecified odds-ratio grid for all "
        "five metrics. It is presented only as a precision diagnostic and not as evidence for the ecological hypothesis.\n"
    )
    text = insert_once(text, power_marker, power_add, "Power-precision")

    if "## Figure captions" not in text:
        text += "\n\n## Figure captions\n"

    base_captions = {
        "**Figure 1. Geographic context of the 34 focal species.**": (
            "\n**Figure 1. Geographic context of the 34 focal species.** The world map shows the broader exact "
            "GBIF occurrence subset retained for the same 34 focal taxa, with symbols distinguishing species "
            "classified as within-population coexistence or geographic differentiation. The species strip reports "
            "the number of occupied climate cells used in the checksum-locked primary climatic analysis. The mapped "
            "occurrence archive is used for geographic context and auditability; it is not represented as the exact "
            "occurrence sample that created the frozen climatic metrics and is not morph-labelled.\n"
        ),
        "**Figure 2. Five climatic-niche metrics show the same direction of association with spatial organization.**": (
            "\n**Figure 2. Five climatic-niche metrics show the same direction of association with spatial organization.** "
            "Points are odds ratios from the production binomial models (`among ~ metric_z + effort_z`) and horizontal "
            "intervals are 95% Wald confidence intervals from family-clustered sandwich standard errors. Odds ratios "
            "below one indicate lower odds of geographically structured rather than within-population flower-colour "
            "variation as occupied climatic breadth increases. The figure is selected as the central result because it "
            "displays all five prespecified climatic summaries symmetrically rather than privileging the strongest "
            "unadjusted p-value.\n"
        ),
        "**Figure 3. Species-level climatic metrics underlying the comparative models.**": (
            "\n**Figure 3. Species-level climatic metrics underlying the comparative models.** Each point represents "
            "one of the 34 focal species (20 within-population, 14 geographically structured). Climatic metrics are "
            "standardized for visualization only so that all five dimensions can be displayed on a common scale; "
            "horizontal bars mark category medians. This is a visualization of the frozen observations, not an "
            "additional inferential analysis.\n"
        ),
    }
    for token, caption in base_captions.items():
        if token not in text:
            text += caption

    figure4_caption = (
        "\n**Figure 4. Leave-one-family-out stability across all five climatic-niche metrics.** Grey points are odds "
        "ratios from 25 unclustered refits per metric, each omitting one represented plant family; diamonds show the "
        "full 34-species point estimates and horizontal grey segments span the deletion range. All 125 deletion "
        "estimates remain below one. Family deletion tests whether the shared direction is concentrated in an "
        "individual represented family; it does not establish phylogenetic independence.\n"
    )
    figure5_caption = (
        "\n**Figure 5. Effect direction and uncertainty across inferential treatments.** For each climatic metric, "
        "the primary family-clustered estimate is shown alongside CR2/Satterthwaite finite-cluster inference, the "
        "median of 100 Open Tree/Grafen polytomy resolutions, and the median dated-phylogeny estimate across "
        "V.PhyloMaker2 scenarios S1-S3. Horizontal intervals are the corresponding 95% intervals; for the dated "
        "phylogeny the plotted interval is the envelope across the three scenario-specific intervals. All point "
        "estimates remain below one, while stricter treatments broaden uncertainty. These treatments are sensitivity "
        "analyses of the same data, not independent tests.\n"
    )
    s1_caption = (
        "\n**Supporting Figure S1. Geographic occurrence context for each of the 34 focal species.** Species-specific "
        "maps show the broader exact GBIF occurrence subset retained for citation and distribution-context quality "
        "control, with panel labels indicating the spatial-organization class and occurrence count. These maps provide "
        "an auditable view of sampled geographic ranges but do not represent morph-specific ranges or the exact primary "
        "occurrence sample used to construct the frozen climatic metrics.\n"
    )
    s2_caption = (
        "\n**Supporting Figure S2. Finite-sample design diagnostic across specified effect sizes.** The two panels show, "
        "for each climatic metric and simulated true odds ratio, the probability that the fitted coefficient is negative "
        "and the probability that the family-clustered Wald test yields p < 0.05 across 3,000 simulations. Simulations "
        "retain the observed 34-species predictor, effort and family structure. The figure describes expected precision "
        "under specified effects and is not evidence for the ecological hypothesis or a post-hoc adequacy criterion.\n"
    )

    support_token = "**Supporting Figure S1. Geographic occurrence context for each of the 34 focal species.**"
    if "**Figure 4. Leave-one-family-out stability across all five climatic-niche metrics.**" not in text:
        if support_token in text:
            text = text.replace("\n" + support_token, figure4_caption + "\n" + support_token, 1)
        else:
            text += figure4_caption
    if "**Figure 5. Effect direction and uncertainty across inferential treatments.**" not in text:
        if support_token in text:
            text = text.replace("\n" + support_token, figure5_caption + "\n" + support_token, 1)
        else:
            text += figure5_caption
    if support_token not in text:
        text += s1_caption
    if "**Supporting Figure S2. Finite-sample design diagnostic across specified effect sizes.**" not in text:
        text += s2_caption

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
- `docs/figures/figure4_leave_one_family_out.png` / `.pdf`
  - 25 leave-one-family-out refits for each of the five metrics; tests concentration in individual families.
- `docs/figures/figure5_inference_method_sensitivity.png` / `.pdf`
  - primary, CR2/Satterthwaite, Open Tree/Grafen and dated-phylogeny estimates compared across all five metrics.
- `docs/figures/figureS1_34_species_distribution_context.png` / `.pdf`
  - species-specific geographic occurrence maps using the broader exact GBIF subset as supporting distribution context/QC.
- `docs/figures/figureS2_power_precision_design.png` / `.pdf`
  - design-based sign-recovery and p < 0.05 probabilities across specified effect sizes; diagnostic only.
- `scripts/make_paper_figures.py`
  - canonical figure generator.
- `docs/FIGURE_PLAN.md`
  - explains why each figure is assigned to the main text or Supporting Information.

Figure 1 and Supporting Figure S1 must retain the distinction between the broader exact GBIF citation subset and the exact occurrence sampling path that produced the frozen primary climatic summaries.
"""
        text = text.replace(marker, section + marker, 1)
    else:
        additions = [
            (
                "figure4_leave_one_family_out",
                "- `docs/figures/figure4_leave_one_family_out.png` / `.pdf`\n"
                "  - 25 leave-one-family-out refits for each of the five metrics; tests concentration in individual families.\n"
            ),
            (
                "figure5_inference_method_sensitivity",
                "- `docs/figures/figure5_inference_method_sensitivity.png` / `.pdf`\n"
                "  - primary, CR2/Satterthwaite, Open Tree/Grafen and dated-phylogeny estimates compared across all five metrics.\n"
            ),
            (
                "figureS2_power_precision_design",
                "- `docs/figures/figureS2_power_precision_design.png` / `.pdf`\n"
                "  - design-based sign-recovery and p < 0.05 probabilities across specified effect sizes; diagnostic only.\n"
            ),
        ]
        generator_marker = "- `scripts/make_paper_figures.py`\n"
        for token, addition in additions:
            if token not in text:
                if generator_marker not in text:
                    raise RuntimeError("Canonical-figure generator marker not found")
                text = text.replace(generator_marker, addition + generator_marker, 1)
    SI.write_text(text, encoding="utf-8")


if __name__ == "__main__":
    update_manuscript()
    update_si()
    print("Updated manuscript figure callouts/captions and Supporting Information index")
