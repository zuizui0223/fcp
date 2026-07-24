#!/usr/bin/env python3
from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
MANUSCRIPT = ROOT / "docs/jbi_manuscript_editorial_revision.md"
SUMMARY = ROOT / "docs/supporting/jbi_table_s16_dated_phylogeny_sensitivity_summary.csv"


def replace_section(text: str, start: str, end: str, replacement: str) -> str:
    pattern = re.compile(
        rf"(?ms)^{re.escape(start)}\n.*?(?=^{re.escape(end)}\n)"
    )
    matches = pattern.findall(text)
    if len(matches) != 1:
        raise AssertionError(f"Expected one section {start!r} -> {end!r}; found {len(matches)}")
    return pattern.sub(replacement.rstrip() + "\n\n", text, count=1)


def replace_once(text: str, old: str, new: str) -> str:
    count = text.count(old)
    if count != 1:
        raise AssertionError(f"Expected one occurrence; found {count}: {old[:100]!r}")
    return text.replace(old, new, 1)


def validate_summary() -> None:
    x = pd.read_csv(SUMMARY)
    if len(x) != 6 or set(x["scenario"]) != {"S1", "S2", "S3"}:
        raise AssertionError("Dated phylogeny summary must contain six S1-S3 rows")
    if not (x["n_species"] == 34).all() or not (x["status"] == "complete").all():
        raise AssertionError("All dated models must contain 34 species and be complete")
    primary = x[x["analysis_label"] == "baseline_primary_occurrences"]
    paginated = x[x["analysis_label"] == "baseline_paginated_quality_filtered"]
    expected = {
        "primary_or": (0.4639, 0.4705),
        "primary_p": (0.1213, 0.1244),
        "paginated_or": (0.3658, 0.3693),
        "paginated_p": (0.0676, 0.0691),
    }
    checks = {
        "primary_or": (primary["odds_ratio"].min(), primary["odds_ratio"].max()),
        "primary_p": (primary["p_value"].min(), primary["p_value"].max()),
        "paginated_or": (paginated["odds_ratio"].min(), paginated["odds_ratio"].max()),
        "paginated_p": (paginated["p_value"].min(), paginated["p_value"].max()),
    }
    for key, (lo, hi) in checks.items():
        elo, ehi = expected[key]
        if lo < elo or hi > ehi:
            raise AssertionError(f"Unexpected {key} range: {(lo, hi)}")


def main() -> None:
    validate_summary()
    text = MANUSCRIPT.read_text(encoding="utf-8")

    abstract = """## Abstract
### Aim

Intraspecific phenotypic variation may occur through local coexistence or geographic differentiation, but comparative studies rarely treat spatial organization as a distinct property. We tested whether documented within-population flower-colour polymorphism and geographically structured flower-colour variation differ in occupied climatic niche breadth.

### Location

Global, literature-derived sample.

### Taxon

Angiosperms with documented natural intraspecific flower-colour variation.

### Methods

We combined audited literature classifications with GBIF occurrences and WorldClim 2.1 data. Binomial models used family-clustered uncertainty, label permutations and family deletion. We repeated the focal model with a paginated, quality-filtered GBIF sample and fitted phylogenetic logistic models using both an Open Tree topology and a time-scaled vascular-plant megaphylogeny.

### Results

The baseline model included 34 species from 25 families. The primary association was negative (odds ratio = 0.426, 95% confidence interval = 0.184–0.985; Wald p = 0.0460), although permutation support was borderline (p = 0.0556). Paginated sampling strengthened the estimate (odds ratio = 0.300, 95% confidence interval = 0.133–0.675; permutation p = 0.0164). Open Tree models retained the negative direction but were imprecise. Time-scaled megaphylogeny models retained all 34 species and were nearly invariant among three placement scenarios: odds ratios were 0.464–0.470 for the primary data and 0.366–0.369 for the paginated data; all phylogenetic confidence intervals included one.

### Main conclusions

The association was stable to stronger occurrence sampling and remained negative under two phylogenetic treatments, but phylogenetic uncertainty intervals included one. The evidence is exploratory and does not test morph-specific tolerance, local adaptation or climatic causation.
"""
    text = replace_section(text, "## Abstract", "**Keywords:**", abstract)

    methods = """### Occurrence-sampling and phylogenetic sensitivity models

The paginated dataset was analysed with the same binomial model, 9,999 label permutations and family-deletion procedure as the primary baseline dataset (Tables S9–S12). For topology-based phylogenetic sensitivity, names were matched to Open Tree Taxonomy without approximate matching; unique score-one matches represented in the synthetic tree were retained. Thirty of 34 species were eligible. We induced the Open Tree topology (Hinchliff et al., 2015), randomly resolved polytomies 100 times, assigned Grafen branch lengths and fitted `phyloglm` logistic MPLE models (Ho & Ané, 2014) with standardised moisture breadth and occurrence effort (Tables S13–S15).

We then fitted a dated-megaphylogeny sensitivity analysis using V.PhyloMaker2 (Jin & Qian, 2019, 2022) and the time-scaled `GBOTB.extended.LCVP` backbone, which incorporates the broad seed-plant phylogeny of Smith and Brown (2018) and earlier time-calibrated plant-tree information (Zanne et al., 2014). Twenty-eight focal species were already represented in the backbone and six were bound to it; no species failed to bind. Using the fixed seed 20260724, we generated placement scenarios S1, S2 and S3 and fitted the same `phyloglm` logistic MPLE formula to each tree for both occurrence datasets. All trees retained 34 species, had finite positive branch lengths and required no branch-length replacement. These models inherit time scaling from the backbone, but the six inserted taxa remain subject to V.PhyloMaker2 placement assumptions (Tables S16–S17).
"""
    text = replace_section(
        text,
        "### Occurrence-sampling and phylogenetic sensitivity models",
        "### Candidate-versus-control comparison",
        methods,
    )

    results = """### Occurrence-sampling sensitivity

The paginated moisture model strengthened to an odds ratio of 0.300 (family-clustered 95% confidence interval = 0.133–0.675; Wald p = 0.00361; permutation p = 0.0164; Table 4; Table S11). Leave-one-family-out odds ratios ranged from 0.247 to 0.359 and remained below one (Table S12). This result reduces concern that the primary direction was generated by the 300-record cap, although the sensitivity sample was itself capped.

### Topology-based phylogenetic sensitivity

Open Tree matching retained 30 species, and all 100 fits completed for each occurrence dataset. The primary-data model gave an odds ratio of 0.592 (95% confidence interval = 0.244–1.434; p = 0.246); the paginated-data model gave 0.472 (0.175–1.272; p = 0.138). Every replicate was negative, but both intervals included one (Table 4; Tables S13–S15).

### Dated-megaphylogeny sensitivity

V.PhyloMaker2 retained all 34 species: 28 were already present in the time-scaled backbone and six were inserted. Results were nearly invariant among scenarios S1–S3. For the primary occurrence data, odds ratios ranged from 0.464 to 0.470, confidence-interval lower limits from 0.176 to 0.180, upper limits from 1.226 to 1.231 and p-values from 0.121 to 0.124. For the paginated data, odds ratios ranged from 0.366 to 0.369, lower limits from 0.124 to 0.127, upper limits from 1.075 to 1.081 and p-values from 0.0677 to 0.0691 (Table 4; Tables S16–S17). Thus, time-scaled branch lengths and retention of all species did not reverse the negative direction, but every interval still included one.
"""
    text = replace_section(
        text,
        "### Occurrence-sampling sensitivity",
        "### Coarse occurrence-cloud alternatives",
        results,
    )

    old_discussion = """The topology-based analysis changed the evidential balance. Negative coefficients persisted under both occurrence designs, but confidence intervals were broad. This may reflect shared evolutionary structure, reduced sample size or uncertainty from a synthetic topology and assumed Grafen branch lengths. The analysis is stronger than family clustering alone but is not equivalent to a dated species phylogeny; it does not establish that ancestry either explains or fails to explain the association.
"""
    new_discussion = """The two phylogenetic analyses changed the evidential balance without reversing the estimated direction. Open Tree models used 30 matched species and broad confidence intervals, potentially reflecting reduced sample size, shared evolutionary structure and uncertainty from a synthetic topology with Grafen branch lengths. The dated-megaphylogeny models retained all 34 species and produced nearly identical estimates under placement scenarios S1–S3. Their paginated-data intervals approached but did not exclude an odds ratio of one. Six species were inserted rather than directly represented in the backbone, so this agreement does not remove phylogenetic uncertainty. Together, the topology-based and time-scaled analyses show that phylogenetic treatment attenuates certainty more than direction; neither supports describing the association as phylogenetically confirmed.
"""
    text = replace_once(text, old_discussion, new_discussion)

    old_final = """Despite these constraints, the study establishes a reproducible framework for comparing the spatial organization of intraspecific phenotypic variation. The audited evidence suggests that documented geographic differentiation may be associated with narrower occupied moisture breadth than documented local coexistence, while simultaneously showing that the strength of this inference depends on classification quality, occurrence sampling and phylogenetic treatment. Morph-labelled locality data and a dated species phylogeny are needed to determine whether the association generalizes and which ecological or historical processes generate it.
"""
    new_final = """Despite these constraints, the study establishes a reproducible framework for comparing the spatial organization of intraspecific phenotypic variation. The audited evidence suggests that documented geographic differentiation may be associated with narrower occupied moisture breadth than documented local coexistence, while simultaneously showing that the strength of this inference depends on classification quality, occurrence sampling and phylogenetic treatment. Morph-labelled locality data and a species-specific molecular phylogeny would refine the estimate and help identify the ecological or historical processes that generate the pattern.
"""
    text = replace_once(text, old_final, new_final)

    reference_anchor = """Ho, L. S. T., & Ané, C. (2014). A linear-time algorithm for Gaussian and non-Gaussian trait evolution models. *Systematic Biology, 63*(3), 397–408. https://doi.org/10.1093/sysbio/syu005
"""
    reference_insert = reference_anchor + """
Jin, Y., & Qian, H. (2019). V.PhyloMaker: An R package that can generate very large phylogenies for vascular plants. *Ecography, 42*(8), 1353–1359. https://doi.org/10.1111/ecog.04434

Jin, Y., & Qian, H. (2022). V.PhyloMaker2: An updated and enlarged R package that can generate very large phylogenies for vascular plants. *Plant Diversity, 44*(4), 335–339. https://doi.org/10.1016/j.pld.2022.05.005
"""
    text = replace_once(text, reference_anchor, reference_insert)

    smith_anchor = """Seabold, S., & Perktold, J. (2010). Statsmodels: Econometric and statistical modeling with Python. In *Proceedings of the 9th Python in Science Conference* (pp. 92–96).
"""
    smith_insert = smith_anchor + """
Smith, S. A., & Brown, J. W. (2018). Constructing a broadly inclusive seed plant phylogeny. *American Journal of Botany, 105*(3), 302–314. https://doi.org/10.1002/ajb2.1019
"""
    text = replace_once(text, smith_anchor, smith_insert)

    zanne_anchor = """Wessinger, C. A., & Rausher, M. D. (2012). Lessons from flower colour evolution on targets of selection. *Journal of Experimental Botany, 63*(16), 5741–5749. https://doi.org/10.1093/jxb/ers267
"""
    zanne_insert = zanne_anchor + """
Zanne, A. E., Tank, D. C., Cornwell, W. K., Eastman, J. M., Smith, S. A., FitzJohn, R. G., McGlinn, D. J., O'Meara, B. C., Moles, A. T., Reich, P. B., Royer, D. L., Soltis, D. E., Stevens, P. F., Westoby, M., Wright, I. J., Aarssen, L., Bertin, R. I., Calaminus, A., Govaerts, R., . . . Beaulieu, J. M. (2014). Three keys to the radiation of angiosperms into freezing environments. *Nature, 506*(7486), 89–92. https://doi.org/10.1038/nature12872
"""
    text = replace_once(text, zanne_anchor, zanne_insert)

    table4 = """### Table 4. Occurrence-sampling and phylogenetic sensitivity

| Occurrence dataset | Model | Species | Odds ratio | 95% CI | Wald p | Permutation p |
|---|---|---:|---:|---:|---:|---:|
| Primary | Family-clustered GLM | 34 | 0.426 | 0.184–0.985 | 0.0460 | 0.0556 |
| Paginated quality-filtered | Family-clustered GLM | 34 | 0.300 | 0.133–0.675 | 0.00361 | 0.0164 |
| Primary | Open Tree topology-based phylogenetic logistic | 30 | 0.592 | 0.244–1.434 | 0.246 | — |
| Paginated quality-filtered | Open Tree topology-based phylogenetic logistic | 30 | 0.472 | 0.175–1.272 | 0.138 | — |
| Primary | V.PhyloMaker2 dated megaphylogeny, S1–S3 | 34 | 0.464–0.470 | 0.176–1.231 | 0.121–0.124 | — |
| Paginated quality-filtered | V.PhyloMaker2 dated megaphylogeny, S1–S3 | 34 | 0.366–0.369 | 0.124–1.081 | 0.0677–0.0691 | — |

*Note.* Open Tree values summarize 100 completed polytomy resolutions of an induced topology with Grafen branch lengths. Dated-megaphylogeny values show the range across V.PhyloMaker2 placement scenarios S1–S3; the confidence interval is the envelope across scenario-specific Wald intervals. Six species were inserted into the time-scaled backbone.
"""
    text = replace_section(text, "### Table 4. Occurrence-sampling and topology-based phylogenetic sensitivity", "## Figure legends and embedded figures", table4)

    supporting = """## Supporting Information
Tables S1–S7 contain the original model matrices and audits. **Table S8** is the primary baseline model dataset; **Table S9** is the paginated model dataset; **Table S10** is the GBIF taxonomic and retention audit; **Tables S11–S12** report paginated robustness and family deletion; **Tables S13–S15** report Open Tree phylogenetic summaries, name resolution and all 200 replicate fits; **Table S16** reports the six fixed-seed dated-megaphylogeny models; and **Table S17** records V.PhyloMaker2 species placement. The Supporting Information index also identifies the GBIF QC manifests, Open Tree topology, dated S1–S3 trees and both phylogenetic provenance manifests.
"""
    text = re.sub(r"(?ms)^## Supporting Information\n.*\Z", supporting.rstrip() + "\n", text, count=1)

    abstract_text = text[text.index("## Abstract"):text.index("**Keywords:")]
    abstract_text = re.sub(r"^### .*?$", "", abstract_text, flags=re.MULTILINE)
    words = re.findall(r"\b[\w’'-]+\b", abstract_text)
    if len(words) > 300:
        raise AssertionError(f"Abstract exceeds 300 words: {len(words)}")
    required = [
        "odds ratios were 0.464–0.470",
        "odds ratios ranged from 0.366 to 0.369",
        "Tables S16–S17",
        "V.PhyloMaker2 dated megaphylogeny, S1–S3",
        "species-specific molecular phylogeny would refine",
    ]
    for phrase in required:
        if phrase not in text:
            raise AssertionError(f"Required manuscript phrase missing: {phrase}")

    MANUSCRIPT.write_text(text, encoding="utf-8")
    print({"abstract_words": len(words), "manuscript": str(MANUSCRIPT)})


if __name__ == "__main__":
    main()
