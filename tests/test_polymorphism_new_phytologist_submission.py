import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MANUSCRIPT = ROOT / "docs" / "POLYMORPHISM_MANUSCRIPT_NEW_PHYTOLOGIST.md"
CANONICAL = ROOT / "docs" / "POLYMORPHISM_MANUSCRIPT.md"
COVER = ROOT / "docs" / "POLYMORPHISM_NEW_PHYTOLOGIST_COVER_LETTER.md"
RGFCA_INTERPRETATION = ROOT / "docs" / "RGFCA_TO_POLYMORPHISM_INTERPRETATION_20260918.md"
FRAME_PROVENANCE = ROOT / "docs" / "POLYMORPHISM_42111_FRAME_PROVENANCE_20260918.md"
SUPPORTING = ROOT / "docs" / "POLYMORPHISM_SUPPORTING_INFORMATION_20260918.md"
VALIDITY = ROOT / "results" / "polymorphism_h2_posthoc_validity_diagnostics_20260922" / "result.json"
HIGHLIGHT = ROOT / "results" / "polymorphism_h2_third_cohort_highlight_validity_20260922" / "result.json"
ADJUDICATION = ROOT / "docs" / "POLYMORPHISM_H2_THIRD_COHORT_HIGHLIGHT_DECISION_ADJUDICATION_20260923.md"
D_TRANSPORT = ROOT / "results" / "polymorphism_fresh_D_transport_20260925" / "result.json"
D_FINITE = ROOT / "results" / "polymorphism_D_finite_sample_sensitivity_20260928" / "result.json"
WHITE_ENV = ROOT / "results" / "polymorphism_white_environment_mechanism_20260925" / "result.json"
WHITE_ENV_OBSERVER = ROOT / "results" / "polymorphism_white_environment_observer_sensitivity_20260925" / "result.json"
BIO5_TRANSPORT = ROOT / "results" / "polymorphism_legacy_white_bio5_replication_20260925" / "result.json"
SILENE_MOLECULAR = ROOT / "results" / "polymorphism_silene_molecular_anchor_20260928" / "result.json"
PAL_WAL_OVERLAP = ROOT / "results" / "polymorphism_pal_wal_h2_overlap_20260928" / "result.json"
LINEAGE_MAP = ROOT / "docs" / "POLYMORPHISM_DATA_LINEAGE_MAP_20260925.md"
WHITE_ENV_PROTOCOL = ROOT / "docs" / "POLYMORPHISM_WHITE_ENVIRONMENT_MECHANISM_PROTOCOL_20260925.md"
WHITE_ENV_SCRIPT = ROOT / "scripts" / "analysis" / "run_white_environment_mechanism_20260925.py"
BIO5_PROTOCOL = ROOT / "docs" / "POLYMORPHISM_LEGACY_WHITE_BIO5_REPLICATION_PROTOCOL_20260925.md"
BIO5_SCRIPT = ROOT / "scripts" / "analysis" / "run_legacy_white_bio5_replication_20260925.py"
H3B_PROTOCOL = ROOT / "docs" / "POLYMORPHISM_H3B_RESERVE_SPAN_PROTOCOL_20260912.md"
H3B_SCRIPT = ROOT / "scripts" / "analysis" / "run_polymorphism_h3b_reserve_span_20260912.R"
H3B_RESULT = ROOT / "results" / "polymorphism_h3b_reserve_span_20260912" / "result.json"
ARCHIVE_ROOT = ROOT / "archive" / "fcp_submission_20260925"
WORLDCLIM_CHECKSUMS = ARCHIVE_ROOT / "worldclim_checksums.txt"
WORLDCLIM_RELEASE_RECEIPT = ARCHIVE_ROOT / "WORLDCLIM_RELEASE_RECEIPT.md"
NP_PROVENANCE_RELEASE_RECEIPT = ARCHIVE_ROOT / "NP_PROVENANCE_RELEASE_RECEIPT.md"


def words(text: str) -> list[str]:
    return re.findall(r"[A-Za-z0-9]+(?:[’'-][A-Za-z0-9]+)*", text)


def test_new_phytologist_front_matter_and_summary_contract() -> None:
    text = MANUSCRIPT.read_text(encoding="utf-8")
    title = text.splitlines()[0].removeprefix("# ").strip()
    assert len(title) <= 130

    summary = text.split("## Summary", 1)[1].split("\n\n---", 1)[0]
    bullets = [line for line in summary.splitlines() if line.startswith("- ")]
    assert len(bullets) == 4
    assert len(words(summary)) <= 200

    keyword_line = next(
        line for line in text.splitlines()
        if line.startswith("**Keywords (alphabetical):**")
    )
    keywords = [x.strip() for x in keyword_line.split("**", 2)[-1].split(":", 1)[-1].split(";") if x.strip()]
    assert 5 <= len(keywords) <= 8


def test_new_phytologist_explains_current_data_architecture_and_necessity() -> None:
    text = MANUSCRIPT.read_text(encoding="utf-8")
    for token in (
        "Study design, data provenance and why each stage was required",
        "Stage 0 — global outcome-blind opportunity frame",
        "Stage 1 — discovery and species-disjoint validation image cohorts",
        "Stage 2 — prospective species-disjoint confirmation cohort",
        "Stage 3 — targeted technical and ecological annotations",
        "A new species- and photo-disjoint cohort was therefore necessary",
        "These are alternative-explanation and mechanism filters",
    ):
        assert token in text

def test_itv_framing_and_reader_chronology_are_explicit() -> None:
    text = MANUSCRIPT.read_text(encoding="utf-8")
    supporting = SUPPORTING.read_text(encoding="utf-8")

    for token in (
        "Species means can erase the structure of intraspecific trait variation",
        "how much",
        "which phenotypic directions",
        "how the variants are arranged geographically",
        "Flower-colour polymorphism is unusually suited to this distributional view",
        "each stage addresses a different inferential failure mode",
        "Phenotype-space generality without a universal geographic map",
        "treating species as distributions rather than mean phenotypes",
    ):
        assert token in text

    assert "An upstream shared-boundary analysis did not yield" not in text

    order = (
        "## S1. Data sources, selection and inferential necessity",
        "## S2. H1 observer-disjoint measurement validity",
        "## S3. Discovery/validation H2 target localization",
        "## S4. Prospective-confirmation H2 chain of custody",
        "## S5. Replicated D–spatial organization",
        "## S6. H3a broad phylogenetic-signal boundary",
        "## S7. H3b sampled-span replication boundary",
        "## S8. Bounded secondary mechanism evidence",
        "## S9. Canonical main-text figures",
    )
    positions = [supporting.index(token) for token in order]
    assert positions == sorted(positions)
    assert "sampling frame -> measurement validity -> geometry discovery -> prospective confirmation" in supporting


def test_frame_provenance_preserves_sampling_boundary() -> None:
    assert FRAME_PROVENANCE.exists()
    frame = FRAME_PROVENANCE.read_text(encoding="utf-8")
    for token in (
        "42,111 species",
        "U100 = 4,730 species",
        "No candidate image pixels or flower-colour outcomes were used",
        "not the denominator for estimating global polymorphism prevalence",
    ):
        assert token in frame

def test_new_phytologist_documents_42111_frame_provenance() -> None:
    text = MANUSCRIPT.read_text(encoding="utf-8")
    frame = FRAME_PROVENANCE.read_text(encoding="utf-8")

    # Main text keeps the inferential frame concise.
    for token in (
        "42,111 unique iNaturalist plant species",
        "4,730 species",
        "POLYMORPHISM_42111_FRAME_PROVENANCE_20260918.md",
    ):
        assert token in text

    # Exact discovery implementation belongs to the provenance record.
    for token in (
        "18 × 9 equal-area global grid",
        "20 metadata-only rounds",
        "3,240 request attempts",
        "No candidate image pixels or flower-colour outcomes were used",
    ):
        assert token in frame


def test_new_phytologist_documents_acquisition_lineage_and_stage_roles() -> None:
    text = MANUSCRIPT.read_text(encoding="utf-8")
    supporting = SUPPORTING.read_text(encoding="utf-8")

    for token in (
        "Research Grade species-rank iNaturalist observations",
        "Selection was colour-blind",
        "deterministic geographic maximin sampling",
        "Table 1. Data sources, lineage and inferential necessity",
        "Discovery image cohort",
        "Species-disjoint validation image cohort",
        "Prospective confirmation cohort",
        "Why it was required",
    ):
        assert token in text

    for token in (
        "Research Grade species-rank iNaturalist records",
        "positional accuracy <=5 km",
        "Observer contribution was capped at two photographs per species",
        "deterministic geographic maximin sampling fixed 100 raw photographs per species",
        "No native-range restriction or explicit captive/wild filter was imposed",
    ):
        assert token in supporting

def test_new_phytologist_required_sections_and_display_items() -> None:
    text = MANUSCRIPT.read_text(encoding="utf-8")
    for heading in (
        "## Summary",
        "## Introduction",
        "## Materials and Methods",
        "## Table 1.",
        "## Results",
        "## Discussion",
        "## Acknowledgements",
        "## Competing interests",
        "## Author contributions",
        "## Data availability",
        "## Supporting Information",
        "## References",
    ):
        assert heading in text

    assert "- Figures: 5" in text
    assert "- Tables: 1" in text



def test_new_phytologist_imports_only_bounded_fresh_D_transport() -> None:
    import json
    import pytest

    text = MANUSCRIPT.read_text(encoding="utf-8")
    canonical = CANONICAL.read_text(encoding="utf-8")
    result = json.loads(D_TRANSPORT.read_text(encoding="utf-8"))

    assert result["overlap_species"] == 136
    assert result["spearman_rho"] == pytest.approx(0.9681064161858759)
    assert result["lin_ccc"] == pytest.approx(0.9720478011581966)
    assert result["calibration_slope"] == pytest.approx(0.96909589220185)
    assert result["absolute_D_change"]["median"] == pytest.approx(0.014761943866009125)

    for manuscript in (text, canonical):
        for token in (
            "136",
            "0.968",
            "0.972",
            "0.969",
            "0.0148",
            "fresh-image",
            "not independent-source replication",
        ):
            assert token.lower() in manuscript.lower()

        # FCP v2 remains a separate measurement-validity study. The current
        # paper imports only the D transport receipt.
        assert "T_white" not in manuscript
        assert "335,994" not in manuscript
        assert "full-pipeline exposure" not in manuscript.lower()


def test_new_phytologist_draft_preserves_frozen_h2_claim() -> None:
    text = MANUSCRIPT.read_text(encoding="utf-8")
    normalized = text.replace("**", "")
    for token in (
        "49,900",
        "377",
        "158 species",
        "0.517",
        "86 species",
        "0.533",
        "p = 0.001",
        "structured-null median of 0.457",
        "exposure-coupled rather than artifact-cleared",
        "same iNaturalist opportunity universe",
        "not an independent-source replication",
        "H2_PROSPECTIVE_WHITE_AXIS_CONFIRMED",
    ):
        assert token.lower() in normalized.lower()



def test_new_phytologist_preserves_postconfirmatory_validity_boundary() -> None:
    text = MANUSCRIPT.read_text(encoding="utf-8")
    result = __import__("json").loads(VALIDITY.read_text(encoding="utf-8"))
    highlight = __import__("json").loads(HIGHLIGHT.read_text(encoding="utf-8"))
    assert result["confirmatory_verdict_changed"] is False
    assert result["inferential_decomposition"]["frozen_observed_W"] == 0.5172457461053418
    assert result["gate_reapplied_structured_null"]["plus_one_upper_p"] == 1 / 300
    assert result["background_white_proxy"]["sample_species"] == 461

    assert highlight["confirmatory_verdict_changed"] is False
    assert highlight["decision"]["state"] == "INDETERMINATE"
    assert highlight["source_identity"]["source_drift_rows"] == 0
    assert highlight["source_identity"]["reacquisition_failed_rows"] == 0
    assert highlight["high_clip"]["rows"] == 2205
    assert highlight["coupling_model"]["odds_ratio"] == __import__("pytest").approx(1.4444932850227639)
    assert highlight["coupling_model"]["or_ci_low"] == __import__("pytest").approx(1.3893320695404507)
    assert highlight["coupling_model"]["or_ci_high"] == __import__("pytest").approx(1.5018445886490097)
    assert highlight["h2_high_clip_sensitivity"]["vector_species"] == 142
    assert highlight["h2_high_clip_sensitivity"]["vector_retention"] == __import__("pytest").approx(142 / 158)
    assert highlight["h2_high_clip_sensitivity"]["observed_W"] == __import__("pytest").approx(0.5034282974026532)
    assert highlight["h2_high_clip_sensitivity"]["structured_null_upper_p"] == __import__("pytest").approx(0.001)
    assert highlight["h2_high_clip_sensitivity"]["support"] is True
    assert ADJUDICATION.exists()

    for token in (
        "increment above a coarse-state-preserving construction baseline",
        "137 (86.7%)",
        "**1.444** (95% CI **1.389–1.502**)",
        "95% CI **1.389–1.502**",
        "142",
        "89.9%",
        "INDETERMINATE",
        "0.968",
        "0.972",
        "exposure-coupled rather than artifact-cleared",
        "not a bitwise numerical reproducer",
    ):
        assert token in text


def test_new_phytologist_reports_bounded_bio5_result_and_failed_transport() -> None:
    import json
    import pytest

    text = MANUSCRIPT.read_text(encoding="utf-8")
    canonical = CANONICAL.read_text(encoding="utf-8")
    env = json.loads(WHITE_ENV.read_text(encoding="utf-8"))
    transport = json.loads(BIO5_TRANSPORT.read_text(encoding="utf-8"))
    observer = json.loads(WHITE_ENV_OBSERVER.read_text(encoding="utf-8"))

    bio5 = next(x for x in env["results"] if x["variable"] == "bio5")
    assert env["eligible_species"] == 281
    assert bio5["median_delta_white_minus_nonwhite_SD"] == pytest.approx(0.0690112924805198)
    assert bio5["wilcoxon_holm_p"] == pytest.approx(0.03544867047368517)
    assert bio5["OR_per_within_species_SD"] == pytest.approx(1.073474538999635)
    assert bio5["p"] == pytest.approx(0.0009188036770296888)
    assert bio5["mechanism_gate_pass"] is True

    assert transport["verdict"] == "LEGACY_BIO5_WHITE_REPLICATION_NOT_SUPPORTED_UNDER_THIS_TEST"
    assert transport["discovery"]["eligible_species"] == 271
    assert transport["discovery"]["species_level"]["wilcoxon_two_sided_p"] == pytest.approx(0.7432522901921289)
    assert transport["reserve"]["eligible_species"] == 260
    assert transport["reserve"]["species_level"]["wilcoxon_two_sided_p"] == pytest.approx(0.054066696426422846)

    paired_bio5 = next(x for x in observer["paired_results"] if x["variable"] == "bio5")
    balanced_bio5 = next(x for x in observer["observer_balanced_results"] if x["variable"] == "bio5")
    assert paired_bio5["n_species"] == 106
    assert paired_bio5["median_delta"] == pytest.approx(0.0)
    assert paired_bio5["wilcoxon_two_sided_p"] == pytest.approx(0.4849619155258311)
    assert paired_bio5["OR_per_within_species_SD"] == pytest.approx(0.7866358392496202)
    assert balanced_bio5["n_species"] == 352
    assert balanced_bio5["median_delta"] == pytest.approx(0.05409791430882366)
    assert balanced_bio5["wilcoxon_two_sided_p"] == pytest.approx(0.0749642018520577)

    for manuscript in (text, canonical):
        for token in (
            "### Post-confirmatory environmental filter and BIO5 transport test",
            "### The prospective BIO5 association is observer-sensitive and does not transport as a common rule",
            "**Holm-adjusted p = 0.0354**",
            "OR = **1.073**",
            "p = **0.000919**",
            "p = **0.743**",
            "p = **0.0541**",
            "median BIO5 contrast was **0.000 SD**",
            "OR = **0.787**",
            "p = **0.0750**",
            "observer conditioning",
            "does not support a common cross-cohort BIO5 rule",
        ):
            assert token in manuscript

    assert "Temperature is therefore not supported as a universal cross-species driver" in text
    assert "some of the within-cohort signal may reflect observer-associated geographic sampling" in text
    assert "rather than in one universal BIO5 coefficient" in text



def test_submission_data_lineage_is_reader_traceable() -> None:
    import json

    text = MANUSCRIPT.read_text(encoding="utf-8")
    lineage = LINEAGE_MAP.read_text(encoding="utf-8")

    for path in (
        LINEAGE_MAP,
        WHITE_ENV_PROTOCOL,
        WHITE_ENV_SCRIPT,
        BIO5_PROTOCOL,
        BIO5_SCRIPT,
        H3B_PROTOCOL,
        H3B_SCRIPT,
        H3B_RESULT,
    ):
        assert path.exists(), f"missing reader-traceable provenance file: {path}"

    for token in (
        "Prospective confirmation cohort",
        "Highlight-validity reacquisition",
        "Climate annotation",
        "POLYMORPHISM_DATA_LINEAGE_MAP_20260925.md",
    ):
        assert token in text

    for token in (
        "5142f7951af0dde5364bb047a566d67e8c479e51",
        "ee854126eed2cfe23e333abe2c28d14df24895a5c52cbc389c060a5a4d6f91f4",
        "0e2ed349122739eecfc725fb2d5e313d284cf30752da91f9a0429cff0eeaa5e6",
        "7e538e5c51c05a7cc47b2fcf53eea92634c8a863",
        "10496492307",
        "10709490106",
        "10292218669",
        "10292662493",
        "10292767459",
        "10292399238",
        "post-confirmatory highlight-validity and environmental analyses",
        "Resolved: expiring Actions artifacts",
        "Resolved: WorldClim provider dependence",
        "fcp-worldclim-2.1-10m-20260925",
    ):
        assert token in lineage

    for archived in (
        ARCHIVE_ROOT / "highlight" / "technical_table.csv.gz",
        ARCHIVE_ROOT / "highlight" / "sealed_join_key.csv",
        ARCHIVE_ROOT / "highlight" / "high_clip_ids.csv",
        ARCHIVE_ROOT / "h3a_tree" / "h3a_s1.tre",
        ARCHIVE_ROOT / "h3a_tree" / "h3a_s2.tre",
        ARCHIVE_ROOT / "h3a_tree" / "h3a_s3.tre",
        ARCHIVE_ROOT / "h3a_covariate" / "sampling_opportunity_preoutcome.csv",
        ARCHIVE_ROOT / "h3a_signal" / "h3a_K_permutation_nulls.csv.gz",
        ARCHIVE_ROOT / "h3b" / "h3b_span_permutation_nulls.csv",
        WORLDCLIM_CHECKSUMS,
        WORLDCLIM_RELEASE_RECEIPT,
    ):
        assert archived.exists(), f"missing permanently frozen archive input: {archived}"

    checksum_text = WORLDCLIM_CHECKSUMS.read_text(encoding="utf-8")
    assert "00513224583665ec0f2f955a4ec252730c4deb2004cce9e793492a3f26df4dcf  wc2.1_10m_bio.zip" in checksum_text
    assert "c72ee7f4f9a0eb4b5f6cd7a003eddc05bda22a2e5666d968ed4e817fd36b9026  wc2.1_10m_srad.zip" in checksum_text
    assert "32e02d85868734c32547da70cf9c620f4739a3cf5dbc4f2b751e16ebb2338ae1  bio/wc2.1_10m_bio_5.tif" in checksum_text

    receipt_text = WORLDCLIM_RELEASE_RECEIPT.read_text(encoding="utf-8")
    for token in (
        "fcp-worldclim-2.1-10m-20260925",
        "49869449",
        "16233364",
        "00513224583665ec0f2f955a4ec252730c4deb2004cce9e793492a3f26df4dcf",
        "c72ee7f4f9a0eb4b5f6cd7a003eddc05bda22a2e5666d968ed4e817fd36b9026",
    ):
        assert token in receipt_text

    h3b = json.loads(H3B_RESULT.read_text(encoding="utf-8"))
    assert h3b["decision"]["verdict"] == "H3B_SAMPLED_SPAN_REPLICATION_NOT_SUPPORTED"
    assert h3b["source_workflow_run"] == 34677468362
    assert h3b["source_artifact_id"] == 10292399238



def test_self_contained_np_provenance_release_is_exposed_to_readers() -> None:
    text = MANUSCRIPT.read_text(encoding="utf-8")
    supporting = SUPPORTING.read_text(encoding="utf-8")
    lineage = LINEAGE_MAP.read_text(encoding="utf-8")

    assert NP_PROVENANCE_RELEASE_RECEIPT.exists()
    receipt = NP_PROVENANCE_RELEASE_RECEIPT.read_text(encoding="utf-8")

    for document in (text, supporting, lineage):
        assert "fcp-np-provenance-20260926" in document
        assert "NP_PROVENANCE_RELEASE_RECEIPT.md" in document

    for token in (
        "fcp-np-provenance-20260926",
        "Asset: fcp-np-provenance-20260926.tar.gz",
        "Source commit: ",
        "Asset bytes: ",
        "Asset SHA256: ",
        "Packaged files: ",
    ):
        assert token in receipt

    source_line = next(line for line in receipt.splitlines() if line.startswith("Source commit: "))
    sha_line = next(line for line in receipt.splitlines() if line.startswith("Asset SHA256: "))
    assert len(source_line.removeprefix("Source commit: ").strip()) == 40
    assert len(sha_line.removeprefix("Asset SHA256: ").strip()) == 64

def test_silene_molecular_anchor_is_quantitative_and_bounded() -> None:
    import json
    import pytest

    text = MANUSCRIPT.read_text(encoding="utf-8")
    molecular = json.loads(SILENE_MOLECULAR.read_text(encoding="utf-8"))

    assert molecular["role"].startswith("structured quantitative extraction")
    assert molecular["source"]["doi"] == "10.3389/fpls.2016.00204"
    assert molecular["source"]["molecular_design"].startswith("mRNA-seq of nine")
    assert molecular["bud_expression"]["shared_significant_locus_in_both_pigmented_vs_white_contrasts"] == "F3h1"
    assert molecular["bud_expression"]["shared_F3h1_fold_changes"]["dark_vs_white"] == pytest.approx(49.0)
    assert molecular["bud_expression"]["shared_F3h1_fold_changes"]["light_vs_white"] == pytest.approx(42.2)
    assert molecular["bud_expression"]["dark_vs_white_significant"][2]["locus"] == "Myb1a"
    assert molecular["bud_expression"]["dark_vs_white_significant"][2]["fold_change_pigmented_over_white"] == pytest.approx(5.1)
    assert molecular["sequence_evidence"]["F3h1_SNPs_in_reported_UTR_CDS_table"] == 0
    assert molecular["sequence_evidence"]["expanded_sequence_survey_individuals"] == 38
    assert molecular["sequence_evidence"]["consistent_colour_differentiating_SNP_after_expansion"] is False
    assert molecular["derived_synthesis"]["causal_variant_identified"] is False

    for token in (
        "### Secondary *Silene littorea* molecular-anchor extraction",
        "### Published molecular data anchor the *Silene* PAL phenotype near F3h1/Myb1a",
        "F3h1 expression was **49.0×** higher",
        "**42.2×** higher in light-pink than white buds",
        "**Myb1a** was **5.1×** higher",
        "**F3h1 had zero SNPs**",
        "expanded sequencing of **38 individuals**",
        "not molecular validation of H2",
        "The frozen extraction is source-derived, not raw-read reanalysis or replication",
        "10.3389/fpls.2016.00204",
    ):
        assert token in text


def test_secondary_pal_wal_result_and_moricandia_interpretation_are_bounded() -> None:
    import json
    import pytest

    text = MANUSCRIPT.read_text(encoding="utf-8")
    pal = json.loads((ROOT / "results" / "polymorphism_silene_decoupling_persistence_20260925" / "result.json").read_text(encoding="utf-8"))
    cross = json.loads((ROOT / "results" / "polymorphism_crossspecies_pal_wal_frequency_20260925" / "result.json").read_text(encoding="utf-8"))

    assert pal["phenotypes"]["PAL"]["positive_frequency_median_percent"] == pytest.approx(15.5)
    assert pal["phenotypes"]["WAL"]["positive_frequency_median_percent"] == pytest.approx(0.21)
    assert cross["PAL"]["lower_bound_median_percent"] == pytest.approx(5.0)
    assert cross["WAL"]["numeric_upper_bound_median_percent"] == pytest.approx(0.1)
    assert cross["WAL"]["numeric_upper_bound_max_percent"] == pytest.approx(1.4)

    for token in (
        "petal anthocyanin-loss (PAL)",
        "Whole-plant anthocyanin-loss (WAL)",
        "HPLC-DAD-MS^n tissue profiling",
        "Within *Silene littorea*, PAL and WAL are not visual labels alone",
        "cross-system entries retain the source table's phenotype classifications",
        "median **15.5%**",
        "median **0.21%**",
        "### Secondary PAL/WAL persistence reanalysis",
        "Published pigment-loss frequency tables",
        "### Biochemically anchored flower-restricted anthocyanin loss reaches higher reported natural frequencies",
        "median PAL lower bound was **5%**",
        "median numeric WAL upper bound was **0.1%**",
        "maintenance filter",
        "ascertained, heterogeneous literature sample",
        "Gómez et al. 2020",
        "Narbona et al. 2026",
        "Lacey 2026",
    ):
        assert token in text

    for reference in (
        "10.1186/s12870-019-2082-6",
        "10.1038/s41467-020-17875-1",
        "10.1002/ajb2.70096",
        "10.1002/ajb2.70106",
    ):
        assert reference in text



def test_pal_wal_bridge_is_explicitly_not_estimable() -> None:
    import json
    import pytest

    overlap = json.loads(PAL_WAL_OVERLAP.read_text(encoding="utf-8"))
    supporting = SUPPORTING.read_text(encoding="utf-8")
    manuscript = MANUSCRIPT.read_text(encoding="utf-8")

    assert overlap["schema"] == "fcp_pal_wal_h2_overlap_feasibility_v2"
    assert overlap["status"] == "complete"
    assert overlap["overlap"]["counts"]["discovery"] == 0
    assert overlap["overlap"]["counts"]["reserve_natural_frequency_context"] == 4
    assert overlap["overlap"]["counts"]["prospective_selected"] == 0
    assert overlap["overlap"]["counts"]["PAL_natural_overlap"] == 2
    assert overlap["overlap"]["counts"]["WAL_natural_overlap"] == 2
    assert overlap["overlap"]["counts"]["WAL_greenhouse_only_overlap"] == 1
    assert overlap["decision"]["pal_vs_wal_bridge_estimable"] is False

    reserve = {row["resolved_species"]: row for row in overlap["overlap"]["reserve_D_eligible_all_registry_contexts"]}
    assert reserve["Gymnadenia rhellicani"]["D"] == pytest.approx(0.65153691467969)
    assert reserve["Silene gallica"]["D"] == pytest.approx(0.234404536862004)
    assert reserve["Delphinium nuttallianum"]["D"] == pytest.approx(0.165925925925926)
    assert reserve["Silene dioica"]["D"] == pytest.approx(0.205515088449532)
    assert reserve["Erythranthe guttata"]["frequency_context"] == "greenhouse_excluded_from_natural_panel"

    for token in (
        "no inferential PAL-versus-WAL comparison of D or H2 geometry is estimable",
        "Gymnadenia rhellicani",
        "Silene gallica",
        "D. nuttallianum",
        "S. dioica",
        "greenhouse-only",
        "two PAL and two natural-context WAL overlaps",
    ):
        assert token in supporting

    for token in (
        "A taxonomy-resolved overlap audit found only **two PAL and two natural-context WAL species**",
        "No PAL-versus-WAL test of D or H2 geometry was therefore estimable",
    ):
        assert token in manuscript


def test_new_phytologist_documents_exact_D_spatial_method_and_methodological_scope() -> None:
    text = MANUSCRIPT.read_text(encoding="utf-8")
    for token in (
        "all retained photograph pairs were used to calculate great-circle geographic distance and flower-colour Jensen–Shannon dissimilarity",
        "rho_i = Spearman(d_geo_ij, d_colour_ij)",
        "999 matched within-species null values",
        "(1 + # {rho_null >= rho_obs}) / 1000",
        "Rank(D) and rank(`rho_i`) are separately residualized",
        "Spearman(d_geo_ij, d_flower_ij - d_background_ij)",
        "It is not the difference between separate flower and background Spearman coefficients",
        "The methodological contribution is architectural rather than a claim to a new standalone statistic",
    ):
        assert token in text

    methods_audit = ROOT / "docs" / "POLYMORPHISM_METHODS_CLASSIFICATION_20260918.md"
    assert methods_audit.exists()
    audit = methods_audit.read_text(encoding="utf-8")
    assert "Standard or widely used components" in audit
    assert "Study-specific design choices" in audit
    assert "The paper should not claim a wholly new statistical method" in audit


def test_new_phytologist_spatial_clue_is_reported_without_causal_upgrade() -> None:
    text = MANUSCRIPT.read_text(encoding="utf-8")
    for token in (
        "Greater D is associated with stronger within-species geographic colour organization",
        "partial rho = **0.0992877**, p = **0.025**",
        "partial rho = **0.1162411**, p = **0.010**",
        "structural rather than causal",
        "two-layer ecological question",
        "Wessinger & Rausher 2012",
        "Lacey 2026",
        "Narbona et al. 2026",
        "shared pigment-network architecture",
    ):
        assert token in text


def test_new_phytologist_cover_letter_exists_and_preserves_claim_boundary() -> None:
    assert COVER.exists(), "New Phytologist cover letter has not been created"
    text = COVER.read_text(encoding="utf-8")
    for token in (
        "New Phytologist",
        "Within-species flower-colour variation shows recurrent achromatic–chromatic geometry across plant species",
        "49,900",
        "377",
        "158",
        "0.517",
        "0.533",
        "same iNaturalist opportunity universe",
        "not an independent-source replication",
        "1.444",
        "INDETERMINATE",
    ):
        assert token.lower() in text.lower()


def test_manuscripts_have_no_control_character_math_corruption() -> None:
    for path in (CANONICAL, MANUSCRIPT):
        text = path.read_text(encoding="utf-8")
        bad = [
            ch for ch in text
            if ord(ch) < 32 and ch not in ("\n", "\r")
        ]
        assert not bad, f"{path.name} contains control characters: {[ord(ch) for ch in bad]}"
        assert "W = mean_i (u_i^T q_white)^2" in text
        assert "p = (1 + #(W_null >= W_obs)) / 1000" in text


def test_cover_letter_answers_three_editor_questions_within_50_words() -> None:
    text = COVER.read_text(encoding="utf-8")
    matches = re.findall(
        r"### Question [123][^\n]*\n\n(.+?)(?=\n\n### Question|\n\n## )",
        text,
        flags=re.S,
    )
    assert len(matches) == 3
    for answer in matches:
        assert len(words(answer)) <= 50


def test_new_phytologist_figure_and_supporting_legends_are_present() -> None:
    text = MANUSCRIPT.read_text(encoding="utf-8")
    for idx in range(1, 6):
        assert f"**Figure {idx}." in text
    for idx in range(1, 9):
        assert f"**Fig. S{idx}." in text
    assert text.index("## References") < text.index("## Supporting Information")


def _section_text(text: str, heading: str, stop_headings: tuple[str, ...]) -> str:
    start_marker = f"## {heading}"
    start = text.index(start_marker) + len(start_marker)
    end = len(text)
    for stop in stop_headings:
        marker = f"## {stop}"
        idx = text.find(marker, start)
        if idx != -1:
            end = min(end, idx)
    return text[start:end]


def _declared_word_count(text: str, label: str) -> int:
    match = re.search(rf"^- {re.escape(label)}: ([0-9,]+) words$", text, flags=re.M)
    assert match, f"missing declared word count for {label}"
    return int(match.group(1).replace(",", ""))


def test_new_phytologist_front_page_counts_are_live() -> None:
    text = MANUSCRIPT.read_text(encoding="utf-8")

    intro = len(words(_section_text(text, "Introduction", ("Materials and Methods",))))
    methods = len(words(_section_text(text, "Materials and Methods", ("Table 1.", "Results"))))
    results = len(words(_section_text(text, "Results", ("Discussion",))))
    discussion = len(words(_section_text(text, "Discussion", ("Acknowledgements",))))
    main = intro + methods + results + discussion

    assert _declared_word_count(text, "Introduction") == intro
    assert _declared_word_count(text, "Materials and Methods") == methods
    assert _declared_word_count(text, "Results") == results
    assert _declared_word_count(text, "Discussion") == discussion
    assert _declared_word_count(text, "Main text (Introduction through Discussion)") == main

    assert discussion / main <= 0.30

    figure_match = re.search(r"^- Figures: ([0-9]+)$", text, flags=re.M)
    table_match = re.search(r"^- Tables: ([0-9]+)$", text, flags=re.M)
    assert figure_match and table_match
    display_items = int(figure_match.group(1)) + int(table_match.group(1))
    assert 6 <= display_items <= 8


def test_new_phytologist_keywords_are_actually_alphabetical() -> None:
    text = MANUSCRIPT.read_text(encoding="utf-8")
    keyword_line = next(
        line for line in text.splitlines()
        if line.startswith("**Keywords (alphabetical):**")
    )
    keywords = [
        item.strip()
        for item in keyword_line.split("**", 2)[-1].split(":", 1)[-1].split(";")
        if item.strip()
    ]
    assert keywords == sorted(keywords, key=str.casefold)


def test_new_phytologist_explains_and_stress_tests_gini_simpson_D() -> None:
    import json
    import pytest

    text = MANUSCRIPT.read_text(encoding="utf-8")
    result = json.loads(D_FINITE.read_text(encoding="utf-8"))

    assert result["rank_stability"]["validation"]["spearman_raw_vs_corrected"] == pytest.approx(0.999978412058831)
    assert result["D_spatial"]["validation_primary"]["corrected"]["p_upper_geometry_preserving_spatial_null"] == pytest.approx(0.025)

    for token in (
        "Gini–Simpson diversity",
        "probability that two observations belong to different colour states",
        "D_{\\mathrm{corr}}",
        "40/39 = 1.0256",
        "### Finite-sample correction does not alter D-based conclusions",
        "0.999977",
        "0.999978",
        "0.099561",
        "0.116155",
    ):
        assert token in text
