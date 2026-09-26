from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ACTIVE = ROOT / ".github" / "workflows"
ARCHIVE = ROOT / "archive" / "workflows"

EXPECTED_ACTIVE = {
    "archive-worldclim-release-20260925.yml",
    "build-np-provenance-snapshot-20260926.yml",
    "freeze-fcp-expiring-artifacts-20260925.yml",
    "legacy-white-bio5-replication-20260925.yml",
    "polymorphism-h3b-reserve-span.yml",
    "polymorphism-manuscript-claim-guard.yml",
    "polymorphism-new-phytologist-submission-guard.yml",
    "polymorphism-publication-figure-tests.yml",
    "regenerate-polymorphism-figures-reader-flow-20260926.yml",
    "repository-layout-guard.yml",
    "third-cohort-highlight-terminal-recovery-20260923.yml",
    "third-cohort-highlight-validity-20260922.yml",
    "white-environment-mechanism-20260925.yml",
    "worldclim-archive-checksum-freeze-20260925.yml",
}

REQUIRED_CURRENT_FILES = (
    "README.md",
    "CURRENT_PAPER_REPRODUCIBILITY.md",
    "docs/POLYMORPHISM_MANUSCRIPT_NEW_PHYTOLOGIST.md",
    "docs/POLYMORPHISM_MANUSCRIPT.md",
    "docs/POLYMORPHISM_CURRENT_CLAIM_LEDGER_20260918.md",
    "docs/POLYMORPHISM_DATA_LINEAGE_MAP_20260925.md",
    "docs/POLYMORPHISM_SUPPORTING_INFORMATION_20260918.md",
    "docs/POLYMORPHISM_PROVENANCE_RELEASE_MANIFEST_20260926.md",
    "archive/fcp_submission_20260925/NP_PROVENANCE_RELEASE_RECEIPT.md",
)

def test_current_paper_entrypoints_are_stable() -> None:
    for rel in REQUIRED_CURRENT_FILES:
        assert (ROOT / rel).exists(), rel

def test_active_actions_surface_is_intentionally_small() -> None:
    names = {p.name for p in ACTIVE.glob("*.yml")}
    assert names == EXPECTED_ACTIVE
    assert not any(name.startswith("jbi-") for name in names)
    assert not any(name.startswith("disttrait-") for name in names)
    assert "global-literature-discovery.yml" not in names
    assert "p500-prospective-location-blind-measurement.yml" not in names
    assert "manuscript-consistency.yml" not in names

def test_historical_workflows_are_preserved_but_inactive() -> None:
    archived = list(ARCHIVE.rglob("*.yml"))
    assert len(archived) >= 75
    assert len(list((ARCHIVE / "jbi").glob("*.yml"))) >= 60
    assert len(list((ARCHIVE / "disttrait").glob("*.yml"))) >= 13
    assert (ARCHIVE / "literature" / "global-literature-discovery.yml").exists()
    assert (ARCHIVE / "p500" / "p500-prospective-location-blind-measurement.yml").exists()

def test_root_routes_reader_to_current_reproducibility_contract() -> None:
    readme = (ROOT / "README.md").read_text(encoding="utf-8")
    for token in (
        "CURRENT_PAPER_REPRODUCIBILITY.md",
        "POLYMORPHISM_MANUSCRIPT_NEW_PHYTOLOGIST.md",
        "POLYMORPHISM_DATA_LINEAGE_MAP_20260925.md",
        "archive/workflows/",
        "NP_PROVENANCE_RELEASE_RECEIPT.md",
    ):
        assert token in readme

def test_provenance_builder_tracks_current_paper_surface() -> None:
    workflow = (ACTIVE / "build-np-provenance-snapshot-20260926.yml").read_text(encoding="utf-8")
    for token in (
        "CURRENT_PAPER_REPRODUCIBILITY.md",
        "docs/POLYMORPHISM_MANUSCRIPT_NEW_PHYTOLOGIST.md",
        "docs/POLYMORPHISM_SUPPORTING_INFORMATION_20260918.md",
        "docs/figures/polymorphism_20260918/**",
        "results/polymorphism_*/**",
        "README.md",
    ):
        assert token in workflow
