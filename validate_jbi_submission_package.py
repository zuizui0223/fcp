#!/usr/bin/env python3
"""Validate the Journal of Biogeography manuscript package without changing it.

The validator checks journal-facing structure, numerical guardrails, supporting-file
provenance, DOI-preparation documentation and unresolved submission fields. It is
intentionally strict: failures identify inconsistent or unsafe claims, whereas
unresolved author-controlled fields are reported as blockers rather than invented.
"""
from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent
DOCS = ROOT / "docs"
MANUSCRIPT = DOCS / "jbi_manuscript_editorial_revision.md"
INDEX = DOCS / "jbi_supporting_information_index.md"
EDITORIAL = DOCS / "jbi_editorial_submission_check.md"
CHECKLIST = DOCS / "jbi_submission_completion_checklist.md"
TITLE_PAGE = DOCS / "jbi_title_page_template.md"
COVER_LETTER = DOCS / "jbi_cover_letter_template.md"
AUTHOR_FORM = DOCS / "jbi_author_confirmation_form.md"
SEARCH_PROVENANCE = DOCS / "jbi_literature_search_provenance.md"
TAXON_IMAGE = DOCS / "jbi_taxon_image_candidate.md"
GBIF_PROTOCOL = DOCS / "jbi_gbif_doi_protocol.md"
ARCHIVE_PROTOCOL = DOCS / "jbi_archive_release_protocol.md"
ZENODO_TEMPLATE = DOCS / "jbi_zenodo_metadata_template.json"
REPORT = DOCS / "jbi_submission_validation_report.json"
GBIF_BUNDLE = DOCS / "supporting/jbi_gbif_doi_bundle"
EXACT_GBIF_SHA = "f25ae0cf2c84c45ae461a932d6c6063edda64591913a2495e4a3da82d573f094"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def require(condition: bool, message: str, failures: list[str]) -> None:
    if not condition:
        failures.append(message)


def section(text: str, start: str, end: str) -> str:
    i = text.find(start)
    if i < 0:
        raise ValueError(f"Missing start marker: {start}")
    j = text.find(end, i + len(start))
    if j < 0:
        raise ValueError(f"Missing end marker: {end}")
    return text[i:j]


def word_count(text: str) -> int:
    return len(re.findall(r"\b[\w’'-]+\b", text, flags=re.UNICODE))


def parse_supporting_index(text: str) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    pattern = re.compile(
        r"^\| S(?P<id>\d+) \| `(?P<path>[^`]+)` \| (?P<rows>\d+) \| `(?P<sha>[0-9a-f]{64})` \|$",
        re.MULTILINE,
    )
    for match in pattern.finditer(text):
        rows.append(
            {
                "id": int(match.group("id")),
                "path": match.group("path"),
                "rows": int(match.group("rows")),
                "sha256": match.group("sha"),
            }
        )
    return rows


def csv_rows(path: Path) -> int:
    with path.open(encoding="utf-8", errors="replace") as handle:
        return max(0, sum(1 for _ in handle) - 1)


def main() -> None:
    failures: list[str] = []
    warnings: list[str] = []

    required_files = [
        MANUSCRIPT,
        INDEX,
        EDITORIAL,
        CHECKLIST,
        TITLE_PAGE,
        COVER_LETTER,
        AUTHOR_FORM,
        SEARCH_PROVENANCE,
        TAXON_IMAGE,
        GBIF_PROTOCOL,
        ARCHIVE_PROTOCOL,
        ZENODO_TEMPLATE,
        DOCS / "figures/moisture_effect_strict_vs_enriched.svg",
        DOCS / "figures/moisture_leave_one_family_out.svg",
        DOCS / "supporting/jbi_opentree_induced_topology.tre",
        DOCS / "supporting/jbi_dated_phylogeny_s1.tre",
        DOCS / "supporting/jbi_dated_phylogeny_s2.tre",
        DOCS / "supporting/jbi_dated_phylogeny_s3.tre",
        DOCS / "supporting/jbi_phylogenetic_sensitivity_manifest.json",
        DOCS / "supporting/jbi_dated_phylogeny_manifest.json",
        GBIF_BUNDLE / "jbi_gbif_exact_occurrence_subset.csv.gz",
        GBIF_BUNDLE / "jbi_gbif_parent_dataset_counts.csv",
        GBIF_BUNDLE / "jbi_gbif_exact_species_counts.csv",
        GBIF_BUNDLE / "jbi_gbif_broad_download_request.json",
        GBIF_BUNDLE / "jbi_gbif_derived_dataset_metadata_template.json",
        GBIF_BUNDLE / "jbi_gbif_doi_bundle_manifest.json",
        ROOT / "prepare_jbi_submission_release.py",
        ROOT / "validate_jbi_gbif_doi_bundle.py",
    ]
    for path in required_files:
        require(path.is_file() and path.stat().st_size > 0, f"Missing or empty file: {path}", failures)

    manuscript = MANUSCRIPT.read_text(encoding="utf-8")
    index_text = INDEX.read_text(encoding="utf-8")
    editorial = EDITORIAL.read_text(encoding="utf-8")
    checklist = CHECKLIST.read_text(encoding="utf-8")
    title_page = TITLE_PAGE.read_text(encoding="utf-8")
    cover_letter = COVER_LETTER.read_text(encoding="utf-8")
    author_form = AUTHOR_FORM.read_text(encoding="utf-8")
    provenance = SEARCH_PROVENANCE.read_text(encoding="utf-8")
    taxon_image = TAXON_IMAGE.read_text(encoding="utf-8")
    gbif_protocol = GBIF_PROTOCOL.read_text(encoding="utf-8")
    archive_protocol = ARCHIVE_PROTOCOL.read_text(encoding="utf-8")
    zenodo_template = ZENODO_TEMPLATE.read_text(encoding="utf-8")

    title_match = re.search(r"^# (.+)$", manuscript, re.MULTILINE)
    running_match = re.search(r"^\*\*Running title:\*\* (.+)$", manuscript, re.MULTILINE)
    require(title_match is not None, "Manuscript title missing", failures)
    require(running_match is not None, "Running title missing", failures)
    title = title_match.group(1).strip() if title_match else ""
    running_title = running_match.group(1).strip() if running_match else ""
    require(len(title) <= 115, f"Title exceeds 115 characters: {len(title)}", failures)
    require(len(running_title) < 40, f"Running title is not below 40 characters: {len(running_title)}", failures)

    abstract = section(manuscript, "## Abstract", "**Keywords:**")
    abstract_no_headings = re.sub(r"^#{2,3} .+$", "", abstract, flags=re.MULTILINE)
    abstract_words = word_count(abstract_no_headings)
    require(abstract_words <= 300, f"Abstract exceeds 300 words: {abstract_words}", failures)
    for heading in ("Aim", "Location", "Taxon", "Methods", "Results", "Main conclusions"):
        require(f"### {heading}" in abstract, f"Abstract heading missing: {heading}", failures)

    keyword_match = re.search(r"^\*\*Keywords:\*\* (.+)$", manuscript, re.MULTILINE)
    require(keyword_match is not None, "Keywords line missing", failures)
    keywords = [x.strip() for x in keyword_match.group(1).split(",")] if keyword_match else []
    require(6 <= len(keywords) <= 10, f"Keyword count outside 6–10: {len(keywords)}", failures)
    require(keywords == sorted(keywords, key=str.lower), "Keywords are not alphabetized", failures)

    required_main_sections = [
        "## Introduction",
        "## Methods",
        "## Results",
        "## Discussion",
        "## References",
        "## Data Accessibility Statement",
        "## Tables",
        "## Figure legends and embedded figures",
        "## Supporting Information",
    ]
    for heading in required_main_sections:
        count = len(re.findall(rf"^{re.escape(heading)}$", manuscript, flags=re.MULTILINE))
        require(count == 1, f"Expected one main section: {heading}; found {count}", failures)

    numerical_guardrails = [
        "odds ratio = 0.426",
        "0.184–0.985",
        "p = 0.0556",
        "odds ratio of 0.300",
        "0.133–0.675",
        "p = 0.0164",
        "odds ratio of 0.592",
        "0.244–1.434",
        "0.472",
        "0.175–1.272",
        "0.464 to 0.470",
        "0.366 to 0.369",
        "from 16 to 19 July 2026",
        "Tables S16–S17",
        "Table S18",
        "Table S19",
        "Appendix S1",
        "58,455-row occurrence subset",
        "389 parent GBIF datasets",
        EXACT_GBIF_SHA,
    ]
    for phrase in numerical_guardrails:
        require(phrase in manuscript, f"Required manuscript guardrail missing: {phrase}", failures)

    forbidden_claims = [
        r"\bdemonstrates?\b",
        r"\bconfirms?\b",
        r"\bclimate drives\b",
        r"\blocal adaptation explains\b",
        r"\bphylogenetically robust\b",
        r"\bindependent of ancestry\b",
        r"\brobust evidence\b",
    ]
    discussion_body = section(manuscript, "## Discussion", "## References")
    for pattern in forbidden_claims:
        require(re.search(pattern, discussion_body, flags=re.IGNORECASE) is None, f"Forbidden claim in Discussion: {pattern}", failures)

    expected_references = [
        "Dalrymple, R. L.",
        "Fick, S. E.",
        "Hinchliff, C. E.",
        "Ho, L. S. T.",
        "Jin, Y., & Qian, H. (2019)",
        "Jin, Y., & Qian, H. (2022)",
        "Koski, M. H.",
        "Narbona, E.",
        "Priem, J.",
        "Rausher, M. D.",
        "Seabold, S.",
        "Smith, S. A.",
        "Trunschke, J.",
        "Wessinger, C. A.",
        "Zanne, A. E.",
    ]
    references = section(manuscript, "## References", "## Data Accessibility Statement")
    for reference in expected_references:
        require(reference in references, f"Expected reference missing: {reference}", failures)

    supporting_rows = parse_supporting_index(index_text)
    require([row["id"] for row in supporting_rows] == list(range(1, 20)), "Supporting index must contain S1–S19 in order", failures)
    supporting_checks: list[dict[str, object]] = []
    for row in supporting_rows:
        path = ROOT / str(row["path"])
        exists = path.is_file()
        actual_rows = csv_rows(path) if exists else None
        actual_sha = sha256(path) if exists else None
        require(exists, f"Supporting table missing: {path}", failures)
        if exists:
            require(actual_rows == row["rows"], f"Row-count mismatch for {path}: {actual_rows} != {row['rows']}", failures)
            require(actual_sha == row["sha256"], f"SHA-256 mismatch for {path}: {actual_sha} != {row['sha256']}", failures)
        supporting_checks.append({**row, "exists": exists, "actual_rows": actual_rows, "actual_sha256": actual_sha})

    require("Appendix S1" in index_text and "jbi_literature_search_provenance.md" in index_text, "Appendix S1 missing from Supporting index", failures)
    require("jbi_table_s18_blinded_classification_review.csv" in index_text, "S18 missing from Supporting index", failures)
    require("jbi_table_s19_rule_classification_key.csv" in index_text, "S19 missing from Supporting index", failures)
    require("jbi_classification_review_protocol.md" in index_text, "Classification review protocol missing from Supporting index", failures)
    require("jbi_classification_rule_audit.md" in index_text, "Classification rule audit missing from Supporting index", failures)
    require("jbi_gbif_exact_occurrence_subset.csv.gz" in index_text, "Exact GBIF archive missing from Supporting index", failures)
    require(EXACT_GBIF_SHA in index_text, "Exact GBIF archive SHA missing from Supporting index", failures)
    require("16 to 19 July 2026" in provenance, "Search-date range missing from provenance appendix", failures)
    require("must not be reconstructed retrospectively" in provenance, "Human-review boundary missing from provenance appendix", failures)

    require(re.search(r"GBIF\s+derived[- ]dataset\s+DOI", checklist, flags=re.IGNORECASE) is not None, "GBIF Derived Dataset DOI blocker missing from submission checklist", failures)
    require("human screening/classification" in checklist.lower(), "Human-review blocker missing from submission checklist", failures)
    require("prepare_jbi_submission_release.py --strict" in checklist, "Strict release gate missing from submission checklist", failures)
    require("Human evidence screening and classification" in author_form, "Human-review section missing from author confirmation form", failures)
    require("CRediT author contributions" in author_form, "CRediT section missing from author confirmation form", failures)
    require("Final submission sign-off" in author_form, "Final sign-off missing from author confirmation form", failures)
    require("Why two GBIF citation records are needed" in gbif_protocol, "GBIF citation boundary missing from protocol", failures)
    require(EXACT_GBIF_SHA in gbif_protocol, "Exact GBIF archive SHA missing from protocol", failures)
    require("python prepare_jbi_submission_release.py --strict" in archive_protocol, "Strict release command missing from archive protocol", failures)
    require("REPLACE_WITH_COMMIT_SHA" in zenodo_template, "Zenodo commit placeholder missing", failures)
    require("CC0" in taxon_image and "Ipomoea purpurea" in taxon_image, "Verified taxon-image candidate missing", failures)
    require("AUTHOR CONFIRMATION REQUIRED" in cover_letter, "Cover-letter author guardrails missing", failures)

    unresolved_counts = {
        "manuscript_not_verified": manuscript.lower().count("not verified"),
        "editorial_not_verified": editorial.lower().count("not verified"),
        "checklist_not_verified": checklist.lower().count("not verified"),
        "title_page_not_verified": title_page.lower().count("not verified"),
        "cover_letter_confirmation_markers": cover_letter.count("AUTHOR CONFIRMATION REQUIRED"),
        "zenodo_replace_markers": zenodo_template.count("REPLACE_"),
    }
    if unresolved_counts["manuscript_not_verified"]:
        warnings.append(f"Anonymized manuscript still contains {unresolved_counts['manuscript_not_verified']} `Not verified` markers")

    manuscript_body = section(manuscript, "## Introduction", "## References")
    manuscript_words = word_count(manuscript_body)
    if manuscript_words > 6000:
        warnings.append(f"Introduction–Discussion word count exceeds 6,000: {manuscript_words}")

    report = {
        "status": "pass" if not failures else "fail",
        "title_characters": len(title),
        "running_title_characters": len(running_title),
        "abstract_words": abstract_words,
        "keywords": keywords,
        "introduction_through_discussion_words": manuscript_words,
        "supporting_tables": supporting_checks,
        "exact_gbif_archive_sha256": EXACT_GBIF_SHA,
        "unresolved_counts": unresolved_counts,
        "warnings": warnings,
        "failures": failures,
    }
    REPORT.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, ensure_ascii=False))
    if failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
