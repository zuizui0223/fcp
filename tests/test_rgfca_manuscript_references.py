"""Internal citation consistency, not automated verification of scientific truth.

Read only the active manuscript and its bounded primary-reading notes. Never open
reserve outcomes, images or inference inputs. The bounded reference set must be
extended deliberately alongside a source audit, not by accepting arbitrary DOIs.
"""
from collections import Counter
from pathlib import Path
import re

import pytest

ROOT = Path(__file__).resolve().parents[1]
DOI_LINK = re.compile(r"\[([^\]]+)\]\(https://doi\.org/(10\.\d{4,9}/[^)\s]+)\)", re.I)
# Issue year and final published title, manually checked against primary sources.
REFERENCES = {
    "10.1093/biosci/biab093": (2021, "Observing the Observers: How Participants Contribute Data to iNaturalist and Implications for Biodiversity Science"),
    "10.1098/rsif.2015.0086": (2015, "Range bagging: a new method for ecological niche modelling from presence-only data"),
    "10.1111/2041-210x.12018": (2013, "Dismantling the Mantel tests"),
    "10.1002/ece3.7307": (2021, "Is color data from citizen science photographs reliable for biodiversity research?"),
    "10.1002/aps3.11546": (2023, "A pipeline for the rapid collection of color data from photographs"),
    "10.1016/j.cub.2025.03.035": (2025, "Delayed flowering phenology of red-flowering plants in response to hummingbird migration"),
    "10.1086/739413": (2026, "High-Throughput iNaturalist Image Analysis Reveals Flower Color Divergence in Monarda fistulosa"),
    "10.2202/1544-6115.1585": (2010, "Permutation P-values should never be zero: calculating exact P-values when permutations are randomly drawn"),
    "10.1111/ecog.02881": (2017, "Cross-validation strategies for data with temporal, spatial, hierarchical, or phylogenetic structure"),
    "10.1111/2041-210x.13525": (2021, "Spatial thinning and class balancing: Key choices lead to variation in the performance of species distribution models with citizen science data"),
    "10.1111/nph.18361": (2022, "The ecological implications of intra- and inter-species variation in phenological sensitivity"),
    "10.1038/s41586-020-2649-2": (2020, "Array programming with NumPy"),
    "10.1038/s41592-019-0686-2": (2020, "SciPy 1.0: fundamental algorithms for scientific computing in Python"),
    "10.25080/majora-92bf1922-00a": (2010, "Data Structures for Statistical Computing in Python"),
    "10.1109/mcse.2007.55": (2007, "Matplotlib: A 2D graphics environment"),
}


def validate_citations(manuscript, audit_text):
    if manuscript.count("\n## References\n") != 1:
        raise ValueError("Exactly one reference section required")
    body, bibliography = manuscript.split("\n## References\n")
    cited = DOI_LINK.findall(body)
    refs = DOI_LINK.findall(bibliography)
    ref_counts = Counter(doi.lower() for _, doi in refs)
    if set(ref_counts) != set(REFERENCES) or set(ref_counts.values()) != {1}:
        raise ValueError("Reference set or unique DOI census differs from audited set")
    if {doi.lower() for _, doi in cited} != set(ref_counts):
        raise ValueError("Uncited reference or citation without a reference")
    audited_dois = {doi.lower() for _, doi in DOI_LINK.findall(audit_text)}
    if not set(ref_counts) <= audited_dois:
        raise ValueError("Missing primary-source reading record")
    rows = [row for row in bibliography.splitlines() if row.startswith("- ")]
    if len(rows) != len(refs):
        raise ValueError("Each reference must have one entry")
    for row in rows:
        row_links = DOI_LINK.findall(row)
        if len(row_links) != 1:
            raise ValueError("Each reference must have one DOI")
        label, doi = row_links[0]
        if label.lower() != doi.lower():
            raise ValueError("Reference DOI label and target differ")
        year, title = REFERENCES[doi.lower()]
        if f"({year})" not in row or title not in row:
            raise ValueError("Final publication title or year mismatch")
    for label, doi in cited:
        if str(REFERENCES[doi.lower()][0]) not in label:
            raise ValueError("In-text publication year mismatch")
    return len(ref_counts)


@pytest.fixture
def documents():
    manuscript = (ROOT / "docs/RGFCA_MANUSCRIPT.md").read_text(encoding="utf-8")
    audits = "\n".join((ROOT / "docs" / name).read_text(encoding="utf-8") for name in (
        "RGFCA_IMAGE_ECOLOGY_LITERATURE_AUDIT.md", "RGFCA_STATISTICAL_LITERATURE_AUDIT.md",
        "RGFCA_SCIENTIFIC_SOFTWARE_AUDIT.md"))
    return manuscript, audits


def test_all_core_references_are_cited_once_in_list_and_audited(documents):
    assert validate_citations(*documents) == 15


def test_orphan_reference_rejected(documents):
    manuscript, audits = documents
    manuscript = manuscript.replace("[Drake (2015)](https://doi.org/10.1098/rsif.2015.0086)", "Drake (2015)")
    with pytest.raises(ValueError, match="Uncited"):
        validate_citations(manuscript, audits)


def test_duplicate_reference_rejected(documents):
    manuscript, audits = documents
    row = next(row for row in manuscript.splitlines() if row.startswith("- Drake,"))
    with pytest.raises(ValueError, match="unique DOI"):
        validate_citations(manuscript + "\n" + row, audits)


@pytest.mark.parametrize("before,after", [
    ("10.1086/739413", "10.1101/2025.05.21.655392"),
    ("(2026)", "(2025)"),
    ("Delayed flowering phenology of red-flowering plants in response to hummingbird migration",
     "A massive community-science dataset reveals convergent evolution of delayed flowering phenology"),
])
def test_old_preprint_or_wrong_publication_metadata_rejected(documents, before, after):
    manuscript, audits = documents
    with pytest.raises(ValueError, match="Reference set|title or year"):
        validate_citations(manuscript.replace(before, after), audits)


def test_missing_reading_record_rejected(documents):
    manuscript, audits = documents
    with pytest.raises(ValueError, match="reading record"):
        validate_citations(manuscript, audits.replace("10.1086/739413", "missing"))


def test_same_publication_can_be_cited_twice_in_body(documents):
    manuscript, audits = documents
    manuscript = manuscript.replace("## References", "[Drake (2015)](https://doi.org/10.1098/rsif.2015.0086).\n\n## References")
    assert validate_citations(manuscript, audits) == 15


def test_audit_does_not_close_replication_or_submission(documents):
    manuscript, audits = documents
    for required in ("not submission-ready", "replication and", "are pending",
                     "not a systematic review", "software/data-provider bibliography",
                     "not scientific validity"):
        assert required in manuscript
    assert "No reserve outcome file was opened" in audits
    assert "colour-assignment Methods were **not verified**" in audits
    assert "No Access" in audits
