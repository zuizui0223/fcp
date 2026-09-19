from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MANUSCRIPT = ROOT / "docs" / "DISTTRAIT_METHODS_MANUSCRIPT_V0_1_20260919.md"
AUDIT = ROOT / "docs" / "DISTTRAIT_METHODS_LITERATURE_AUDIT_20260919.md"


def test_disttrait_methods_manuscript_has_audited_literature_and_clean_notation() -> None:
    raw = MANUSCRIPT.read_bytes()
    assert b"\\r" not in raw

    text = raw.decode("utf-8")
    required = (
        "Bolnick et al. 2011",
        "Siefert et al. 2015",
        "Carmona et al. 2016",
        "Carmona 2019",
        "Bird et al. 2014",
        "Di Cecco et al. 2021",
        "Laitly et al. 2021",
        "Phipson & Smyth 2010",
        "Legendre & Fortin 2010",
        "Guillot & Rousset 2013",
        "DerSimonian & Laird 1986",
        "Viechtbauer 2010",
        "Röver et al. 2015",
        "Gower 1971",
        "Martorelli et al. 2025",
        "Shokri et al. (2022)",
        "`D = 1 - sum_k(p_k^2)`",
        "`rho_i = Spearman(d_geo[j,k], d_trait[j,k])`",
        "Cochran `Q` and `tau^2`",
        "## References",
    )
    for token in required:
        assert token in text

    doi_tokens = (
        "10.1016/j.tree.2011.01.009",
        "10.1111/ele.12508",
        "10.1016/j.tree.2016.02.003",
        "10.1002/ecy.2876",
        "10.1016/j.biocon.2013.07.037",
        "10.1093/biosci/biab093",
        "10.1002/ece3.7307",
        "10.2202/1544-6115.1585",
        "10.1111/j.1755-0998.2010.02866.x",
        "10.1111/2041-210x.12018",
        "10.1016/0197-2456(86)90046-2",
        "10.18637/jss.v036.i03",
        "10.1186/s12874-015-0091-1",
        "10.2307/2528823",
        "10.5281/zenodo.16537297",
        "10.1242/jeb.244842",
    )
    for doi in doi_tokens:
        assert doi in text

    assert AUDIT.exists()
    audit = AUDIT.read_text(encoding="utf-8")
    assert "do not cite Mantel critiques as proof" in audit
    assert "do not cite random-effects literature as proof" in audit


def test_disttrait_methods_manuscript_preserves_claim_boundaries() -> None:
    text = MANUSCRIPT.read_text(encoding="utf-8")
    for token in (
        "should not be described as a wholly new statistical family",
        "does not guarantee universal type-I error control",
        "could not identify whether the association was",
        "biological or selection-induced",
        "metric choice is part of the estimand",
        "Benchmark performance",
        "claims remain supported by frozen repository receipts",
    ):
        assert token in text
