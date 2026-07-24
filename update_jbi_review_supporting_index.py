#!/usr/bin/env python3
"""Add or refresh S18-S19 and their provenance in the JBI package."""
from __future__ import annotations

import hashlib
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent
INDEX = ROOT / "docs/jbi_supporting_information_index.md"
MANUSCRIPT = ROOT / "docs/jbi_manuscript_editorial_revision.md"
S18 = ROOT / "docs/supporting/jbi_table_s18_blinded_classification_review.csv"
S19 = ROOT / "docs/supporting/jbi_table_s19_rule_classification_key.csv"
PROTOCOL = ROOT / "docs/jbi_classification_review_protocol.md"
AUDIT = ROOT / "docs/jbi_classification_rule_audit.md"
BUILDER = ROOT / "build_jbi_blinded_classification_review.py"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def data_rows(path: Path) -> int:
    return max(0, len(path.read_text(encoding="utf-8").splitlines()) - 1)


def main() -> None:
    for path in (INDEX, MANUSCRIPT, S18, S19, PROTOCOL, AUDIT, BUILDER):
        if not path.is_file() or path.stat().st_size == 0:
            raise SystemExit(f"Missing or empty review-supporting file: {path}")

    if data_rows(S18) != 34 or data_rows(S19) != 34:
        raise SystemExit(f"S18/S19 must each contain 34 data rows: {data_rows(S18)}, {data_rows(S19)}")

    index = INDEX.read_text(encoding="utf-8")
    index = re.sub(
        r"(?m)^\| S18 \|.*\n|^\| S19 \|.*\n",
        "",
        index,
    )
    s17_pattern = re.compile(r"(?m)^(\| S17 \|[^\n]+\n)")
    if len(s17_pattern.findall(index)) != 1:
        raise SystemExit("Expected exactly one S17 row in Supporting Information index")
    additions = (
        f"| S18 | `docs/supporting/{S18.name}` | {data_rows(S18)} | `{sha256(S18)}` |\n"
        f"| S19 | `docs/supporting/{S19.name}` | {data_rows(S19)} | `{sha256(S19)}` |\n"
    )
    index = s17_pattern.sub(r"\1" + additions, index, count=1)

    old_provenance = (
        "Tables S1–S7 derive from workflow run `29972327794`. Tables S8–S15 derive from successful Open Tree phylogenetic-sensitivity run `30067762848`, artifact `8586932030`, digest `sha256:a3ce368fa0dc42bcc26edfca7f09286a8bfe8b609d1b9e58fc75b6f096baf16f`. Tables S16–S17 derive from the fixed-seed dated-megaphylogeny run `30076757379`, artifact `8590190840`, digest `sha256:8f11f59a12758f67124647f719fcc79532651c0512f9e0c199a6afa80d178a68`."
    )
    new_provenance = old_provenance + (
        " Tables S18–S19 are generated from the frozen S6 manifest and the resolved evidence queue by "
        "`build_jbi_blinded_classification_review.py`; S18 hides the rule label and remains blank for "
        "independent review, whereas S19 is the post-review comparison key."
    )
    if old_provenance not in index:
        if "Tables S18–S19 are generated" not in index:
            raise SystemExit("Supporting provenance paragraph changed unexpectedly")
    else:
        index = index.replace(old_provenance, new_provenance, 1)

    additional_anchor = "Additional source-backed files:\n"
    if index.count(additional_anchor) != 1:
        raise SystemExit("Additional source-backed files anchor missing")
    review_lines = (
        f"\n- `docs/jbi_classification_review_protocol.md`: blinded review instructions, operational labels and adjudication requirements; SHA-256 `{sha256(PROTOCOL)}`.\n"
        f"- `docs/jbi_classification_rule_audit.md`: plural-expression code-parity audit documenting zero frozen-label changes; SHA-256 `{sha256(AUDIT)}`.\n"
    )
    index = re.sub(
        r"\n- `docs/jbi_classification_review_protocol\.md`:[^\n]*\n",
        "\n",
        index,
    )
    index = re.sub(
        r"\n- `docs/jbi_classification_rule_audit\.md`:[^\n]*\n",
        "\n",
        index,
    )
    index = index.replace(additional_anchor, additional_anchor + review_lines, 1)
    INDEX.write_text(index, encoding="utf-8")

    manuscript = MANUSCRIPT.read_text(encoding="utf-8")
    current = re.search(r"(?ms)^## Supporting Information\n.*\Z", manuscript)
    if current is None:
        raise SystemExit("Supporting Information section missing from manuscript")
    replacement = """## Supporting Information
Tables S1–S7 contain the original model matrices and audits. **Table S8** is the primary baseline model dataset; **Table S9** is the paginated model dataset; **Table S10** is the GBIF taxonomic and retention audit; **Tables S11–S12** report paginated robustness and family deletion; **Tables S13–S15** report Open Tree phylogenetic summaries, name resolution and all 200 replicate fits; **Table S16** reports the six fixed-seed dated-megaphylogeny models; and **Table S17** records V.PhyloMaker2 species placement. **Table S18** is a 34-row blinded human-review template that hides the frozen rule label, and **Table S19** is the separate rule-label and source-diagnostic key to be opened only after first-pass review. **Appendix S1** reconstructs the automated literature-search chronology from repository commits and states which human-review details remain unrecorded. The Supporting Information index also identifies the classification review protocol and rule audit, GBIF QC manifests, the exact occurrence and DOI-preparation bundle, Open Tree topology, dated S1–S3 trees and both phylogenetic provenance manifests.
"""
    manuscript = manuscript[: current.start()] + replacement
    MANUSCRIPT.write_text(manuscript, encoding="utf-8")

    print(
        {
            "status": "pass",
            "s18_rows": data_rows(S18),
            "s18_sha256": sha256(S18),
            "s19_rows": data_rows(S19),
            "s19_sha256": sha256(S19),
            "review_protocol_sha256": sha256(PROTOCOL),
            "rule_audit_sha256": sha256(AUDIT),
        }
    )


if __name__ == "__main__":
    main()
