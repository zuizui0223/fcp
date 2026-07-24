#!/usr/bin/env python3
from pathlib import Path

path = Path("docs/jbi_manuscript_editorial_revision.md")
text = path.read_text(encoding="utf-8")

old = """## Data Accessibility Statement

Analysis code, source-level evidence fields, frozen classification manifests, correction logs and model outputs are maintained in the public GitHub repository `zuizui0223/fcp`. The exact files supporting the manuscript tables are indexed in `docs/jbi_supporting_information_index.md`.

**Not verified:** before submission, archive the exact code and data release in a permanent repository and insert its DOI here. A citable GBIF derived-dataset DOI for the occurrence data must also be added. The GitHub repository alone should not be treated as the final preservation record.
"""
new = """## Data Accessibility Statement

Analysis code, source-level evidence fields, frozen classification manifests, correction logs and model outputs are maintained in the public GitHub repository `zuizui0223/fcp`. The exact files supporting the manuscript tables are indexed in `docs/jbi_supporting_information_index.md`. The 58,455-row occurrence subset used for the paginated sensitivity analysis is frozen as `jbi_gbif_exact_occurrence_subset.csv.gz` (SHA-256 `f25ae0cf2c84c45ae461a932d6c6063edda64591913a2495e4a3da82d573f094`) together with contribution counts for 389 parent GBIF datasets, a broad-download request and Derived Dataset registration metadata. The distinction between the broad GBIF retrieval and the locally capped and coordinate-deduplicated exact subset is documented in `docs/jbi_gbif_doi_protocol.md`.

**Not verified:** before submission, archive the exact code and data release at a permanent DOI; submit the prepared broad occurrence request through an authenticated GBIF account; register the frozen exact subset as a GBIF Derived Dataset; and insert the resulting release, GBIF download and GBIF Derived Dataset DOI citations. The GitHub repository alone should not be treated as the final preservation record.
"""
if text.count(old) != 1:
    raise SystemExit(f"Expected one Data Accessibility block, found {text.count(old)}")
text = text.replace(old, new, 1)

old_support = """**Table S17** records V.PhyloMaker2 species placement. **Appendix S1** reconstructs the automated literature-search chronology from repository commits and states which human-review details remain unrecorded. The Supporting Information index also identifies the GBIF QC manifests, Open Tree topology, dated S1–S3 trees and both phylogenetic provenance manifests.
"""
new_support = """**Table S17** records V.PhyloMaker2 species placement. **Appendix S1** reconstructs the automated literature-search chronology from repository commits and states which human-review details remain unrecorded. The Supporting Information index also identifies the GBIF QC manifests, the exact occurrence and DOI-preparation bundle, Open Tree topology, dated S1–S3 trees and both phylogenetic provenance manifests.
"""
if text.count(old_support) != 1:
    raise SystemExit(f"Expected one Supporting Information sentence, found {text.count(old_support)}")
text = text.replace(old_support, new_support, 1)

required = [
    "58,455-row occurrence subset",
    "f25ae0cf2c84c45ae461a932d6c6063edda64591913a2495e4a3da82d573f094",
    "389 parent GBIF datasets",
    "docs/jbi_gbif_doi_protocol.md",
    "GBIF Derived Dataset DOI citations",
]
for phrase in required:
    if phrase not in text:
        raise SystemExit(f"Missing required accessibility phrase: {phrase}")

path.write_text(text, encoding="utf-8")
print("GBIF accessibility statement integrated")
