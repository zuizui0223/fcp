#!/usr/bin/env python3
from pathlib import Path

path = Path("docs/jbi_manuscript_editorial_revision.md")
text = path.read_text(encoding="utf-8")

old = """No explicit language filter was coded, but the English search phrases and dependence on indexed title and abstract metadata may have produced language and database-coverage bias. **Not verified:** the original execution dates of all discovery and follow-up searches, the identities and number of manual screeners, and the procedure for resolving disagreements are not recorded in the latest manuscript materials. These details must be supplied before submission.
"""
new = """Repository history records automated global discovery and targeted follow-up or enrichment activity from 16 to 19 July 2026 (Appendix S1). The first global output was committed on 16 July; deferred-candidate follow-up and evidence aggregation occurred on 16–17 July; and targeted enrichment of confirmed cases with ambiguous spatial evidence was committed on 19 July. No explicit language filter was coded, and the English search phrases and dependence on indexed title and abstract metadata may have produced language and database-coverage bias. **Not verified:** the repository does not record the identities or number of human screeners, independent duplicate screening, inter-reviewer agreement or a formal disagreement-resolution procedure. These details require author confirmation and must not be reconstructed retrospectively.
"""
if text.count(old) != 1:
    raise SystemExit(f"Expected one literature-method paragraph, found {text.count(old)}")
text = text.replace(old, new, 1)

old_support = """Tables S1–S7 contain the original model matrices and audits. **Table S8** is the primary baseline model dataset; **Table S9** is the paginated model dataset; **Table S10** is the GBIF taxonomic and retention audit; **Tables S11–S12** report paginated robustness and family deletion; **Tables S13–S15** report Open Tree phylogenetic summaries, name resolution and all 200 replicate fits; **Table S16** reports the six fixed-seed dated-megaphylogeny models; and **Table S17** records V.PhyloMaker2 species placement. The Supporting Information index also identifies the GBIF QC manifests, Open Tree topology, dated S1–S3 trees and both phylogenetic provenance manifests.
"""
new_support = """Tables S1–S7 contain the original model matrices and audits. **Table S8** is the primary baseline model dataset; **Table S9** is the paginated model dataset; **Table S10** is the GBIF taxonomic and retention audit; **Tables S11–S12** report paginated robustness and family deletion; **Tables S13–S15** report Open Tree phylogenetic summaries, name resolution and all 200 replicate fits; **Table S16** reports the six fixed-seed dated-megaphylogeny models; and **Table S17** records V.PhyloMaker2 species placement. **Appendix S1** reconstructs the automated literature-search chronology from repository commits and states which human-review details remain unrecorded. The Supporting Information index also identifies the GBIF QC manifests, Open Tree topology, dated S1–S3 trees and both phylogenetic provenance manifests.
"""
if text.count(old_support) != 1:
    raise SystemExit(f"Expected one Supporting Information paragraph, found {text.count(old_support)}")
text = text.replace(old_support, new_support, 1)

required = [
    "from 16 to 19 July 2026",
    "Appendix S1",
    "independent duplicate screening",
    "must not be reconstructed retrospectively",
]
for phrase in required:
    if phrase not in text:
        raise SystemExit(f"Missing required phrase: {phrase}")

path.write_text(text, encoding="utf-8")
print("literature provenance integrated")
