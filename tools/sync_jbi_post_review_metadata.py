#!/usr/bin/env python3
from __future__ import annotations

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MANUSCRIPT = ROOT / "docs/jbi_manuscript_editorial_revision.md"
TITLE_PAGE = ROOT / "docs/jbi_title_page_template.md"
EDITORIAL = ROOT / "docs/jbi_editorial_submission_check.md"
CHECKLIST = ROOT / "docs/jbi_submission_completion_checklist.md"
AUTHOR_FORM = ROOT / "docs/jbi_author_confirmation_form.md"
VALIDATOR = ROOT / "validate_jbi_submission_package.py"


def section(text: str, start: str, end: str) -> str:
    i = text.find(start)
    j = text.find(end, i + len(start))
    if i < 0 or j < 0:
        raise SystemExit(f"Missing section markers: {start!r}, {end!r}")
    return text[i:j]


def word_count(text: str) -> int:
    return len(re.findall(r"\b[\w’'-]+\b", text, flags=re.UNICODE))


def replace_regex(text: str, pattern: str, replacement: str, label: str) -> str:
    updated, count = re.subn(pattern, replacement, text, count=1, flags=re.MULTILINE)
    if count != 1:
        raise SystemExit(f"Expected one {label}; found {count}")
    return updated


def main() -> None:
    manuscript = MANUSCRIPT.read_text(encoding="utf-8")
    abstract = section(manuscript, "## Abstract", "**Keywords:**")
    abstract_words = word_count(re.sub(r"^#{2,3} .+$", "", abstract, flags=re.MULTILINE))
    main_words = word_count(section(manuscript, "## Introduction", "## References"))

    title = TITLE_PAGE.read_text(encoding="utf-8")
    title = replace_regex(
        title,
        r"^- Structured abstract: \*\*\d+ words\*\* excluding headings and keywords\.$",
        f"- Structured abstract: **{abstract_words} words** excluding headings and keywords.",
        "title-page abstract count",
    )
    title = replace_regex(
        title,
        r"^- Main text from Introduction through Discussion: \*\*[\d,]+ words\*\*\.$",
        f"- Main text from Introduction through Discussion: **{main_words:,} words**.",
        "title-page main-text count",
    )
    title = title.replace("included in the audited baseline dataset", "included in the baseline-unambiguous dataset")
    TITLE_PAGE.write_text(title, encoding="utf-8")

    editorial = EDITORIAL.read_text(encoding="utf-8")
    editorial = re.sub(
        r"validated length is \d+ words",
        f"validated length is {abstract_words} words",
        editorial,
    )
    editorial = re.sub(
        r"Introduction–Discussion contains [\d,]+ words",
        f"Introduction–Discussion contains {main_words:,} words",
        editorial,
    )
    editorial = editorial.replace("Tables S1–S17", "Tables S1–S19")
    editorial = editorial.replace("S1–S17 row counts", "S1–S19 row counts")
    editorial = editorial.replace(
        "| Manual screening documentation | △ Needs author confirmation | The repository does not identify human screener names or number, independent duplicate screening, agreement statistics or a formal disagreement-resolution procedure. |",
        "| Classification review state | △ Rule-derived; human review pending | Labels are source-traceable and rule-derived; the resolved queue remains `unreviewed`. S18 is a blank blinded 34-species review sheet, S19 is the separate key, and the evaluator computes agreement without changing analysis inputs. |",
    )
    editorial = editorial.replace(
        "| Classification | ✓ Pass | Four classes, source audit, freeze rule and strict manifest are documented. |",
        "| Classification transparency | ✓ Pass | Four operational classes, deterministic rules, the frozen manifest, 34/34 rule reproduction, source matching, plural-expression audit and the human-review boundary are documented. |",
    )
    editorial = editorial.replace(
        "| Supporting provenance | ✓ Pass | Analysis tables, literature chronology, placement audit, Open Tree topology, three dated trees, exact GBIF bundle and manifests are indexed or directly referenced. |",
        "| Supporting provenance | ✓ Pass | Analysis tables, literature chronology, S18/S19 review packet, classification protocol and rule audit, placement audit, phylogenies, exact GBIF bundle and manifests are indexed or directly referenced. |",
    )
    editorial = editorial.replace(
        "The repository does not establish whether screening was single-reviewer, duplicated or formally adjudicated.",
        "The repository establishes that the current labels are deterministic screening labels and that the resolved queue remains unreviewed; independent human verification must still be completed or explicitly reported as not performed.",
    )
    EDITORIAL.write_text(editorial, encoding="utf-8")

    checklist = CHECKLIST.read_text(encoding="utf-8")
    checklist = checklist.replace("Supporting Tables S1–S17", "Supporting Tables S1–S19")
    checklist = checklist.replace("Tables S1–S17", "Tables S1–S19")
    checklist = checklist.replace(
        "| Human screener names and number of screeners | ✗ Author confirmation required | Complete Section 1 of `docs/jbi_author_confirmation_form.md`. |",
        "| Repository-supported classification state | ✓ Verified | Current labels are deterministic, source-traceable screening labels; the resolved queue remains `unreviewed`. |",
    )
    checklist = checklist.replace(
        "| Independent duplicate screening or agreement assessment | ✗ Author confirmation required | State what was actually done; explicitly say when it was not performed. |",
        "| Blinded classification review packet | ✓ Prepared | S18 hides the current label; S19 is the separate key. `evaluate_jbi_completed_classification_review.py` computes agreement and Cohen's kappa without changing analysis inputs. |",
    )
    checklist = checklist.replace(
        "| Disagreement-resolution procedure | ✗ Author confirmation required | Describe the actual procedure; do not create a retrospective protocol. |",
        "| Completed human verification and adjudication | ✗ Author action required | Complete one or preferably two independent S18 copies, freeze them before opening S19, evaluate agreement and document adjudication. |",
    )
    checklist = checklist.replace(
        "| Supporting Information index | ✓ Prepared and SHA-256 checked |\n| Tables S1–S19 | ✓ Prepared |",
        "| Supporting Information index | ✓ Prepared and SHA-256 checked |\n| Tables S1–S19 | ✓ Prepared |\n| S18/S19 classification-review packet | ✓ Prepared; not completed | S18 is blank and blinded; S19 reproduces 34/34 frozen labels from the retained queue text and all 34 classification sources match the queue best source. |",
    )
    CHECKLIST.write_text(checklist, encoding="utf-8")

    author = AUTHOR_FORM.read_text(encoding="utf-8")
    pattern = re.compile(r"(?ms)^## 1\. Human evidence screening and classification\n.*?(?=^## 2\. Final authorship\n)")
    if len(pattern.findall(author)) != 1:
        raise SystemExit("Expected one human evidence review section in author form")
    review_section = """## 1. Human evidence screening and classification

### Repository-supported state

- Current labels were generated by deterministic text-screening rules.
- Mixed and unclear rows were routed to a manual-review queue and excluded from binary models.
- The resolved queue retains `review_status = unreviewed` for all rows.
- Completed independent human screening, duplicate review, agreement statistics and formal adjudication are not documented in the repository.
- S18 provides a blank 34-species blinded review sheet; S19 contains the separate rule-label key.
- Current code reproduces all 34 frozen labels from the retained queue text, and all 34 classification source identifiers match the queue best source.
- The `polymorphic populations` plural-rule correction changed zero frozen labels and zero model inputs.

### Confirm any review performed outside the repository

- Was any full-source human verification completed outside the repository? [Yes/No]
- Reviewer name(s) or initials: [Name(s) / Not applicable]
- Species or stages reviewed: [Details / Not applicable]
- Independent duplicate review: [Performed / Not performed / Partly performed]
- Inter-reviewer agreement calculated: [Yes/No]
- Disagreement resolution: [Actual procedure / Not applicable]
- Evidence record or dated notes available: [Location / None]

Do not describe off-repository review as completed unless a reviewer, scope and date can be documented.

### Recommended completion route

1. Copy S18 once for each reviewer.
2. Keep S19 closed until first-pass labels are frozen.
3. Complete all 34 labels, reviewer identity or initials, date and concise notes.
4. Run `evaluate_jbi_completed_classification_review.py` on the completed sheet(s).
5. Adjudicate reviewer disagreements and consensus disagreements with the frozen rule label.
6. Record accepted changes in the correction log and rerun every analysis; do not overwrite S6 or S8 automatically.

### Manuscript wording decision

Select the factually correct final statement:

- [ ] `Spatial labels were generated by deterministic screening rules; completed independent human verification was not performed.`
- [ ] `One investigator verified all 34 classifications against the cited sources; independent duplicate review and inter-reviewer agreement assessment were not performed.`
- [ ] `Two or more investigators independently verified all 34 classifications; raw agreement and Cohen's kappa were [values], and disagreements were resolved by [procedure].`
- [ ] Other verified wording: [text]

"""
    author = pattern.sub(review_section, author, count=1)
    AUTHOR_FORM.write_text(author, encoding="utf-8")

    validator = VALIDATOR.read_text(encoding="utf-8")
    insertion = """
    require(
        f"Structured abstract: **{abstract_words} words**" in title_page,
        f"Title-page abstract count is stale; expected {abstract_words}",
        failures,
    )
    require(
        f"Main text from Introduction through Discussion: **{manuscript_words:,} words**" in title_page,
        f"Title-page main-text count is stale; expected {manuscript_words:,}",
        failures,
    )
    require(
        f"validated length is {abstract_words} words" in editorial,
        f"Editorial abstract count is stale; expected {abstract_words}",
        failures,
    )
    require(
        f"Introduction–Discussion contains {manuscript_words:,} words" in editorial,
        f"Editorial main-text count is stale; expected {manuscript_words:,}",
        failures,
    )
"""
    marker = "    if manuscript_words > 6000:\n        warnings.append(f\"Introduction–Discussion word count exceeds 6,000: {manuscript_words}\")\n"
    if "Title-page abstract count is stale" not in validator:
        if marker not in validator:
            raise SystemExit("Validator word-count insertion marker missing")
        validator = validator.replace(marker, marker + insertion, 1)
    VALIDATOR.write_text(validator, encoding="utf-8")

    print({"abstract_words": abstract_words, "main_words": main_words})


if __name__ == "__main__":
    main()
