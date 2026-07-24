#!/usr/bin/env python3
"""Validate the limited, source-backed author metadata prefill.

This validator protects the verified name, division, institution, general postal
address and institutional email while ensuring that the prefill is not silently
converted into an unverified final author order, corresponding-author designation,
laboratory assignment or ORCID.
"""
from __future__ import annotations

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent
DOCS = ROOT / "docs"
TITLE_PAGE = DOCS / "jbi_title_page_template.md"
AUTHOR_FORM = DOCS / "jbi_author_confirmation_form.md"
COVER_LETTER = DOCS / "jbi_cover_letter_template.md"
CHECKLIST = DOCS / "jbi_submission_completion_checklist.md"

VERIFIED_NAME = "ZHANG RUIQI"
VERIFIED_DIVISION = "Division of Forest and Biomaterials Science"
VERIFIED_AFFILIATION = "Graduate School of Agriculture, Kyoto University"
VERIFIED_POSTAL = "Kitashirakawa Oiwake-cho, Sakyo-ku, Kyoto 606-8502"
VERIFIED_EMAIL = "zhang.ruiqi.77h@st.kyoto-u.ac.jp"


def fail(message: str) -> None:
    raise SystemExit(message)


def require_text(text: str, phrase: str, label: str) -> None:
    if phrase not in text:
        fail(f"{label} missing verified phrase: {phrase}")


def main() -> None:
    paths = (TITLE_PAGE, AUTHOR_FORM, COVER_LETTER, CHECKLIST)
    for path in paths:
        if not path.is_file() or path.stat().st_size == 0:
            fail(f"Missing or empty author-prefill file: {path}")

    title_page = TITLE_PAGE.read_text(encoding="utf-8")
    author_form = AUTHOR_FORM.read_text(encoding="utf-8")
    cover_letter = COVER_LETTER.read_text(encoding="utf-8")
    checklist = CHECKLIST.read_text(encoding="utf-8")

    for label, text in (
        ("title page", title_page),
        ("author form", author_form),
        ("cover letter", cover_letter),
        ("completion checklist", checklist),
    ):
        require_text(text, VERIFIED_NAME, label)
        require_text(text, VERIFIED_DIVISION, label)
        require_text(text, VERIFIED_AFFILIATION, label)
        require_text(text, VERIFIED_POSTAL, label)
        require_text(text, VERIFIED_EMAIL, label)

    require_text(title_page, "final author order", "title page")
    require_text(title_page, "any additional authors", "title page")
    require_text(title_page, "laboratory", "title page")
    require_text(title_page, "Corresponding author", "title page")
    require_text(author_form, "[Confirm order]", "author form")
    require_text(author_form, "Not verified", "author form")
    require_text(author_form, "research student", "author form")
    require_text(author_form, "No ORCID linked to the institutional email was found", "author form")
    require_text(cover_letter, "CONFIRM THAT THE FOLLOWING VERIFIED PERSON IS THE CORRESPONDING AUTHOR", "cover letter")
    require_text(checklist, "Verified identity prefill", "completion checklist")
    require_text(checklist, "Verified division and affiliation prefill", "completion checklist")
    require_text(checklist, "Verified institutional postal address", "completion checklist")
    require_text(checklist, "Verified institutional email prefill", "completion checklist")
    require_text(checklist, "Final author order and any additional authors", "completion checklist")

    combined = "\n".join((title_page, author_form, cover_letter))
    orcids = re.findall(r"\b\d{4}-\d{4}-\d{4}-\d{3}[\dX]\b", combined)
    emails = re.findall(r"[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\.[A-Za-z]{2,}", combined)
    unexpected_emails = sorted({email for email in emails if email.lower() != VERIFIED_EMAIL.lower()})
    if orcids:
        fail(f"Unverified ORCID-like value found in author-prefill files: {orcids}")
    if unexpected_emails:
        fail(f"Unexpected email-like value found in author-prefill files: {unexpected_emails}")
    if not emails:
        fail("Verified institutional email is absent from author-prefill files")

    print(
        {
            "status": "pass",
            "verified_name": VERIFIED_NAME,
            "verified_division": VERIFIED_DIVISION,
            "verified_affiliation": VERIFIED_AFFILIATION,
            "verified_postal": VERIFIED_POSTAL,
            "verified_email": VERIFIED_EMAIL,
            "laboratory_confirmed": False,
            "final_author_order_confirmed": False,
            "corresponding_author_confirmed": False,
            "email_prefilled": True,
            "orcid_prefilled": False,
        }
    )


if __name__ == "__main__":
    main()
