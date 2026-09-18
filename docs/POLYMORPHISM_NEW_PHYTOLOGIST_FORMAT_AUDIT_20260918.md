# New Phytologist format audit — 2026-09-18

Official source checked on 2026-09-18:

- New Phytologist Author Guidelines: https://nph.onlinelibrary.wiley.com/hub/journal/14698137/about/author-guidelines

This is an editorial-format audit only. It changes no scientific result, cohort, threshold, null model or claim boundary.

## Current official requirements relevant to this manuscript

For a Full Paper, the current author guidelines state that:

- the title should be concise, approximately 130 characters including spaces;
- the manuscript should contain Summary, Introduction, Materials and Methods, Results, Discussion, Acknowledgements, Competing interests, Author contributions, Data availability, References and Supporting Information;
- the Summary should be <=200 words and organized in four bullet points;
- 5–8 keywords should be supplied;
- Full Papers are usually in the region of 6,500–7,500 words;
- Full Papers usually contain 6–8 display items;
- line and page numbering should be used for review, with continuous line numbering;
- full author names, affiliations, corresponding-author email and ORCIDs are required;
- the cover letter should answer three editor questions, with a maximum of 50 words per answer:
  1. What hypotheses or questions does this work address?
  2. How does this work advance our current understanding of plant science?
  3. Why is this work important and timely?

Initial submission is otherwise free-format.

## Audit of the current derivative

File:

- `docs/POLYMORPHISM_MANUSCRIPT_NEW_PHYTOLOGIST.md`

Current measured status:

| Requirement | Current state | Audit |
|---|---:|---|
| Title length | 111 characters | PASS |
| Summary | 149 words | PASS |
| Summary bullets | 4 | PASS |
| Keywords | 6 | PASS |
| Required scientific sections | present | PASS |
| Figures | 5 | |
| Tables | 1 | |
| Total display items | 6 | PASS |
| Introduction | 499 words | |
| Materials and Methods | 1,675 words | |
| Results | 864 words | |
| Discussion | 798 words | |
| Main text | 3,836 words | EDITORIAL RISK: shorter than usual Full Paper range |
| Discussion share of main text | ~25% | PASS relative to <=30% style preference |
| Continuous line numbering | not yet applied | PENDING FILE CONVERSION |
| Author list / affiliations / correspondence | placeholders | EXTERNAL INPUT REQUIRED |
| ORCIDs | not yet inserted | EXTERNAL INPUT REQUIRED |
| Permanent archive DOI/version | placeholder | PENDING ARCHIVE |
| SI legends in final journal form | not yet complete | PENDING |

## Cover-letter audit

File:

- `docs/POLYMORPHISM_NEW_PHYTOLOGIST_COVER_LETTER.md`

The three required editor questions are now explicit.

Current answer lengths:

- Question 1: 36 words;
- Question 2: 36 words;
- Question 3: 40 words.

All are within the 50-word maximum.

Administrative header fields remain placeholders and therefore the cover letter is scientifically formatted but not yet administratively final.

## Main-text length decision

The current 3,836-word main text is substantially shorter than the journal's usual 6,500–7,500-word Full Paper range.

This is not treated as a reason to add generic prose. The defensible route is:

1. retain the concise conceptual Introduction;
2. expand Materials and Methods only where the frozen repository contains replication-critical details currently omitted from the journal derivative;
3. expand Results only where display-item interpretation needs explicit numerical context;
4. expand Discussion only where evidence boundaries, source dependence and biological alternatives need substantive treatment;
5. do not create new H3 predictor analyses or new mechanistic claims to increase length.

A presubmission enquiry is useful if the manuscript remains substantially shorter than the usual Full Paper range after replication-critical Methods are restored.

## Claim boundary

Journal adaptation must not alter the frozen scientific ceiling in:

- `docs/POLYMORPHISM_CURRENT_CLAIM_LEDGER_20260918.md`
- `docs/POLYMORPHISM_H2_THIRD_COHORT_RESULT_AND_MANUSCRIPT_CLAIM_FREEZE_20260917.md`

In particular, no format change can turn species-disjoint same-source confirmation into independent-source replication, a global prevalence estimate, or a mechanistic/adaptive claim.
