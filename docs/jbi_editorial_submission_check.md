# Journal of Biogeography editorial submission check

## Verified source

- Repository: `zuizui0223/fcp`
- Branch reviewed: `main`
- Latest manuscript reviewed: `docs/jbi_manuscript_draft.md`
- Reviewed manuscript blob SHA: `91d061710613790091b0489e4035246bbdffd6ee`
- Editorial revision: `docs/jbi_manuscript_editorial_revision.md`

No unverified citation, figure number, table number, Supporting Information number, author detail, funding statement or archive DOI has been invented. Unresolved items are marked `Not verified` in the editorial revision.

## JBI format check

| Item | Status | Finding |
|---|---|---|
| Title | ✓ Pass after revision | Revised to 74 characters and made non-causal. |
| Running title | ✓ Pass after revision | Added; 35 characters. |
| Abstract structure | ✓ Pass | Uses Aim, Location, Taxon, Methods, Results and Main conclusions. |
| Abstract length | ✓ Pass | Below 300 words excluding headings. |
| Keywords | ✓ Pass | Seven keywords. |
| Introduction | △ Needs revision | Fully reorganized, but literature citations remain unverified. |
| Methods | △ Needs revision | WorldClim source, resolution, variables, metrics and model are verified; screening, GBIF cleaning, software versions and permutation details remain incomplete. |
| Results | △ Needs revision | Focal audited values are retained exactly; complete secondary results and the 20-specification matrix are not yet in the manuscript. |
| Discussion | ✓ Pass after revision | Causal and confirmatory overstatement removed. |
| Acknowledgements | ✗ Missing | Not verified. |
| Author contributions | ✗ Missing | Not verified. |
| Conflict of interest | ✗ Missing | Not verified. |
| Data Accessibility | △ Needs revision | GitHub availability is verified; permanent archive DOI is not verified. |
| Tables | ✗ Missing | Final tables and numbering are not verified. |
| Figure legends | ✗ Missing | Final legends and numbering are not verified. |
| References | ✗ Missing | No verified reference list is present. |
| Supporting Information | ✗ Missing | Final numbering and cross-references are not verified. |

## Statistical language

Preferred focal wording:

> Geographically structured variation was negatively associated with realized moisture niche breadth (odds ratio = 0.426, family-clustered 95% confidence interval = 0.184–0.985; clustered Wald p = 0.0460), whereas permutation support was borderline (two-sided p = 0.0556).

The manuscript should not describe this as a confirmed, robust or causal effect because:

1. five metrics were evaluated at four thresholds;
2. the clustered interval narrowly excludes one;
3. the permutation p-value is 0.0556;
4. the estimate weakens in the broader evidence set; and
5. family deletion is not a phylogenetic correction.

Use `was associated with`, `suggestive`, `evidence-sensitive`, `consistent with`, `may reflect` and `support was borderline`. Avoid `demonstrates`, `confirms`, `climate drives`, `local adaptation explains` and unqualified `robust evidence`.

## Editor Check

**Provisional decision: Major revision.**

The question is suitable for Journal of Biogeography and the distinction between local coexistence and geographic differentiation is potentially publishable. The source-level classification audit is a substantive strength. The manuscript is not submission-ready because the literature foundation, reproducible screening methods, complete GBIF methods, software and permutation details, full multiplicity results, figures, tables, references, declarations and persistent data archive are incomplete.

## Likely reviewer concerns, in priority order

1. Post-analysis selection of moisture breadth from 20 specifications.
2. Lack of a phylogenetic comparative model.
3. Classification error and under-documentation of local coexistence.
4. Non-random, research-effort-dependent literature sampling.
5. Mismatch between species-level realized niches and morph-level biological interpretation.
6. Borderline and method-dependent statistical support.
7. Incomplete literature-discovery and adjudication methods.
8. GBIF cleaning and the operational definition of occupied climate cells.
9. Small-sample cluster-robust inference with 34 species and 25 families.
10. Interpretation of attenuation in the broader evidence set.
11. Negative-control analyses cannot exclude morph-specific environmental sorting.
12. The novelty claim requires a documented and cited literature review.

## Submission blockers

- Add and verify all citations and the complete reference list.
- Document literature search, screening and classification adjudication.
- Document GBIF retrieval and cleaning exactly.
- Report software versions, permutation algorithm, seed and diagnostics.
- Supply the full 20-specification results.
- Add a phylogenetic model or explicitly justify its absence.
- Finalize tables, figures, legends and cross-references.
- Verify authors, affiliations, acknowledgements, funding, CRediT roles and conflicts.
- Archive the exact code and data release with a persistent DOI or reviewer-private link.
- Prepare the separate non-anonymized title page and required taxon image.