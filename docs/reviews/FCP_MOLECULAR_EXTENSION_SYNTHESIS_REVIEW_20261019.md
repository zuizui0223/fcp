# FCP molecular extension and integrated-claim review — 2026-10-19

## Scope and reviewers

The root integrator independently checked the deposited-expression arithmetic and source/output hashes. The molecular analyst then reviewed the root-written integrated Japanese report, particularly its biological interpretation and distinct evidence denominators. Neither check establishes species orthology, causal alleles, natural-population representativeness, pigment flux or fitness.

## Independent deposited-data calculation

Read the original gzip/TSV bytes with the Python CSV parser, preserving duplicate headers and column order. Checked every derived source-column index, literal header and occurrence index against the unmodified table, then independently calculated the reported ratios from the extracted source values.

| Check | Result |
|---|---:|
| Per-column FPKM values checked | 36 |
| Maximum difference from source FPKM | 0 |
| Gene-by-colour log2 ratios checked | 12 |
| Maximum difference in log2 ratio | 4.8634 × 10⁻¹² |
| Ambiguous Matthiola white value rows retaining empty sample accessions | 6 |
| Distinct drought sample accessions | 18 |
| Drought phenotype-by-treatment cells | 6, each with 3 submitted samples |
| Molecular input hashes checked | 82 |
| Molecular output hashes checked | 17 |

The small ratio difference is consistent with the derived CSV's twelve-significant-digit serialization. The builder's current hash also matches its final receipt. The focused molecular tests cover five checks; the root's combined new-analysis and dependency suite contains 37 passing tests.

## Biological claims and denominators

- Original, expanded and union title-screen counts are distinct: 131/15, 682/65 and 801/78 runs/projects. The six environmental species in the screen union are not the seven environmental species with curated T2+ evidence.
- Combined curation is 14 taxa, 24 contrast records and 39 contrast-by-organ rows. Shared alleles, cultivar derivatives, repeated studies and tissues are not counted as independent origins.
- Orychophragmus RNA material is **white-violet versus violet**, with source-reported two biological replicates per group, not pure-white knockout material. The integrated text retains this distinction.
- The expression calculations use deposited Arabidopsis-reference assignments and disclose the inconsistent original reference-processing description. No orthology, count-based differential-expression, flux, resource-saving or adaptive claim is made from those ratios.
- Matthiola's repeated W1 column names remain unresolved. The conditional descriptive means do not assert two verified independent white plants or identify the alleles from the separate 2021 cultivar study.
- Antirrhinum's Eluta evidence concerns near-white regions with retained coloured patterning. The RNA genotype/year confounding is retained without discarding the independent gene-intervention evidence.
- The drought dataset's independent plant/pool unit and pre/post-treatment morph definition remain unknown. The six-cell sample layout does not establish six randomized biological treatment groups or a demonstrated genotype-by-environment adaptation. Submitted abstract effect sizes are not promoted to results.
- Species-level literature annotation is allowed without matching the photographed individual. The lack of a representative mechanism-by-environment regression reflects the sparse, heterogeneous, selectively documented seven-species set.

## Integrated prose correction

The reviewer flagged “distant places” as too broad for the unbuffered spatial-cell comparisons: adjacent cells need not be distant. The report opening and subsection heading now say **held-out spatial cells, including the 50 km buffer conditions**, matching the actual design. The B and C comparisons explicitly retain their different model structures, training units and loss definitions.

The draft's molecular and tradeoff claims otherwise matched the available evidence. Temperature association, pigment-production costs, floral heat/water balance and long-term polymorphism maintenance remain distinct claims with distinct evidence requirements.

## Verification commands

The root directly read the original matrices, derived column/ratio CSVs, all input/output receipt entries and drought sample map, using `csv`, `gzip`, `hashlib`, `math` and `pandas`. The executable extraction and its tests are in `scripts/analysis/build_fcp_molecular_expansion_20261019.py` and `tests/test_fcp_molecular_expansion_20261019.py`.

The complete executed analysis and scoped limitations are in [the integrated report](../FCP_COMMON_UNITS_TRANSPORT_MOLECULAR_20261019.md). This review is a scientific-data and claim check, not a full repository test or a new GitHub Actions result.
