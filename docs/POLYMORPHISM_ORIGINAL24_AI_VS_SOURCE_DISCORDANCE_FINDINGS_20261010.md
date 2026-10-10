# AI-only blinded original-photo pilot: comparison with frozen computer colour labels (2026-10-10)

## New executed evidence; not a human reannotation

The original 24 open-data image thumbnails in this PR passed **24/24 exact prior source-image byte-SHA256 equivalence**. Before unsealing the historical original algorithm colours, a **single AI image-only pass** assigned a coarse colour to 22 photographs and abstained on two. Its image-only CSV was frozen before this comparison, SHA256:

- AI-only blind preliminary photo screen: \`f232f3e710793ae2ae60821638c73f74cd28af3f3d9851bf628046cb2f8189be\`
- Original photo-ID pilot with exact historical image bytes: \`06191b36eca9862f6241b5110ae13aee0ff7e7414a5caaaa1025daffca33eadd\`
- Original separately sealed 405-photo source machine-colour key, recovered from preexisting Actions artifact: \`ef0275fedf7ab4e3e0aa5309fd80c9a9aee97de733b79a93bc2137a04f6f4b5c\`
- Original 405-case *human reviewer-facing* blinded queue: \`8cc44a71dfb64601e23f68b68c25e8cbc534a332e0a17ff51a51bd8737a35c27\`

The sealed colour categories, species, study cohort and case-group membership were consulted only AFTER the AI-only, image-only draft was frozen. Frozen AI-source alignment was verified by original photo IDs, blinded case IDs and contact-sheet order. Five local synthetic/source-hash tests passed and independently checked rejection of tampering with either original label key or blind AI source.

## Exploratory observed AI agreement with previous automated labels

| Item | Exact selected photographs |
|---|---:|
| Total original-source blind pilot | **24** |
| AI assigned one of the four original source colour states | **22** |
| AI assigned source-identical coarse colour | **16** |
| AI and previous pipeline assign different coarse colour | **6** |
| AI abstained because flower organ/colour could not be read from the image | **2** |
| Candidates for independently blinded expert organ-and-colour adjudication (mismatch + abstention) | **8** |

This is **agreement between two non-independent-source measurement methods, not classification accuracy**. No human ground truth is present. The 24 original cases are a source target panel with 9 high-impact photo IDs, 9 same-species/source-colour controls, and 6 random source photos, covering just 23 species.

### Source-composition-dependent *screening priority*, NOT white-flower error frequency

| Previous algorithm's assigned state | n photos | AI agrees | AI disagrees | AI abstains |
|---|---:|---:|---:|---:|
| White | **10** | 4 | 4 | 2 |
| Yellow/orange | 7 | 5 | 2 | 0 |
| Red/pink | 2 | 2 | 0 | 0 |
| Blue/purple | 5 | 5 | 0 | 0 |

The 6 mismatches were concentrated in ambiguous visual-support cases (3), but also occur in two initially marked clear cases and one moderate case. Both AI abstentions had initially been marked unclear. The 8 review triggers were present across all selected design groups (high impact 3/9, matched 3/9, random 2/6). **The difference by previous algorithm state suggests a useful *white-label review priority* within these selected photos, not an empirical white-specific misclassification estimate or a global prevalence statement.**

Possible sources of conflict include a different focal taxon in a mixed flowering photograph, petal vs bract vs background, flower petals versus yellow centres, image scale, lighting and AI's own misreading. Without two independently blinded human botanical experts identifying the target flower organ and floral species (where possible), neither AI nor old source computer classification has adjudicated primacy.

## Boundary protecting future independent human review

- **Never provide human reviewers the source machine-colour key, the AI-only source colour guesses, the identified 8 mismatches, a risk-ranked photo order, or the high-impact/matched/random design group.** Maintain the established original case-ID order or independently randomize all 405 case IDs without using the source colours.
- Keep the publicly committed result strictly to aggregate counts. The row-level comparison CSV and pre-sealed AI image-only source CSV are stored in an internal, access-controlled companion audit, not committed or exposed to human reviewer files.
- For each of two independent human reviewers, collect source-blinded case-ID, whether the focal flowering plant/corolla is actually visible, which anatomical organ is coloured, multiple potentially flowering taxa, colour category or abstain, photo quality, and confidence before unsealing. If target species identity cannot be inferred from the image alone, add a **second stage with target species name disclosed but source machine colour still concealed**, frozen separately to avoid confusing colour and target-taxon identification.
- The eight candidate cases should be examined, but the whole 24 pilot must be reviewed unranked, and the eventual full 405 sample should be treated as selected for this audit rather than a global random survey. Sampling weights and external validation require separate design.

### Explicit inference firewall

The 24 AI-source disagreements must **never** be inserted as true corrected labels into the New Phytologist 50km ecological analysis. The earlier 0.3–0.5% adversarial label-swap scenario was outcome-informed, and the present AI images do not adjudicate which real photos, if any, were misclassified. No case is ground-truth corrected, and no original H1/H2, geographic source photo colours, 405 human evaluation record or unopened future 2,000+730 taxon source cohorts were changed.

[Source-scoped aggregate numerical receipt](../results/polymorphism_original24_ai_blind_source_discordance_20261010/aggregate_receipt.json)

**Reproduction:** The companion standalone Python audit accepts the **previously sealed blind AI image-only source CSV**, original 24 manifest and immutable archival 405 sealed-key ZIP, requires their literal SHA256, verifies ID alignment and case order, and writes a public-eligible aggregate and a separately marked **INTERNAL UNBLINDED DO NOT SHOW REVIEWERS** rowwise verification file. The selected AI CSV and rowwise original label comparisons are intentionally not in this open PR, so the original full 405-blinded human review has not been compromised.


## Completed stage-2 species identity supplement

A separate, SHA-verified local generator, [target-taxon-only stage-2 code](../scripts/analysis/build_fcp_stage2_taxon_only_review_20261010.py), has also been executed against the original archived sealed key. It produced exactly **405 original source case/photo identities with 169 distinct target taxon names**. The derived reviewer-visible columns are strictly:

- audit_case_id
- photo_id
- focal_taxon_name_no_prior_colour

No original algorithm colour, AI preliminary colour, high-impact/matched/random source sampling group, collection location or cohort is included. **Keep this 405-case second-stage file inaccessible to either reviewer until both stage-1 image-only annotations have been completed and frozen.** Otherwise knowledge of the focal species' typical colour may contaminate the truly image-blind first stage. This is a species-aware verification of floral organ identity, not a colour-code adjudication and not a genuine matched plant individual identifier. The reviewed source 405 panel remains untouched and human reannotations remain zero.

The original taxon-only row records have intentionally not been committed to the repository and are supplied in a separately marked internal ZIP for the researcher to release at the correct review stage.
