# FCP 27 candidate photographic modes — blind original-photo botanical verification

Date 2026-10-23. This is **human validation preparation, not automatic confirmation of botanical morph identity**. Work is stacked on FCP nine-hue PR #260 and mode/shape PR #262; do not modify original fixed 447,500 photo source or published hypothesis ledger.

## Evidence and purpose

On the **strict 1,761 species with >=40 original classifiable photos and >=30 source photo pairs within 50 km**, 27 species have reproducible **photo-PC1 density multimodality candidates**, 43 are unstable, and 1,691 did not meet the mode detector gate. The 27 are not known to have two genetic flower-colour morphs, nor is there any confident claim of true three-mode phenotypes. Their stronger geographic similarity is exploratory; it weakens after adjustment for species-wide continuous JSD, photo sample size and near-pair opportunity (mode vs non-mode +0.03263, genus-bootstrap CI [−0.00038,+0.06684]).

The crucial next step is **original image botanical review**. Blind observers to previous mode classification, photo-PC1 scores, location, latitude, photo dates, observers, taxonomic names and candidate/control status. Image content itself can reveal some species; masking cannot ensure perfect taxonomic blindness. Avoid making any claim that prepared blank forms constitute a completed review.

## Matching and review sampling

- Take **all 27 source-verified recurrent photo-mode candidate taxa**.
- Select **27 distinct matched controls from source taxa not resolved as multimodal**, excluding the 43 unstable candidates. Controls must have >=40 original geo-valid classified photos, to avoid using species where the same mode-detection method had inadequate evidence.
- One-to-one Hungarian minimum-cost assignment matches original classified photo count, number of 50-km photograph-pair opportunities, species-wide continuous nine-hue Jensen–Shannon divergence, dominant Hellinger PC1 variance fraction; include penalties for different rough colour-label counts and source cohorts. This is a **QC matched set**, not causal matching, and controls should be described as *not detected*, not genetically monomorphic.
- Per taxon sample 3 images near the 12th, 3 near the 50th and 3 near the 88th percentiles of the **same original unsupervised Hellinger-PC1 score**. Select distinct photos and original observers. Central selections are density-centre photographs, not proven intermediate colour forms.
- If available, sample 2 additional original photos with no accepted four-state colour label, to estimate segmentation, organ and quality failure among missing label cases. No missing colours are imputed.
- Expected total across 54 species: 486 classified plus up to 108 unclassified images.
- Produce anonymous stable bundle IDs, individually shuffled anonymous image IDs and a reviewer-specific random order. Case/control status, exact original photo IDs, source URLs, original photographic coarse labels and PC1 strata are retained *only in the separate PI map*.
- Present image previews if **exact original authenticated source photo URL** is present and retrieval succeeds. If a URL is absent/unreachable, the sheet shows "image not available", with incomplete status recorded. No photograph is invented or replaced to repair source failure. Image previews have EXIF metadata stripped; no exact location or observer identity appears in reviewer packets.
- Keep a source license and photographer attribution ledger in PI mapping, and verify image reuse/license/credit conditions before redistributing any images beyond internal blind reviewers.

## Two independent annotators

Use two blinded independent annotators. For each anonymous photo, record:

1. target flower and *petal region* visible? (yes/no/uncertain)
2. unambiguous single biological colour target? (yes/no/uncertain)
3. dominant *perceived* floral colour (white, cream, yellow, orange, red, pink, magenta, purple, blue, green, other, unassessable)
4. flower stage (bud, freshly open, ageing/senescent, damaged, uncertain)
5. lighting, overexposure, white balance, or shade problematic? (yes/no/uncertain)
6. multiple flowers/parts (including sepals, bracts, stamens) or background might create false hues? (yes/no/uncertain)
7. horticultural/cultivation or taxon identification concerns (yes/no/uncertain)
8. reviewer confidence and notes.

Reviewers must not consult mapping, previous model scores or each other's decisions until both first-pass annotations are frozen. Use separate nonpublic controlled delivery for the reviewer packets; **GitHub Actions artifacts in a public repository are not a secure blinded-secret vault** even when outputs reside in separately named artifacts. The distribution manager must safeguard the PI map and confirm attribution rights. De-identification is methodological, not cryptographic.

## Statistical acceptance after human review (not done in this stage)

Require an agreed target petal/flower identity and stable fully open flower condition for colour-mode claims. Evaluate inter-reviewer agreement; retain "uncertain" rather than forcing disagreement into a colour class. Refit the exact PC1 density and spatial models after separately flagging (i) valid floral organ + good colour conditions and (ii) suspect/invalid images. Report mode persistence and its effect on matched control false-discovery rates, but do **not** relabel taxon as natural genetic polymorphism until independent individual/offspring/genotype or repeated field evidence is available.

The historical strict <=5km same-year/month separate-observer photo null remains unsupported and must remain visible.

## Source and reproducibility

- Source photos: immutable release fcp-all4730-original-method-actual-colour-results-20261017, original photo CSV SHA256 8e7f10f0427927208485b2da6fa9c46bb681d6f3611e53b8dfdc606cfc8dabe6
- Mode source CI #39923312111: original 1761 species, 27 recurrent candidates, 43 unstable, 1691 nonresolved.
- Continuous geographic CI #39918785955: 1761 taxa and original 50km reference.
- Code: scripts/analysis/build_fcp27_original_blind_photo_review_20261023.py
- CI: .github/workflows/fcp27-original-photo-blind-review-20261023.yml
- Output: 2 reviewer annotation sheets and image/contact-sheet indexes, a PI-only mapping, source/photo retrieval QC and manifest.

No human review, genotyping, pigment chemistry, experimental fitness or universal adaptation result is claimed by creation of these files.
