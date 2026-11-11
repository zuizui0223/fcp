# FCP original-405 whole-photo QC and 50-km matched control audit (2026-10-10)

**Status:** retrospective photograph-level *technical-quality* sensitivity completed. **0/405 human reannotation cases have been independently scored**. No botanical morph misclassification estimate and no change to the frozen New Phytologist findings.

## Verification and source

The exact source original-photo verification was green 405/405 SHA256 matches, as recorded in `results/polymorphism_adversarial_label_swaps_posthoc_20261010/original405_source_photo_byte_identity_20261010.json` and workflow [38028308416](https://github.com/zuizui0223/fcp/actions/runs/38028308416). The locally analyzed, downloaded immutable review queue artifact ZIP (Actions 38027546572; artifact 11660917466) had SHA256 `a47f82c80273daa4b727930ecb7c13f2957f45e750d30391d9d23f318eb66811`. Full pixel-QC input was the separate 405-image byte-recovery numerical artifact: source photo ID, image SHA256, whole-frame clipping/channel-spread and image byte length. No images or human classifications were redistributed or changed.

## Control-distance diagnosis

There are 140 result-influential source colour photos, 125 originally selected same-species and same-original-colour matched controls, and 140 nonoverlapping originally sampled random controls (405 original photographs, all verified). **Only 99/125 matched control pairs are within 50 km of their index case.** The remaining 26 are >50 km and the farthest is 6,489 km away; 15 index cases have no matched control. For a *local* photographic QC comparison, using all 125 without recognizing this distance failure would conflate locally differing photographic conditions. Distances: <=10km 21; <=25km 71; <=50km **99**; <=100km 101; <=250km 109; <=500km 115; <=1000km 118.

## Actual within-50-km paired numerical result

With 99 matched original same-species/same-source-colour pairs in **62 unique species**, the estimator is the *species-equal average of within-species mean paired differences* (influential photo minus matched control). Species-cluster bootstrap 9,999 draws and sign-flip 9,999 draws, deterministic seed 20261010. Four two-sided exploratory p-values were Holm-adjusted as one family.

| Whole-photo technical indicator | High-impact minus same-species, same-original-colour nearby control | species-bootstrap CI95 | Holm p |
|---|---:|---:|---:|
| fraction of bright-clipped photo pixels | −0.000085 | [−0.001350,+0.001228] | 0.8964 |
| fraction of dark-clipped photo pixels | −0.000528 | [−0.001460,+0.000441] | 0.5768 |
| mean whole-photo RGB channel spread | −0.017555 | [−0.037117,+0.002072] | 0.3868 |
| log original image byte size | −0.113728 | [−0.264722,+0.035407] | 0.4548 |

No tested **whole-photo technical feature** had a significant effect under this exploratory procedure, but that cannot rule out flower-specific misclassification, inaccurate floral organ identification, or geographic errors. The four features were not computed inside verified flower petal ROIs. The matched cases are outcome-selected photo cases; they are not a representative botanical error prevalence study.

The composition-preserving adversarial 70 colour-label swaps concentrate on **62 species** (Discovery 24, Validation 20, Third 18) and often target species with only approximately 30–50 original 50-km local pairs. This exposure structure may amplify a small hypothetical number of influential edits; actual independent photo adjudication is still required before treating label vulnerability as observed error.

## Reader protection and next decision

- DO NOT expose the source algorithm colour labels, case kinds or matched-control key to reviewers before both blinded reannotations are frozen. Use only the existing reviewer-facing HTML and blinded 405-ID queue.
- Keep the 99 near-controls, 26 >50km controls, and 15 unmatched index cases separate in any subsequent real-error analysis.
- 405/405 image byte equivalence verifies *identity*, not botanical organ/pigment biology.
- Until paired independent review export files exist, label-error rate, a corrected 50-km colour-depletion effect, genuine population morph frequency, adaptation and selection all remain **unidentified**.
- Source frozen H1/H2, main manuscript, and untouched prospective 2,000+730 species remain unchanged.

Reproducibility: the locally executed Python code, synthetic tests (5/5 passed), frozen source ZIPs and the complete numerical JSON are collected in the companion ChatGPT deliverable `FCP_405画像品質と近隣対照_追加監査_20261010.zip`. This Markdown receipt reports group aggregates only; no sealed key/photo outcome is committed.
