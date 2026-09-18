# P500 prospective terminal postmortem

Date: 2026-09-16 JST

This document records the terminal state of the one authorized prospective P500 flower-colour measurement/H2 execution. It is a postmortem and claim-ceiling record. It does **not** authorize a rerun, replacement, rescue analysis, threshold change, axis refit, or retrospective promotion of a recomputation to prospective evidence.

## 1. Frozen execution identity

- repository: `zuizui0223/fcp`
- execution branch: `analysis/p500-prospective-execution-gate-20260915`
- exact execution head: `3bf8692db814ae5e8ed81bdfa11b2af2e6042960`
- PREOPENING run: `34952750200` — success
- prospective location-blind measurement run: `34953307374` — terminal workflow conclusion `failure`
- authorization: `docs/POLYMORPHISM_H2_P500_PROSPECTIVE_MEASUREMENT_AUTHORIZATION_20260915.json`
- authorization scope: one bounded run only; no post-outcome rerun, replacement, rescue, or threshold change

The workflow-level `failure` must not be interpreted as an H2 biological negative. The failure occurred after the complete measurement/support stage had passed and after the H2 calculator had run, but before the H2 result was durably serialized and uploaded.

## 2. Measurement and support stage: closed PASS

The location-blind measurement and reassembly completed successfully under the frozen design.

- terminal partitions: **256 / 256**
- frozen rows: **49,900 / 49,900**
- frozen species: **499 / 499**
- unique measurement IDs: **49,900**
- duplicate measurement IDs: **0**
- classifiable rows: **25,162**
- nonclassifiable rows: **24,738**
- measurement-evaluable species under frozen `n_classifiable >= 40`: **373**
- minimum measurement-evaluable species: **250**
- support decision: **PASS**
- replacement rows: **0**
- replacement species: **0**
- persisted image pixels: **false**

Sealed measurement-result artifact:

- name: `p500-prospective-measurement-result-v1`
- artifact id: `10412947874`
- digest: `sha256:0c769f2d657f3201c253cb4a537ee01c998bd0a780f49066cf0ced1e1e9d6872`

The `reassemble-and-support-gate` job (`104517747363`) completed successfully. Thus P500 is not `NOT_EVALUABLE` because of inadequate photo support. The prospective exercise demonstrates that the frozen location-blind measurement pipeline can carry a much larger independent high-depth cohort through complete acquisition, colour measurement, reassembly, and the predeclared support gate.

This is a measurement-feasibility / transport result. It is not by itself evidence that the white-versus-nonwhite H2 geometry transports.

## 3. H2 stage: calculation executed, durable outcome lost

The downstream `h2-prospective-test` job (`104518017021`) reached the frozen H2 executor after support PASS.

Inspection of the exact frozen script and job log establishes the execution order. Before the exception, the executor had already:

1. constructed the primary 0.10 and strict 0.20 species-delta vector tables;
2. run the frozen structured-null calculations with the predeclared seeds and 999 replicates;
3. assigned primary and strict threshold result objects;
4. assigned the primary H2 verdict in memory;
5. constructed the `H2_COMPLETE` stage receipt; and
6. successfully called `validate_stage_transition` for the H2 stage.

It then failed while constructing the final JSON result object at:

```python
"measurement_result": str(args.measurement_result.relative_to(ROOT)),
```

The workflow passed `args.measurement_result` as a repository-relative `Path`, while `ROOT` was absolute. `pathlib.Path.relative_to()` therefore raised:

```text
ValueError: 'results/.../measurement_result.json' is not in the subpath of '/home/runner/work/fcp/fcp' OR one path is relative and the other is absolute.
```

The subsequent terminal validation and artifact-upload steps were skipped. Querying the run for `p500-prospective-h2-result-v1` returns no artifact.

Therefore the primary H2 vector N, observed W, structured-null p, and CONFIRMED/NOT_CONFIRMED decision were **not durably recorded** in the prospective execution record.

## 4. Scientific terminal classification

The admissible terminal interpretation is:

> **P500 measurement/support PASS; prospective H2 outcome not durably evaluable because of a post-calculation serialization failure.**

This is neither `H2_PROSPECTIVE_WHITE_AXIS_CONFIRMED` nor `H2_PROSPECTIVE_WHITE_AXIS_NOT_CONFIRMED`. It is also not a measurement-support failure.

A deterministic recomputation from the sealed measured table could in principle regenerate numbers, but doing so after observing that the original run failed would be a recovery/reanalysis. Under the frozen one-bounded-run contract, such a recomputation cannot be promoted to the untouched prospective confirmation that P500 was designed to provide. Accordingly, no P500 H2 N/W/p/verdict is recovered here.

## 5. Claim-ceiling consequence

P500 changes one part of the evidence state and leaves another unchanged.

### Raised / newly demonstrated

The response-blind expansion successfully transported the **measurement pipeline and support rule** beyond the original 369/363 high-depth cohorts: 499 frozen species and 49,900 frozen photo rows were processed without replacement, and 373 species passed the predeclared per-species measurement-support requirement.

This strengthens the operational claim that the high-depth photo measurement design is practically deployable at larger scale. It does not estimate global polymorphism prevalence and does not establish that every species passing photo support has biologically resolved polymorphism.

### Not raised

The untouched prospective claim for transport of the fixed white-versus-nonwhite axis remains unresolved because no sealed H2 terminal result exists. The biological H2 evidence therefore remains the existing high-depth discovery/reserve evidence:

- primary 10%: discovery N=152, W=0.514625, p=0.001; reserve N=129, W=0.514586, p=0.001;
- strict 20%: discovery N=75, W=0.542355, p=0.001; reserve N=65, W=0.510517, p=0.008.

Those results support the current targeted white-versus-nonwhite geometry claim, with the existing status boundary that `q_white` was isolated after the initial H2 geometry had been opened.

## 6. Manuscript wording after P500

The strongest defensible paper-level statement remains centered on measurement plus geometry:

> Species-level flower-colour polymorphism can be measured reproducibly from high-depth citizen-science photographs. A separate 499-species, 49,900-photo prospective expansion successfully passed the frozen location-blind measurement and support pipeline without replacement, although a post-calculation serialization failure prevented the prospective H2 white-axis outcome from being durably recorded. In the two original species-disjoint high-depth cohorts, the dominant recurrent geometry of within-species colour variation lies along a white-versus-nonwhite axis rather than a general hue axis; broad phylogenetic signal and the discovery association with sampled geographic span do not replicate.

The prospective P500 execution should appear as a transparent validation/postmortem result, not as confirmatory biological evidence for or against H2.

## 7. Future-only repair

A postmortem branch, `analysis/p500-postmortem-no-rerun-20260916`, was created from the exact execution SHA. On that branch only, the serializer was made robust to both relative and absolute CLI paths and regression tests were added. This repair is for future studies only and does not reopen run `34953307374`.

Any future untouched confirmation of the fixed white/nonwhite target must use a genuinely new, independently authorized cohort/study with a new pre-opening chronology. Replaying P500, replacing its species, rerunning only H2 on its sealed measurements, or relabeling a recovered result as prospective confirmation is outside the frozen authorization.