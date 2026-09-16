# H2 Third-Cohort Preopening Frame Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Freeze an outcome-blind third-cohort candidate universe that excludes all legacy high-depth species and every species selected into P500, without opening any flower-colour outcome.

**Architecture:** Reuse the immutable prospective P500 selection artifact as the parent lineage. The parent P100 file already represents U100 after exclusion of all 1,000 legacy high-depth species; derive the third-cohort universe by subtracting all 500 prospectively selected P500 species, including the single species that later failed to reach 100 fresh metadata rows. Fail closed on source hashes, column allowlists, U100 depth, identity uniqueness, and overlap.

**Tech Stack:** Python 3 standard library, pytest, CSV/JSON/SHA256.

**Spec:** `docs/POLYMORPHISM_H2_THIRD_COHORT_PREOPENING_FREEZE_20260916.md`

## Global Constraints

- Biological colour outcomes, morph, D, palette coordinates, H2 W/p/verdict, H3 predictors, and the P500 measured table remain unopened.
- Parent selection artifact is workflow `34708044962`, artifact `10302477571`, digest `sha256:7e43763e3f1d9fa3ff108154dc6ef9443430b70fdb62532e2f6845582d607abb`.
- P100 SHA256 is `1473aad680fe2fa84903c5e11ee104fd0eceef828957f16c2ef1a35f9dd6993c`.
- P500 SHA256 is `f54d07fb2338a20a3f92808b58f0ebc948ab0d5af3cb01fb26702f861184b7a4`.
- Chain must remain U0 42,111 -> U100 4,730 -> legacy-removed P100 3,730 -> subtract all P500 selected 500 -> third candidate universe 3,230.
- This task freezes the candidate universe only; it does not select the final third-cohort sample and does not authorize biological opening.

---

### Task 1: Candidate-frame derivation

**Files:**
- Create: `tests/test_build_polymorphism_h2_third_cohort_frame_20260916.py`
- Create: `scripts/analysis/build_polymorphism_h2_third_cohort_frame_20260916.py`
- Create: `results/polymorphism_h2_third_cohort_candidate_frame_20260916/candidate_pool.csv`
- Create: `results/polymorphism_h2_third_cohort_candidate_frame_20260916/result.json`

**Interfaces:**
- Consumes: unpacked `p100_outcome_blind_pool.csv` and `p500_frozen_selection.csv` from parent artifact 10302477571.
- Produces: three-column outcome-blind candidate pool and machine-readable lineage/firewall receipt.

- [x] **Step 1: Write failing tests** for removal of all P500-selected species, deterministic sorting, source-column firewall, and U100 minimum capacity.
- [x] **Step 2: Run tests and verify RED** because the builder does not yet exist.
- [x] **Step 3: Implement the minimal builder** with hash/count fingerprints and fail-closed validation.
- [x] **Step 4: Run tests and verify GREEN**.
- [x] **Step 5: Execute against the immutable parent artifact** and verify 3,230 rows, zero P500 overlap, minimum capacity 100, and candidate SHA256 `7fc0337074b55a91c0b50773a8d5c3ef82e877074c8b3cacc21f9af58e0ba77e`.

### Task 2: Freeze chronology

**Files:**
- Modify: `docs/POLYMORPHISM_H2_THIRD_COHORT_PREOPENING_FREEZE_20260916.md`

**Interfaces:**
- Consumes: Task 1 result receipt.
- Produces: explicit closed candidate-universe gate and next-step boundary.

- [ ] **Step 1:** Record the resolved legacy exclusion discrepancy: the actual prospective selection receipt used 500 discovery + 500 reserve species, not only the 369/363 D-evaluable subsets described in older prose.
- [ ] **Step 2:** Record exclusion of all 500 P500 selected species rather than only the 499 pixel-opened species.
- [ ] **Step 3:** Record the 3,230-species candidate universe and SHA256.
- [ ] **Step 4:** Keep the next gate closed: final third-cohort sample selection and synthetic end-to-end serialization/artifact qualification must occur before biological opening.
