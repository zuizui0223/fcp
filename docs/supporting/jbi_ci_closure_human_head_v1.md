# JBI Chapter 1 CI closure record

Date: 2026-08-30  
Branch: `analysis/jbi-global-colour-boundaries`  
PR: #20

## Purpose

This record returns the pull-request head to a human-authored commit after canonical-figure generation. GitHub does not recursively start ordinary workflow runs from commits created with the workflow `GITHUB_TOKEN`; those bot-originated check suites can therefore appear as `action_required` even when the generating job itself succeeded.

## Completed validation before this closure commit

- the held-out 720-photograph evaluation and frozen continuous-colour representation are complete;
- the Stage-A species-conditioned random-labelling result is supported (`p = 0.0113`);
- the Stage-B primary shared-transition concentration result is not supported (`p = 0.0906`);
- the Chapter-1 spatial test suite and manuscript/result guards pass;
- the legacy frozen 34-species reproduction lane completed its dedicated tests, five models, phylogenetic sensitivities, 3,000-replicate precision analysis, CR2 analysis, numerical regression checks, artifact upload and canonical-figure integration in workflow run `33308881615`;
- the legacy workflow no longer discovers unrelated repository tests;
- canonical 34-species PDF generation now uses a fixed `SOURCE_DATE_EPOCH` and a two-render SHA-256 equality check, preventing creation timestamps from generating false figure diffs.

## Interpretation boundary

This is CI/provenance documentation only. It changes no data, photograph measurement, statistical result, figure content or biological conclusion. The authoritative scientific state remains:

> Continuous flower colour is spatially organized within species in the six-species held-out analysis, while a universal cross-species global transition boundary is not confirmed at the frozen primary support.
