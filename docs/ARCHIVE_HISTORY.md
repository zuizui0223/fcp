# Historical / exploratory work removed from the active paper path

The repository began as a mathematical phase-theory project and later accumulated several exploratory empirical branches. The current Journal of Biogeography paper no longer depends on those implementations.

Git history and closed pull requests remain the authoritative archive. In particular:

- early one-population, two-patch, network, temporal and finite-population phase-theory scripts/docs are historical conceptual work, not inputs to the 34-species paper;
- candidate/control experiments and coarse range/turnover analyses were exploratory reviewer checks and are not the primary paper pipeline;
- unreviewed expanded species analyses (including the 51/78/107-species development tracks) are not used as primary evidence;
- PR #8 contains the later systematic-search expansion and remains relevant as search-provenance infrastructure;
- closed PRs #9, #10, #12 and #13 preserve exploratory phylogenetic/extraction/older collinearity development without keeping those scripts in the active root.

The active scientific pipeline is documented in `docs/PIPELINE_34SPECIES.md`. Files are removed from the active tree when they are reproducible from Git history and no longer required by that pipeline.
