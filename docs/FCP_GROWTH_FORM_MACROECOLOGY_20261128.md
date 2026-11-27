# Growth form and life cycle: macroecological results

## Completed 27 November analysis

Run 45980912758 (job 140444264618) passed; prior landscape ZIP and GIFT 3.2 join hashes authenticated. Results stored in results/fcp4475_lifehistory_landscape_verified_20261127.json. Source 4475 taxa; GIFT records available for lifecycle 3832, woodiness 4148, self-fertilization 1024.

Two pooled trait coefficients exclude zero under the approximate 18-comparison Bonferroni percentile intervals:
- Woodiness vs nonwoody, range-wide four-category photo mismatch (3967 species: 1005 woody, 2962 nonwoody): -0.041854; family interval [-0.069860,-0.011690].
- Woodiness x standardized thermal partition at 250 km (1124 species: 306 woody, 818 nonwoody): -0.011537; family interval [-0.021668,-0.001395]. This is a difference in association, not a direct test of a positive slope in each group.

Annual vs perennial diversity +0.023662 and white fraction -0.038517 have nominal intervals excluding zero, but family intervals include zero. Woody vs nonwoody geographic sorting main effects at both scales also lose support after multiplicity correction. Self-fertilization diversity and mean white fraction unresolved; spatial comparisons only 97 taxa and bootstrap unstable, so not evidence of no relationship. No robust within-genus confirmation: woody diversity bootstrap had 3480/3999 valid draws; spatial woody contrasts had only 27 mixed genera, below original criterion. Distinguish macroecological association from causal life history, population genetic or phylogenetic claims.

## Focused 28 November extension

Woodiness and annuality are confounded biological categories. Compare perennial woody vs perennial nonwoody and annual nonwoody vs perennial nonwoody. Existing frozen 4475-taxon GIFT/landscape ledger is the only input. Same species within each pair of categorical and continuous nine-palette outcomes, fixed 250 km spatial scale (matched 50/250 km parent), no observer splitting.

Primary pooled family: 2 comparisons x (2 diversity contrasts + 2 spatial main contrasts + 2 thermal interaction contrasts) = 12 coefficients, approximate Bonferroni 95% intervals from 9999 fixed-design whole-genus Rademacher score draws (HC2 leverage correction). Mixed-strategy genera allow a separately reported genus fixed-effect contrast using nuisance projection; this avoids the rank instability of resampling rare covariate levels and does not constitute full phylogenetic adjustment. The estimation procedure differs from the previous random-genus percentile bootstrap; do not interpret differences as caused solely by holding life cycle fixed. Unobserved and variable traits remain explicit in the all-taxon ledger. Within-genus sparse groups may remain unestimable.

Question: does the lower within-species photo-colour diversity and weaker thermal-landscape association of woody plants persist after holding perenniality constant, or are they associated with life cycle? Generation time, actual longevity, breeding system, pigment physiology and selection are not measured. Pending genuine-source workflow; original manuscript and prior frozen results remain unchanged.
