# FCP reproductive assurance vs photographic white-colour variation
## Source-identified cross-species feasibility contract, 2026-10-08

### Primary question
Can the hypothesized trade-off between coloured flowers attracting pollinators and white flowers gaining reproductive assurance from autonomous selfing be assessed in the 1,499-species FCP photographic frame, using DIRECTLY MEASURED reproductive traits from independent published data?

FCP currently confirms geographically structured photographed flower-colour variation, not heritable local flower-colour polymorphism or morph-specific fitness. A species' mean reproductive system cannot by itself establish a cost/benefit equilibrium between white and coloured forms within that species. This analysis tests cross-dataset source coverage before any regression.

### Immutable FCP source
Exactly three species-disjoint original measured photographic cohorts, 500/500/499 nominal species. Photo SHA256: discovery ee854126eed2cfe23e333abe2c28d14df24895a5c52cbc389c060a5a4d6f91f4; validation 0e2ed349122739eecfc725fb2d5e313d284cf30752da91f9a0429cff0eeaa5e6; third 57630fc9f281bce94a0c40a70aaf7bce879dde93d6154175adcd021e8f5c1186. Require >=40 classifiable photo records/species with valid geocoordinates. A photographed white+chromatic species must have >=5 on each side; not a wild genotype label.

### Three independent organism reproductive studies in the public island repository
External data repository: zuizui0223/island, pinned to SHA 92f007797fbbcb3c6743a4ffb2628fa609329499.

- Rodger et al. 2021, Science Advances, DOI 10.1126/sciadv.abd3524; original 1528-row pollinator-exclusion study (Figshare 10.6084/m9.figshare.14607882.v1). Only original measured auto-prefixed fruit and seed outputs represent autonomous fruit/seed production; observed numeric zero is a measured zero. Snapshot gzip SHA256 fa745c578f3537933fafedc1d36b4ea266348cd83d7f6cbb231c253b0f348d3f.
- Goodwillie, Kalisz and Eckert 2005, Annual Review, DOI 10.1146/annurev.ecolsys.36.091704.175539. Only genetic multilocus outcrossing fraction Mean-tm from original natural-population observations PopType=0, not experimental/crop/source summary proxy. Outcrossing rate is not direct bagged fruit production.
- Razanajatovo et al. 2016, Nature Communications, DOI 10.1038/ncomms13313. Accept directly measured Autofertility_index_FS or Autofertility_index_SFL (0..1). Separate Self-compatibility_index_FS/SFL does NOT substitute for autofertility. Keep index 0 as valid measured outcome.

These original source tables contain source-level taxonomic and population limitations. No assumed inferred self-compatibility, flowering-syndrome proxy or genus-level guessed value can enter this primary feasibility ledger. Island's own source policy explicitly distinguishes reported measurements from inferred proxies.

### Taxon and photo opportunity matching
Match exact species binomial (first genus + specific epithet) from the source. No fuzzy synonym fill or automatic reported/proxy reconciliation. If Rodger source exposes multiple plausible taxonomic column names, report HOLD_SCHEMA until unambiguously resolved, not whichever yields the most favourable overlaps.

Each original source contributes its own species with directly measured trait, distinct from the total number of articles, populations or photos. For each direct trait source separately, summarize how many of the source's unique taxa overlap all three original FCP high-depth photo groups (each >=40 classifiable photos), and how many show >=5 photographed white plus >=5 nonwhite photo labels. Do not infer genome colour variant, true plant population or actual morph-fitness payoff from shared scientific species name.

### Fixed comparison gate
An exploratory cross-species trait–photo colour association is allowed only when a single independent source has at least 80 unique FCP species with directly measured trait, each of the three disjoint photo cohorts has at least 20 such species, and each cohort has at least 8 white+chromatic photo-eligible species among them. If not, report HOLD: INSUFFICIENT DIRECT REPRODUCTIVE TRAIT OVERLAP. A fail is an identifiability limit, not evidence that reproductive assurance and flower colour have no ecological relationship.

No cross-source pooling of incompatibly defined autonomous reproduction, genetic outcrossing and self-compatibility scores. No general selection, stable population morph, pigment cost or fertility benefit claim can follow from a source-species photo match.

### Deliverables
Separate source-backed JSON and per-species opportunity CSV, synthetic tests enforcing the trait-family distinctions, and reproducible GitHub Actions that retrieve and verify the external repository commit and the original FCP image measurement SHA values. The current manuscript H1/H2/frozen geography remains unchanged.

### Source-confirmed Rodger name mapping correction (before first complete result)
The public island source implementation src/island_v2/rodger_2021_autofertility_checkpoint.py explicitly normalizes the original Rodger S3 genus.species field by replacing underscores with spaces. This primary-source mapping is now fixed by column name in our code; generic taxon fields are not tried and no result-favourable column search is allowed. The checkpoint also treats any nonmissing negative or invalid pollinator-exclusion measurement as invalid for that source row. The feasibility script has been aligned to this measured-evidence rule; observed numerical zero remains valid evidence of absent autonomous production.
