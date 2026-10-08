# FCP NEXT-2000 — independent species confirmation contract
## Frozen BEFORE acquisition of any new flower-colour pixels | 2026-10-08

### 1. The prior 1,499-species work is discovery, not the confirmation set

Treat ALL outcome-exposed FCP work through 2026-10-08 as hypothesis formation: the original discovery+validation 1,000 high-depth species; the unsuccessful one-shot P500 500 selected species (499 measured but all 500 identities permanently excluded); and the later third high-depth species selection of 500 (499 biological measurements, all 500 permanently excluded). All exploratory PR #127–#132 and any environmental/pollinator/white-vs-coloured mechanism suggestions are TRAINING, not independent confirmation.

**This next study has one focal confirmation target:** within-species photographed flower-colour states are more geographically homogeneous among <=50-km conspecific photo pairs than expected from each newly sampled species' exact sampled visible four-colour composition **and calendar month**. This is a comparative ecological regularity about *geographic organization* of sampled phenotypes; it is NOT genetic balancing selection, a worldwide white allele fitness trade-off, or species-polymorphism prevalence.

### 2. Provenance and untouched candidate numbers

A metadata-only iNaturalist opportunity census, made without flower-colour outcomes, identified 42,111 species in 162 equal-area grid cells. Its observer-capped >=100-photo capacity layer was U100=4,730 species, not all angiosperms. See the pre-existing documented frame at docs/POLYMORPHISM_42111_FRAME_PROVENANCE_20260918.md and original Step-8B canonical census SHA256 1141bfa1ea5f1f2ffa995e6695483a0fc1ec86add7a00e7c76d6da447e02b22c.

Exclude from those U100 opportunities:
- former discovery and reserve: 1,000 high-depth species (already removed in the historically frozen P100=3,730);
- P500 previously prospectively *selected*: 500, even its single failed measurement;
- previous third-cohort *selected*: 500, even its single failed fresh photo opportunity.

Thus U100 - all historical allocations = **2,730 species**, not 4,730, not 3,230, and not a new reinterpretation of the old 1,499 photo species. The historical pre-third candidate-pool size must be exactly 3,230, SHA256 7fc0337074b55a91c0b50773a8d5c3ef82e877074c8b3cacc21f9af58e0ba77e.

Original parent frozen metadata selection files, with strictly permitted fields [inat_taxon_id, species, after_observer_cap]:
- p100_outcome_blind_pool.csv from Actions workflow run 34708044962, artifact 10302477571; SHA256 1473aad680fe2fa84903c5e11ee104fd0eceef828957f16c2ef1a35f9dd6993c;
- p500_frozen_selection.csv from same artifact; SHA256 f54d07fb2338a20a3f92808b58f0ebc948ab0d5af3cb01fb26702f861184b7a4;
- third selected species manifest from immutable commit 7e538e5c51c05a7cc47b2fcf53eea92634c8a863, SHA256 16db6a2fc265f0ab20e7dae4e6f9058b907f5c19af3ddab5b6077797dd1def59.

If any source SHA/count/identity differs, terminate without any image opening or taxon replacement.

### 3. Single outcome-blind split, two separate roles

From the 2,730 truly untouched species, use ONLY its taxon ID and canonical original scientific name to compute SHA256(SALT | taxon_id | species), sort ascending with taxon ID tie break.

Fixed salt: FCP_NEXT2000_SPATIAL_20261008_FIXED_V1.

Freeze:
- the FIRST **2,000** identity-hashed species for the new confirmatory cohort;
- the remaining **730** in an unexamined reserve for a separate, future-only independent replication (NEVER choose a replacement after seeing any pigment result or failed species);
- canonical 3,230 parent and 2,730 residual opportunity ledgers with source and output SHA256;
- exact 2,000/730 selected manifest byte hashes and taxon overlap = zero with every previously selected organism.

The 2,000 new species could yield up to 200,000 new photographic records at 100 selected photos/species, and potential total photographed species 1,499 + 2,000 = **3,499**. This is an upper design target, NOT an assertion that 2,000 new taxa will yield >=40 classifiable photos or that the new measurement model will pass. The 2,730 are historical metadata opportunities; at 2026-10-08 the fresh iNaturalist API may return fewer.

Taxon sampling is unbiased *with respect to previously measured floral colour* inside the qualified U100 opportunity frame. It is not a probability sample of all global flowering plants, because geographic citizen-science observation and photo availability are nonrandom.

### 4. Observation-metadata first; no hidden image access

After committing the actual selected manifest and its SHA to Git, a *distinct prospective stage* may attempt fresh iNaturalist metadata for the 2,000 species in fixed, bounded batches, using the already tested source-photo constraints: research-grade, flowering annotation, georeferenced unobscured location, <=5,000-m positional accuracy, supported licenses, observer cap 2, geographic maximin opportunity sampling, **exactly 100 photos/species** where possible. Source and photo IDs must not overlap any previous discovery/validation/P500/third measured or candidate-ID pools. Use the historical 2026-09 metadata query as technical starting design but freeze final API query, request pagination, seed and retry cadence BEFORE collecting the first new candidate photos. A good historical sample opportunity does NOT authorize a silent species substitution, repeated favourable metadata draws, capacity-gate relaxation, or random re-query after opening biological outcomes.

The first fresh metadata-only attempt is a feasibility gate; image pixels and morphological colour values remain unopened. No biological claim can be made from this gate. Materialize exact observation/photo-ID ledgers with SHA256 and observer identity, species counts, missingness and regional spatial opportunity. Reserve 730 stay untouched; if <1,000 of the 2,000 can supply 100 new ID-distinct photo rows, report a terminal capacity failure and request a genuinely new preregistration for another study, NOT a rescue selection from reserve species.

### 5. Confirmatory hypothesis and frozen analysis (after metadata and technical gates)

**Single primary biological estimand:** mean across equal-weight new species of 50-km local four-colour pairwise discordance depletion relative to that species' own month-constrained composition-preserving 199-permutation null.

Use the previously outcome-exposed and stabilized *calendar month* implementation as the PRE-FROZEN ANALYSIS SOURCE for the NEW sample only:
- script: scripts/analysis/run_polymorphism_specieswide_space_season_20261008.py
- immutable historical implementation commit: 969c67d274868b5dc4e8c2c8c9fa5d1f1beb5527
- labels: white, yellow_orange, red_pink, blue_purple as the existing fixed palette;
- eligible biological species: >=40 globally classifiable georeferenced photo records, >=40 with valid observation dates, >=30 local conspecific <=50km photo pairs;
- label permutations: 199, fixed seed and taxon ID; all species' monthly colour counts held fixed;
- species bootstrap: 1,999 replicates; no photo-count weighting;
- matched-null no-swap strata: retain geographically evaluable species with ZERO *identifiable* excess, document n_identifiable separately; do NOT treat this as biological zero.

**Pre-fixed confirmation gate** in the new, disjoint 2,000-species sample:
1. >=300 geographically evaluable species and >=250 month-null identifiable species (among new photo specimens only);
2. species-equal mean observed-minus-month-null spatial-depletion excess >0;
3. plus-one upper-tail matched null P <= 0.05 and 95% species bootstrap lower endpoint >0;
4. independently photographed (two nonempty, different observer IDs) local pair mode has >0 mean effect, P <= 0.05 and a strictly positive bootstrap lower endpoint. For the observer-disjoint set, >=240 null-identifiable species required.

No primary statistic may be redefined after photo opening, or swapped to year×month, alternative band/radius, individual white coefficients or a favourable lineage if the primary fails. Confirmatory classification must be PASS, FAIL_ADMISSIBLE, INSUFFICIENT_COVERAGE or TECHNICAL_FAILURE with each of the exact gates reported.

**Secondary diagnostics not allowed to rescue primary:** year×month null, 25/100/250km thresholds, nonwhite-only analysis, alternate continuous nine-colour JSD, genus-balanced and leave-one-genus-out checks, regional 0–30/30–60/60+ latitude replication, observer/multiple-photo duplicate audits. These are falsifications or generality diagnostics only, with full opportunity counts and uncertainty. Do not report a negative latitude-band association as proof the true effect is exactly zero.

**White colour biological validity restriction:** high-clipping light/exposure sensitivity discovered in previous H2 work means the coarse white label cannot be interpreted as absence of anthocyanin. Independent blinded crop/organ/exposure verification must be validated before any novel claim about evolutionary gain/loss of the white morph. The primary all-four photographic colour topology remains a source-conditioned phenotype-space result.

### 6. Technical execution firewall and scientific status

Before any new photo pixels are opened:
- run the full source-photo measurement, reassembly, JSON serialization, CI gate/figure and upload workflow ON SYNTHETIC DUMMY DATA;
- freeze exact runtime libraries, trained segmentation model, palette/lighting quality gate, source execution SHA, a single bounded biological acquisition/measurement authorization and artifact paths;
- prove the outcome receiver can durably serialize the final result and upload it before biological execution, to avoid repeating P500's late serialization failure;
- validate one-run metadata/photo ID exclusion and disjoint species/photographs across all prior cohorts.

The new 2,000-species selection and metadata-only opportunity alone **does NOT** constitute independent biological confirmation, and does not authorize computationally expensive photo retrieval by itself. Freeze capacity and the real image processing budget before execution. No new confirmatory effect size, P value, white pigment interpretation, balancing selection, reproductive fitness or global genetic polymorphism prevalence is claimed until a full independent, validated first-and-only measurement run is complete.

### 7. Version and evidence separation

The existing New Phytologist manuscript and all old H1/H2, third cohort and PR #127–#132 results are **discovery** and immutable. The new prospective results must live under this next2000 project path with a distinct artifact, a prereg chronological ledger and an untouched selection receipt. The reserved 730 taxon identities and any new photos must NEVER be searched for alternative hypotheses until they receive a separate preregistration and authorization.
