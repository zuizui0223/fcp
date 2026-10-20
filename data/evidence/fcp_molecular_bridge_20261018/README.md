# FCP molecular evidence bridge — frozen 2026-10-18

This is a bounded archive-metadata screen plus manually verified primary-paper evidence. It is not an exhaustive census of flower-color mechanisms, a raw RNA-seq analysis, or a test of universal molecular adaptation.

## Inputs and reproduction

`curated_contrasts.json` contains 18 contrast records from 11 taxa. The output grain is **contrast × tissue**; sources or repeated sequencing runs do not create independent species. `independence_group` records shared experiments, and `independent_allele_origin_status` distinguishes independent orchid alleles from the shared poppy allele. `source_verification.csv` documents primary DOI, PMCID, sections and verification depth. Missing replication counts are explicitly unextracted, rather than zero.

`ena_title_screen_runs.tsv` is the unedited ENA Portal API response to the query recorded below and in `retrieval_manifest.json`; `ena_*_runs.tsv` are accession-specific follow-ups. `epmc_*.json` hold official Europe PMC DOI metadata. The manifest contains URLs, retrieval times, sizes and SHA256 hashes. `screen_decisions.json` assigns a disposition to **every** returned project. Only a manually reviewed primary paper can upgrade a metadata hit to a mechanism record.

Run the offline bridge from the repository root:

```bash
python scripts/analysis/build_fcp_molecular_bridge_20261018.py
```

A fresh live metadata snapshot can be obtained into a new directory with:

```bash
python scripts/analysis/fetch_fcp_molecular_bridge_metadata_20261018.py --outdir /tmp/new-fcp-ena-snapshot
```

A refresh is not expected to reproduce archived bytes: public database metadata can change. The frozen-input bridge is deterministic. No sequence reads are retrieved by either script.

## Search universe and omissions

ENA `read_run`, limit=0, fields and complete URL in the manifest:

```
tax_tree(33090) AND library_source="TRANSCRIPTOMIC" AND
((study_title="*flower*" AND
(study_title="*white*" OR study_title="*polymorph*")) OR
(study_title="*petal*" AND study_title="*white*"))
```

The query searches submitted **study titles**, not the complete literature, sample attributes, genome projects, functional genetics or every transcriptomic study for all 4,475 species. It returned 131 runs, 15 projects and 13 submitted taxon labels. Genus/family names, infraspecific names and hybrids are not guessed into roster taxa. Lack of a hit has no implication that a species lacks data or a mechanism. Titles lacking these English terms and mislabelled studies can be missed; the original Silene and Parrya seeds are obtained by explicit accession follow-up, outside this denominator.

## Evidence schema and claim limits

- T0: trait/ecophysiology; no pigment-gene cause identified.
- T1: chemical pigment or organ phenotype; no causal gene.
- T2: expression/metabolite association; no causal variant established.
- T3: variant association or mapped candidate. The `T3_mapped_candidate_variant` Orychophragmus record does **not** identify a validated causal nucleotide.
- T4: gene intervention evidence. The Gymnadenia RNAi was limited; it is not allele-specific rescue.

These are assay descriptions, not interchangeable confidence scores. `mechanism_class`, `causal_variant_identified`, `functional_intervention`, `loss_scope` and source notes must be read together. PAL/WAL is assigned to a particular measured contrast, never imputed to every white individual or to a species/genus. Anthocyanins and carotenoids remain distinct; Vinca chemistry is not inferred from RGB. “Whole-plant loss reported” is distinguished from chemical quantification of every organ.

The roster match uses only trimming, whitespace collapse and casefolding, followed by exact full-name matching. Taxon ID joins then connect the existing 4,475-species roster to the existing 1,574-species environment panel. No synonym or authority-stripping routine is applied. A matched molecular species provides a species-level annotation; **no FCP photo individual or observed white morph is genotyped by this bridge**. The strict matched-mechanism/FCP-individual count is therefore zero.

## Accession audit requiring attention

As retrieved, SRP033277/PRJNA229760 has an unrelated gull-fecal bacterial-community title and includes both Silene and gut-metagenome runs. The raw response is preserved unchanged. Only scientific_name=`Silene littorea` and strategy=`RNA-Seq` are selected: 18 runs, nine experiments and nine BioSamples; each sample has two technical runs. The other 14 runs are excluded. Nine color×stage samples are not biological replication within each color×stage cell.

PRJNA415931 reports six experiments but only two BioSamples. Its primary Primula paper reports three biological libraries per color. This is retained as a source-versus-metadata ambiguity; six sequencing libraries are not automatically reclassified as either six independent plants or two biological replicates.

## Outputs

`results/fcp_molecular_bridge_20261018/` includes the 31-row evidence table, all 4,475 coverage rows, exact matched/unmatched curated taxa, complete ENA project decisions, deduplicated run and library audits, a 42-row species×environment-variable bridge, summary counts and an execution receipt. Environment differences are inherited descriptive white-minus-nonwhite means. No regression, causal mediation, resource-cost estimate, universality test or mechanism-frequency estimate is produced.
