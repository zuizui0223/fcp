# White/pigmented proximal-mechanism literature audit — 2026-09-25

Status: **exploratory mechanism audit specified before the targeted OpenAlex results are opened.**

This line is separate from the frozen New Phytologist manuscript. The white/non-white phenotype and the BIO5 result are already known, so this is not a prospective test of phenotype association. Its purpose is narrower: determine whether direct molecular studies of the same FCP species repeatedly implicate a common pigment-control architecture.

## Target species

Build the target list from the frozen third-cohort measurement artifact and the frozen response-blind highlight technical seal.

Retain species with:
- globally classifiable rows in the four frozen states;
- valid near-clip information;
- exact frozen high-clip exclusion;
- >=5 measured white rows;
- >=5 measured non-white rows.

This should reproduce the 281-species primary white-environment panel. If it does not, stop.

## Search source

OpenAlex Works API, queried by exact canonical species name.

For every target species, use two fixed searches:

1. pigment mechanism:
   `"<species>" (anthocyanin OR flavonoid OR pigment OR pigmentation OR MYB OR bHLH OR WD40 OR DFR OR CHS OR CHI OR F3H OR F3'H OR F3'5'H OR ANS OR UFGT)`

2. white/pigmented contrast:
   `"<species>" ("white flower" OR "white-flowered" OR "white morph" OR acyanic OR albino OR purple OR pink OR red OR blue) (gene OR transcriptome OR metabolome OR expression OR anthocyanin OR flavonoid OR pigment)`

No query is retuned after results are inspected.

## Candidate eligibility

A work is a candidate only if title or abstract contains:
- the exact target species binomial;
- a floral-display term (flower, floral, petal, corolla, tepal, labellum, bract);
- a pigment/molecular term.

Horticultural/cultivar-only, transgenic-only and artificial-mutagenesis studies are flagged and cannot support the natural-polymorphism synthesis unless the abstract also explicitly states natural/wild population or naturally occurring morph evidence.

## Mechanism coding

Automated coding is navigation-only. Each candidate receives fixed flags:

- `anthocyanin_pathway`: anthocyanin / flavonoid biosynthesis evidence;
- `MYB_MBW_regulation`: MYB, bHLH, WD40/WDR or MBW regulation;
- `structural_gene`: CHS, CHI, F3H, F3'H, F3'5'H, DFR, ANS, UFGT or related pathway enzyme;
- `expression_regulation`: differential expression, transcriptome, expression, transcriptional regulation;
- `metabolite_evidence`: HPLC, metabolomics, anthocyanin/flavonoid content;
- `loss_or_reduction_language`: loss, absent, reduced, suppression, downregulation, inactivation, acyanic;
- `white_contrast`: white/acyanic/albino contrasted with a pigmented state;
- `temperature_link`: heat, temperature, warm, thermal;
- `natural_context`: wild, natural population, polymorphism, morph, geographic variation.

## Evidence tiers

- **Tier A — direct natural molecular mechanism:** exact species, natural/wild morph contrast, white/pigmented contrast, and a gene/regulatory/pathway mechanism.
- **Tier B — direct natural biochemical mechanism:** exact species, natural/wild morph contrast, white/pigmented contrast, and pigment/metabolite evidence but no causal gene/regulator.
- **Tier C — mechanistic candidate:** molecular/pigment evidence but natural white/pigmented contrast not established in the abstract.
- **Excluded/artificial:** cultivar/breeding/transgenic/mutagenesis only without natural context.

Automated output does not promote a work to Tier A/B solely from keywords. Tier A/B require manual verification of the title/abstract evidence snippet in a second, explicit adjudication file.

## Comparative question

The primary descriptive synthesis is:

1. How many of the 281 target species have any direct natural white/pigmented molecular/biochemical evidence?
2. Across how many plant families?
3. Among manually verified Tier A/B species, what fraction implicate anthocyanin/flavonoid loss or downregulation?
4. Among Tier A species, how often is the implicated change regulatory (MYB/MBW/expression) versus a structural enzyme?
5. Is direct temperature-sensitive pigment regulation documented in any of the target species?

No p-value is required for the literature evidence proportions. This is an evidence map, not a random sample of all angiosperms.

## Interpretation boundary

A repeated anthocyanin/MYB pattern would support a **shared proximal developmental/genetic route** to white/pale states. It would not establish:
- that white is evolutionarily derived in every species;
- that BIO5 causes those molecular changes;
- that the same mutation recurs;
- that pollinators are irrelevant;
- prevalence among all flowering plants.

The mechanistic model to evaluate is:

`temperature / other ecological filters -> pigment regulatory state -> anthocyanin production -> white versus pigmented phenotype`

with the ecological input and molecular route treated as separate inferential layers.
