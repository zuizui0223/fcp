## Same-species bridge evidence

Three systems sharpen the heat–pathway synthesis.

### Parrya nudicaulis — strongest natural population bridge

The natural purple–white polymorphism is geographically structured across Alaska. White-flowered frequency increases with growing-season temperature. The same system shows that white petals lack anthocyanins because flux is blocked near the entry to the anthocyanin biosynthetic pathway, with petal-specific regulatory downregulation near CHS rather than a shared coding knockout.

This is the closest existing natural population analogue to the FCP working hypothesis:

`
warmer environment -> lower pigment-pathway output -> white endpoint
`

It still does not prove temperature is the causal selective agent because the population-level temperature association is observational.

### Moricandia arvensis — direct experimental heat × pigment-network bridge

Moricandia provides the strongest experimental bridge because the same system contains both a temperature-dose response in floral anthocyanin and a public floral transcriptome.

#### Temperature dose in published Source Data

The published `Figure 3F-J` Source Data contain a subgroup of 15 individuals followed through:

`
14.17 C spring -> 23.75 C mild summer -> 28.75 C hotter summer
`

The mild- and hotter-summer regimes use the same stated summer photoperiod (16/8 h day/night).

After averaging the two flower measurements within individual-period:

- mean cyanidin at 23.75 C = **1.893**;
- mean cyanidin at 28.75 C = **0.531**;
- hot/mild ratio = **0.281**, or a **71.9% reduction** in mean anthocyanin;
- **13/15** individuals decreased;
- paired Wilcoxon one-sided `p = 0.000427`;
- directional sign-test `p = 0.00369`.

The result is reversible in a separate 14-individual sequence:

`
23.75 C mild summer -> 14.17 C spring return
`

where mean cyanidin rises from **1.236** to **4.708** (**3.81x**; one-sided Wilcoxon `p = 0.000305`).

A third 14-individual sequence going from spring-like to hot conditions shows cyanidin falling from **5.690** to **0.179**, with **14/14** individuals decreasing (`p = 6.10e-5`).

The period-2-to-period-3 change differs strongly between mild->hot and mild->spring-return sequences (Mann-Whitney `p = 3.72e-5`), arguing against a generic third-period decline.

This is post-publication/post-preflight evidence, and the sequential mild->hot subgroup can retain period/order carryover. It therefore supports a within-species mechanistic bridge rather than a new independent confirmation.

#### Matching pigment-network response

Reanalysis of all public supplementary workbooks and Source Data shows that every normalized FCP white-state node class represented in the Moricandia significant-DEG table is downregulated in summer-white flowers:

- **MYB90**: logFC -4.056, FDR 4.81e-11;
- **CHS**: logFC -1.293, FDR 0.034;
- **U75C1/U78D2**: four significant rows, median logFC -1.5315;
- **TTG1 / WD40**: logFC -1.411, FDR 0.0279.

Additional anthocyanin-pathway context genes **PAL, 4CL and DFR** are also significantly lower.

Thus Moricandia supports the chain:

`
greater thermal severity
        -> lower floral anthocyanin
        -> lower output / regulation at FCP-overlapping pigment-network nodes
        -> white / achromatic floral phenotype
`

The RNA-seq spring-versus-mild-summer comparison changes both temperature and photoperiod, so the transcriptomic step cannot be assigned to temperature alone. The separate mild-versus-hot cyanidin contrast is valuable precisely because those two summer regimes share the stated 16/8 h photoperiod.

### Ipomoea purpurea — critical counterexample

The same anthocyanin network is environmentally temperature-sensitive in natural flowering-season data, but the sign is not a simple universal whitening response: principal floral anthocyanin was positively associated with preceding temperature, with transcript regulation likely involving myb1.

This is useful rather than inconvenient. It shows why a global BIO5 main effect can fail even when the underlying pigment network is genuinely thermosensitive.

## Revised mechanistic model

The strongest formulation is:

`
(genetic/regulatory state) × (thermal environment)
        -> anthocyanin/flavonoid pathway flux
        -> white/pale versus pigmented floral phenotype
`

The expected cross-species signal is **not** a universal positive BIO5 coefficient. The expected biological generality is that temperature and genotype can perturb overlapping nodes of the same pigment network, with sign and magnitude depending on genotype, developmental stage, light and metabolic context.

This formulation reconciles three otherwise awkward results:

1. the global FCP phenotype-space axis recurrently aligns white versus nonwhite;
2. direct natural molecular systems converge on multiple parts of the pigment network;
3. BIO5 does not replicate as a universal cross-species predictor.

The next decisive experiment is therefore not another global temperature regression. It is a **natural colour genotype/morph × temperature factorial experiment** measuring colour, anthocyanin chemistry and the same pathway nodes in the same flowers.
