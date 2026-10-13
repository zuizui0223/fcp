# Within-line distance-conditioned colour clumping: prospective-analysis specification relative to this new endpoint

Date 2026-10-13 (JST). **Source-exposed, retrospectively designed exploratory analysis**, not an external pre-registration. Already observed the PR182 and PR183 within-line source outcomes and site distributions before defining this test. Freeze this file in the GitHub branch **before computing the new distance-conditioned statistics or p-values**. Do not upgrade its p-values to confirmatory selection evidence.

## Fixed ecological question

Within the original Takahashi et al. (2016) *Sisyrinchium* M-species field quadrats, does within-site, within-recorded-line **purple/white clumping among quadrats** persist after preserving not only every line's M-species purple/white inventory but also the inventory of quadrats that are **near versus far from the water margin**? Distinguish a coarse water-margin gradient from additional plot-scale colour heterogeneity. Positive answer would show unresolved fine-scale structure **conditional on the chosen two-bin coarse distance model**, not environmental adaptation, reproductive interference, genetics, or a demonstration that distance never matters.

## Source, observations and scope

Original field census archived at Dryad DOI 10.5061/dryad.j8q34; valid public same-DOI mirror Zenodo record 5003124. Source ZIP `data.zip` length 127502, MD5 `008ec196aa5acb769dbfd36e835e32ba`, expected internal `data/frequency.txt` SHA256 `a12cdf46483a6d8b78bd00110c196a25f3c8849967ae8a2c37d6aee3d29f516e`. Original data 512 quadrats, 16 (year,population) strata, 13 named sites, 2011/2013. Rows are **quadrat** counts, NOT independent evolutionary populations. Fields `v,w` = M-purple/M-white; `l,s` = *other* separate white-flowered species and must never count toward M-white. Quadrats with n=v+w=0 are `NO_TARGET_PLANTS`, still included in original source ledger but not in within-M colour denominators. Other species counts are not used for this fixed clumping endpoint. Three quadrats have recorded negative `distance=-1`, which should be classified as **UNKNOWN_DISTANCE**, never treated as 1 metre closer to shore; those bins are frozen and unshuffled between distance categories. Do not impute precise coordinates, meter distances for source labels, or genets.

## Fixed eligibility (not tuned to result)

Start with the **same 11 siteyears** eligible for PR183 secondary within-line clumping (source `results/result.json['full']['results']['clumping']['names']`), with ≥20 M individuals and ≥5 of each colour at siteyear level, and ≥2 lines with ≥2 occupied quadrats and both M morphs. Keep a source ledger of all 16 original siteyears, including held strata. For new inference require at least **two independent recorded lines**, each containing >=1 distance bin with >=2 occupied quadrats and both colour labels in that bin; require >=30 total M individuals among informative bins. A line with a single plot or bins containing a single morph have no within-bin randomization opportunity. Failing siteyear is `HOLD_DISTANCE_BIN_SUPPORT`, not null evidence. Eligibility for the new endpoint is evaluated before its magnitude or p-values are inspected.

## Primary distance strata and exact conditional null

Within each original **(year,population,line)**, use the M-present quadrats with nonnegative recorded water-edge distance. Let `u` be the sorted **unique distance values** in the line, not sorted observed morph rates. Split distance into two classes by the *median unique-distance rank*: `NEAR` if `distance <= u[(len(u)-1)//2]`, otherwise `FAR`. This keeps distance ties intact, so a common distance value is never arbitrarily divided based on plot row order or flower colours. If the line has only 1 unique nonnegative distance value all its plots are NEAR. If source `distance<0` classify as `UNKNOWN_DISTANCE` even if zero other plots have that value. Never bin by source morph counts or tune the cutpoint on clumping strength.

Let each siteyear have plot counts `n_j=v_j+w_j`, `v_j`, total `N,V`, `p=V/N`. Observed phenotypic Gini-Simpson between-plot fraction

`B = sum_j n_j*(v_j/n_j-p)**2/(N*p*(1-p))`, for p strictly between 0 and 1, excluding n=0 plots.

Under the null **conditionally randomize individual M purple labels only among quadrats within the same (year,population,line,NEAR/FAR/UNKNOWN_DISTANCE)** while retaining *exact* `N_h,V_h` per line-bin and `n_j` of each plot, and without changing independent other-white-species counts `l,s`. The analytic null expectation is

`E(B | line-bin inventories) = B_between_bins + sum_h N_h*p_h*(1-p_h)*(m_h-1)/(N_h-1)/(N*p*(1-p))`,

where `B_between_bins = sum_h N_h*(p_h-p)**2/(N*p*(1-p))`, `p_h=V_h/N_h`, and m_h is the number of M-occupied quadrats in the bin; 0 contribution for N_h<=1 or m_h<=1. Positive `Delta_distance = B-E(B | line-bin inventories)` represents **sub-bin patchiness beyond this two-class distance conditioning**. This is a source-scale phenotype allocation, *not genetic F_ST*, climate explained fraction, or adaptation.

Primary study statistic: arithmetic mean of `Delta_distance` across eligible **siteyears**, each exactly one vote. Use `numpy.random.Generator.multivariate_hypergeometric` for 9,999 draws; seed **2026101312**, calculate one-sided plus-one `p=(1+#(null >= observed))/10000`. A positive p-bound only describes conditional label exchangeability *at fixed local water-distance bins*; real plants, genotypes and habitats were not randomized.

## Two fixed sensitivity checks and interpretive stops

1. **Adjacent-distance-pair null (finer conditioning):** for each line, sort unique nonnegative distance values and create bins of consecutive *two unique distance ranks* `[u[0],u[1]], [u[2],u[3]],...` with ties in one bin; unknown distance stays isolated. Use the same analytic formula and `R=9999`, seed **2026101313**, with its own support gate ≥2 informative lines and ≥30 M individuals in multi-quadrat mixed bins. A null/no-power outcome is not evidence against prior PR183 clumping; a positive outcome would show heterogeneity down to this more restricted distance scale, which still varies among sites.
2. **Latest year per named population:** use only the later year when the same `population` occurs in both 2011/2013 among the primary eligible siteyears. Require ≥3 independent named populations, otherwise HOLD. Reuse the **same** simulated draws restricted to selected siteyear columns, report number and names. This prevents double-voting a repeated study site and is not an independent sample.

Posthoc audit after frozen main and sensitivity results: leave one siteyear out to test fragility; sign count of siteyear Delta, occupied plot and eligible-bin coverage. Never change main endpoint or p-value after seeing results. No multi-species generality claims from this one *Sisyrinchium* study; earlier *Boechera* clumping was dominated by Judd Falls.

## Required reproducibility

Source bytes and Zenodo MD5+size, internal frequency SHA256, 512 key checks, two bin ledgers (all recorded plot keys), observed B, B_between_bin, analytic null E, Delta for every eligible/HOLD siteyear; exact integer simulation inventory conservation; independent Fraction calculation without importing main analysis; tiny hand-enumerated cases and deterministic seeds; original PR183 unchanged. Save execution package/report and remote checked result. Full repository tests and newly collected biological data are out of scope. Existing FCP manuscript H1/H2/PR170-183 unchanged.