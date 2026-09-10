# Polymorphism Step 6b — opportunity decomposition result

**Span-specific verdict: `SPAN_ADJUSTED_EFFECTS_POSITIVE_BUT_NOT_TWO_TRANCHE_HOLM_SUPPORTED`.**

## Span-only D–spatial adjustment

- discovery raw D: partial rho = **0.103241**, raw p = **0.0489476**, Holm p = **0.0978951**
- reserve raw D: partial rho = **0.101872**, raw p = **0.0500975**, Holm p = **0.0978951**
- two-tranche span-only gate: **False**

## D_unbiased + span-only sensitivity

- discovery: partial rho = **0.102523**, Holm p = **0.0976951**
- reserve: partial rho = **0.102140**, Holm p = **0.0976951**
- two-tranche D_unbiased span-only gate: **False**

## >=10% threshold + span-only

- discovery adjusted spatial-rho difference = **0.015875**, Holm p = **0.111494**
- reserve adjusted spatial-rho difference = **0.011662**, Holm p = **0.184641**
- two-tranche threshold gate: **False**

## n_classifiable-only sensitivity

- discovery raw D: partial rho = **0.156401**, p = **0.00214989**
- reserve raw D: partial rho = **0.077892**, p = **0.139543**
- discovery D_unbiased: partial rho = **0.156182**, p = **0.00244988**
- reserve D_unbiased: partial rho = **0.078191**, p = **0.135443**
- rho(D, n): discovery **-0.537139**, reserve **-0.555901**
- rho(D_unbiased, n): discovery **-0.539973**, reserve **-0.558194**

## Genus clustering decomposition

- discovery raw D: gain=0.203129, p=0.00289986 ; span-only: gain=0.211267, p=0.00149993 ; n-only: gain=0.179873, p=0.00319984 ; span+n: gain=0.161050, p=0.00679966
- discovery D_unbiased raw: gain=0.203940, p=0.00234988 ; span-only: gain=0.212332, p=0.00154992
- reserve raw D: gain=0.154462, p=0.0218989 ; span-only: gain=0.152109, p=0.0234988 ; n-only: gain=0.085055, p=0.114394 ; span+n: gain=0.084398, p=0.112944
- reserve D_unbiased raw: gain=0.154258, p=0.020999 ; span-only: gain=0.154105, p=0.020749

Raw photo opportunity is fixed at 100 photos/species in both tranches. `n_classifiable` is therefore post-classification yield rather than raw sampling effort. Span-only survival does not erase the separate sensitivity to classifiable yield.
