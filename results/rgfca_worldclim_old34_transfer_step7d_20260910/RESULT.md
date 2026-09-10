# RGFCA Step 7D — historical WorldClim transfer

- discovery pure C/S complete: **76 C* / 10 S***
- reserve pure C/S complete: **80 C* / 5 S***

## Directional hypotheses

### H1a_mean_bio4
- discovery delta(S-C) **380.1107145501992**, raw p **0.99845**, Holm **1**, pass **False**
- reserve delta(S-C) **220.5267575073242**, raw p **0.955802**, Holm **1**, pass **False**
- two-tranche recurrent: **False**

### H1b_mean_bio15
- discovery delta(S-C) **-18.091435519667233**, raw p **0.0709965**, Holm **0.354982**, pass **False**
- reserve delta(S-C) **-4.100082559585573**, raw p **0.40313**, Holm **1**, pass **False**
- two-tranche recurrent: **False**

### H3a_spatial_niche_turnover
- discovery delta(S-C) **0.09010429645772322**, raw p **0.299435**, Holm **1**, pass **False**
- reserve delta(S-C) **0.08787555585224349**, raw p **0.387431**, Holm **1**, pass **False**
- two-tranche recurrent: **False**

### H3b_regional_centroid_sep
- discovery delta(S-C) **0.21004545359491666**, raw p **0.343683**, Holm **1**, pass **False**
- reserve delta(S-C) **-0.3679761030343257**, raw p **0.966752**, Holm **1**, pass **False**
- two-tranche recurrent: **False**

### H3c_regional_gaussian_overlap
- discovery delta(S-C) **0.037906063228258735**, raw p **0.640668**, Holm **1**, pass **False**
- reserve delta(S-C) **0.21331585608255715**, raw p **0.974101**, Holm **1**, pass **False**
- two-tranche recurrent: **False**

## H0 total niche-size negative control

- discovery delta(S-C) **9.37182419530179**, p **0.000249988**
- reserve delta(S-C) **9.661780091846277**, p **0.00109995**

BIO4/BIO15 are climatological seasonality, not dynamic year-specific climate. H8 remains unopened.
