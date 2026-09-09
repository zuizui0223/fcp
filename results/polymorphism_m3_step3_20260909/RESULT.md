# Polymorphism M3 Step 3 — result

**Frozen diagnosis: `FLOWER_BIOLOGICAL_SUPPORTED`.**

## Discovery (369 species; M3 learned independently in reserve)

- rho(D, mean M3 full-12): **0.255991**, permutation p **4.99975e-05**
- rho(D, SD M3 full-12): **-0.015573**
- rho(D, mean M3 after removing green/brown/black): **0.300195**, p **4.99975e-05**
- partial rho(D, mean M3 | flower nuisance fraction): **0.271723**, p **4.99975e-05**
- rho(D, flower nuisance fraction): **0.065295**

## Reserve component diagnostic (363 species; diagnostic, not independent)

- rho(D, mean flower-only M3): **0.316098**
- rho(D, mean background-only M3): **0.007629**
- rho(D, mean flower-minus-background M3): **0.322102**
- partial rho(D, flower M3 | background M3): **0.316098**
- rho(D, flower nuisance fraction): **0.156996**

## Interpretation boundary

- diagnosis: **FLOWER_BIOLOGICAL_SUPPORTED**
- direct flowering-stage inference: **not allowed** (no frozen phenological-stage label)
- Step 1 Claim 2 remains rejected regardless of this diagnosis.
