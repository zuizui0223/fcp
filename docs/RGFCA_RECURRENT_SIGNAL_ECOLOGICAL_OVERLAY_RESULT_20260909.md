# RGFCA recurrent-signal ecological overlay — exploratory result

Date: 2026-09-09 JST

Protocol freeze: `106ffd9986a5a2316a0909df14fb26f0dee28296`.

## Result

The prespecified primary ecological response was the flower-only component `S_F`. Across 200 recurrence realizations, none of M1–M3 showed calibrated support for a within-species association with absolute latitude or island share. All flower-only 2.5–97.5% realization intervals crossed zero and all two-sided exploratory p-like calibrations were >= 0.587.

The matched flower-minus-background component `S_D` showed a recurrent negative latitude slope for M1 (median -0.024616 per 10 degrees; 2.5–97.5% -0.042506 to -0.003676), but its p-like calibration was 0.0746. Decomposition identifies the stronger and more clearly calibrated pattern in the matched background component `S_B`: M1 background latitude slope median +0.025454 (0.006535 to 0.044687; p-like 0.0249) and M2 +0.016099 (0.002759 to 0.028515; p-like 0.0448). Because `S_D = S_F - S_B`, the negative matched M1 latitude tendency is therefore primarily consistent with a latitudinal background-colour gradient rather than a flower-only gradient in this reserve.

No island-share association was calibrated for the flower-only component. The strongest island-share diagnostic was M2 background (median +0.020707; p-like 0.169), also not calibrated against the species-conditioned null.

## Primary flower-only rows

| Mode | Predictor | Median slope | 2.5% | 97.5% | Fraction >0 | p-like | Informative species median |
|---|---|---:|---:|---:|---:|---:|---:|
| M1 | abs_latitude | 0.000832 | -0.009798 | 0.013901 | 0.575 | 0.9104 | 326 |
| M1 | island_share | -0.002459 | -0.028221 | 0.022833 | 0.430 | 0.8657 | 117 |
| M2 | abs_latitude | 0.003482 | -0.010676 | 0.018918 | 0.645 | 0.7662 | 326 |
| M2 | island_share | 0.004000 | -0.025480 | 0.032654 | 0.605 | 0.8358 | 117 |
| M3 | abs_latitude | 0.000468 | -0.007356 | 0.008044 | 0.570 | 0.8905 | 326 |
| M3 | island_share | -0.005734 | -0.021905 | 0.007055 | 0.225 | 0.5871 | 117 |

## Verification

- eligible photographs: 24,885; eligible species: 500; median admitted species per realization: 326; median states: 1,877;
- exact upstream mode-reproduction check maximum absolute numerical difference: `4.440892098500626e-16`;
- realization-level slope identity maximum error: observed `1.0408340855860843e-16`, null `1.9905951886833861e-16`;
- a complete independent rerun produced byte-identical CSV outputs.

## Claim boundary

This is exploratory, photo-derived, within-species evidence. It does not establish adaptation, pollinator attraction, causal island effects, or a universal latitude rule. The strongest ecological overlay in this stage is a background-colour latitude gradient, not a demonstrated flower-only gradient. Measurement-validity limitations retained elsewhere in the RGFCA branch remain unchanged.
