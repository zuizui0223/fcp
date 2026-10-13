# Sisyrinchium M-flower colour: distance-block-conditioned microspatial sorting

New ecological question: does M-purple/M-white quadrat clumping within observed source lines persist after conditioning on an additional source-distances-within-lines inventory? Original field data from Takahashi, Takakura & Kawata (2016), PLOS ONE DOI `10.1371/journal.pone.0164381`, original Dryad `10.5061/dryad.j8q34`, publicly mirrored by Zenodo record 5003124 (same DOI, CC0). Research is **retrospective and exploratory**, using exactly the same source already used in PR #183. No novel source individuals or causal mechanism.

## Recompute from the ZIP

1. Use Python 3.13 + `numpy==2.3.5`, `matplotlib==3.10.8`, `pytest==9.0.2` (`requirements-lock.txt`).
2. `python src/distance_block.py --source source/frequency.txt --source-zip source/Sisyrinchium_data.zip --parent source/parent_line_result.json --output results`
3. `python src/verify_independent.py --root .` (standalone exact `Fraction` implementation, no import of primary script; rebuilds all full 9999 conditioned simulations)
4. `pytest -q tests/test_distance_block.py`
5. `python src/make_figures.py --root .`

On GitHub, the source ZIP is deliberately **not committed to the repository**. The CI re-downloads the original from its Zenodo mirror using the source-preserving **parent recovery script**, verifies original ZIP MD5 and embedded `frequency.txt` SHA256, and copies the exact frozen PR #183 `result.json` (SHA locked) before replaying the analysis. The stand-alone downloadable conversation ZIP **does include the authenticated original source**, input ledger and all code/results. Never use invented local plant counts.

## New exploratory outcome

11 parent line-eligible siteyears; the pre-declared distance-bin power gate retains 10 siteyears (9 named sites), holding Hirose Y13 with only 27 M individuals in bins supporting within-bin shuffles. Original 512 quadrats, violet M 5081, white M 8723, other species L 4463 / S 29. Three quadrats with `distance=-1` are retained as `UNKNOWN_DISTANCE` and not used in the near/far ordering.

Distance *within recorded line* median-unique-distance split: mean observed `B=0.394421118760634`, source-matched expected `E0(B)=0.33060691238523`, **excess=0.0638142063754042**, plus-one p=`0.0001` after 9,999 hypergeometric draws. **All 10 siteyears positive**, latest-year-only 9 named sites mean excess `0.0626068174272312`, p=0.0001; leave-one-out range 0.0542448–0.0695993, all p=0.0001. Adjacent-distance-**pair** sensitivity: excess=0.063773430707215, p=0.0001, 10/10 positive; pair-bin grouping is *not always a refinement* of median near/far grouping, so the two estimates are **alternative distance-blockings**, not a causal gradient of stricter adjustments.

Restricting the previous PR #183 line-only clumping benchmark to the **same ten siteyears** gives excess +0.206335656841248; conditional near/far bins reduce this to +0.0638142063754042. **The 69.1% difference of two conditional null references is not an environmental explained fraction or water-distance-mediated causal effect.** Residual within-bin distances, humidity, other species, genotype/population history, and spatial autocorrelation remain unrestricted. P values benchmark label exchangeability within source-defined bins, not a randomized experiment or phylogenetically independent population sampling. No inference that L-species interference causes this pattern. These are not genetic F_ST.

## Scope and evidence

- `docs/PROTOCOL.md` is frozen before new scores (GitHub commit `329099093265e46bca0c82b51804cd061a0a29a9`). The previous result was already exposed and used to design this stronger null, so this is not prospectively registered in the scientific preregistration sense.
- `results/result.json` all 16 siteyear holds; `results/quadrat_distance_bin_ledger.csv` all 512 original quadrats; `siteyear_distance_bin_ledger.csv` and `quadrat_distance_bin_ledger.csv` are audit ledgers.
- `src/verify_independent.py` separately reconstructs exact analytical nulls and 9,999 random placements per block, with hand-enumerated simple controls.
- Main FCP H1/H2 and parent PR170–183 field/photo outcomes unchanged. PR #182 Boechera result is fragile to removal of Judd Falls; do not promote this Sisyrinchium result to a cross-species generality claim.
