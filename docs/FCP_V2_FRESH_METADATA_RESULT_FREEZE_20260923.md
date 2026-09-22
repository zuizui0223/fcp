# FCP v2 fresh-metadata terminal result freeze — 2026-09-23

Status: **PASS / PRE-PIXEL TERMINAL DENOMINATOR FROZEN**

The one-shot metadata-only run completed successfully before any FCP v2 image pixel or biological colour outcome was opened.

## Frozen terminal denominator

- Panel P queue: 300 species
- Panel P species with exactly 100 fresh eligible rows: **299**
- Panel P terminal species: **200**
- Panel P terminal rows: **20,000**
- Panel N queue: 300 species
- Panel N species with exactly 100 fresh eligible rows: **299**
- Panel N terminal species: **200**
- Panel N terminal rows: **20,000**
- combined terminal species: **400**
- combined terminal rows: **40,000**
- metadata request errors: **0 / 600**
- previous observation IDs excluded: **228,362**
- previous photo IDs excluded: **228,362**
- queue extension: forbidden
- queue reordering: forbidden
- metadata rerun: forbidden
- species replacement after pixel opening: forbidden

The terminal 200 species in each panel are the first 200 full-capacity species in the already frozen 300-species queue. No biological outcome entered this choice.

## Outcome firewall at freeze

All are false:

- image pixels opened;
- flower colour opened;
- morph opened;
- palette opened;
- D opened;
- q_white/W opened;
- spatial colour outcome opened;
- H3 predictors opened.

## Durable artifact

Original workflow run: **35798104597**

Artifact:

- ID: **10725735599**
- name: `fcp-v2-fresh-metadata-terminal-20260923`
- digest: `sha256:247c09eee70fe56e381536d96ce090293e11d9f2924f27f45274850076be40b0`

The original workflow job concluded failure only because the branch advanced while the one-shot metadata draw was running, so its final result commit was non-fast-forward. The metadata draw, durable-result validation and artifact upload had already completed successfully. No metadata request was rerun. The terminal receipt was recovered from that immutable artifact and committed separately.

## Important file fingerprints inside the artifact

- `authorized_metadata_400x100.csv.gz`: `5cedb8581dd5da60006df8e9c939c318e94d21ff536b06c61d40d587f8800e0f`
- `terminal_species_manifest.tsv`: `ea7a24d8a7abc7f1dcfa7e388acc2eef61307ba29a63f6828b823cfd99b50668`
- Panel P terminal species: `9781a277541273268f37d8fbb03e9ff019a68deebd3c89597dcc217338e36a9f`
- Panel N terminal species: `665d156d912cd2190a813c1c3c67affe72dd32b8cd496328beaa60a5fc6956a7`

## Next authorized gate

The next step may open image pixels **only for response-blind technical measurement**.

Before biological colour outcomes are allowed:

1. build blind measurement IDs from the frozen 40,000 rows;
2. freeze the 20-images/species heavy-counterfactual subset from metadata identity only;
3. run technical exposure/background/ROI diagnostics;
4. seal the complete technical table and technical-quality strata;
5. read back and hash-validate that seal.

No morph, palette, D, q_white/W or spatial colour outcome may be calculated before that seal exists.
