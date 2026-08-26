#!/usr/bin/env python3
"""Build an outcome-independent literature-attention diagnostic for the C/S freeze.

Attention is the number of v2.2 title/abstract records containing the exact accepted
canonical species name. C/S evidence labels, query membership, and old spatial labels
are not used. This is deliberately a conservative name-mention count: synonym-heavy
species can be undercounted, so it is a sensitivity covariate rather than a completeness
measure.
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import pandas as pd


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument('--dataset', required=True)
    p.add_argument('--records', required=True)
    p.add_argument('--out', required=True)
    p.add_argument('--qc-out', required=True)
    args = p.parse_args()

    species = pd.read_csv(args.dataset)
    records = pd.read_csv(args.records)
    text = (records['title'].fillna('').astype(str) + ' ' + records['abstract'].fillna('').astype(str)).str.lower()

    rows = []
    for name in species['canonical_name'].astype(str):
        pattern = r'(?<![a-z])' + re.escape(name.lower()) + r'(?![a-z])'
        count = int(text.str.contains(pattern, regex=True).sum())
        title_count = int(records['title'].fillna('').astype(str).str.lower().str.contains(pattern, regex=True).sum())
        rows.append({
            'canonical_name': name,
            'n_v22_exact_name_records': count,
            'n_v22_exact_name_title_records': title_count,
            'attention_definition': 'exact accepted canonical binomial mention in v2.2 title/abstract',
        })

    out = pd.DataFrame(rows).sort_values('canonical_name').reset_index(drop=True)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, index=False)
    qc = {
        'species': int(len(out)),
        'records_scanned': int(len(records)),
        'species_with_zero_exact_name_records': int((out.n_v22_exact_name_records == 0).sum()),
        'median_exact_name_records': float(out.n_v22_exact_name_records.median()),
        'max_exact_name_records': int(out.n_v22_exact_name_records.max()),
        'semantic_guard': 'Outcome-independent exact-name attention diagnostic; synonym-heavy taxa may be undercounted and zero is not absence of literature.',
    }
    Path(args.qc_out).write_text(json.dumps(qc, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(qc, indent=2))


if __name__ == '__main__':
    main()
