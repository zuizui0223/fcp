#!/usr/bin/env python3
"""Compare all frozen FCP taxon IDs with real mutable Open Data taxonomy.

This is an identity/eligibility audit only; it neither opens photo pixels nor
uses any historical flower colour or H2 result. A mutable root object can NEVER
be the authoritative fixed-date source of a biological sample.
"""
from __future__ import annotations
import argparse
import csv
import json
from collections import Counter
from pathlib import Path

from fcp_pipeline.fcp_allphotos import opportunity_frame, sha256_file
from scripts.acquisition.import_inat_open_data_snapshot import (
    AUTHORITATIVE_42111_SHA256, read_tsv, integer, lineage_source_taxon,
)


def audit(frame_path: Path, taxa_path: Path, output_dir: Path) -> dict:
    if sha256_file(frame_path) != AUTHORITATIVE_42111_SHA256:
        raise ValueError('Original 42,111 source taxon identity is not byte-exact')
    frame = opportunity_frame(frame_path, expected=42111)
    original = dict(zip(frame.inat_taxon_id.astype(int), frame.species.astype(str)))
    valid = set(original)
    exact: dict[int, dict] = {}
    descendants = Counter()
    n_taxa = n_covered = 0
    for row in read_tsv(taxa_path, {'taxon_id', 'ancestry', 'rank', 'name', 'active'}):
        n_taxa += 1
        tid = integer(row['taxon_id'])
        if tid is None:
            continue
        included_taxon = lineage_source_taxon(tid, row['ancestry'], valid)
        if included_taxon is not None:
            n_covered += 1
            descendants[included_taxon] += 1
        if tid in valid:
            if tid in exact:
                raise ValueError(f'Repeated included taxon ID in provider: {tid}')
            exact[tid] = {'provider_name': row['name'], 'provider_rank': row['rank'],
                          'provider_active': row['active']}
    statuses = []
    name_mismatch = []
    missing = []
    rank_mismatch = []
    for tid, src_name in original.items():
        rec = exact.get(tid)
        if rec is None:
            missing.append({'inat_taxon_id': tid, 'fcp_name': src_name})
            status = 'taxon_id_absent_in_provider'
        else:
            if str(rec['provider_name']).strip() != src_name.strip():
                name_mismatch.append({'inat_taxon_id': tid, 'fcp_name': src_name,
                                      'provider_name': rec['provider_name']})
            if rec['provider_rank'].strip().lower() != 'species':
                rank_mismatch.append({'inat_taxon_id': tid, 'fcp_name': src_name,
                                      'provider_rank': rec['provider_rank']})
            status = 'exact_id_present'
        statuses.append({'inat_taxon_id': tid, 'fcp_name': src_name, 'status': status,
                         'provider_name': rec['provider_name'] if rec else '',
                         'provider_rank': rec['provider_rank'] if rec else '',
                         'provider_active': rec['provider_active'] if rec else '',
                         'provider_taxonomic_descendant_nodes': descendants[tid]})
    output_dir.mkdir(parents=True, exist_ok=True)
    with (output_dir/'fcp_taxa_crosswalk_42111.csv').open('w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f,fieldnames=statuses[0].keys());writer.writeheader();writer.writerows(statuses)
    report = {
        'schema': 'fcp_fresh_original_taxonomy_real_provider_audit_v1',
        'role': 'SCHEMA_AND_TAXONOMY_AUDIT_ONLY_NOT_A_FROZEN_SAMPLE',
        'photo_pixels_opened': False, 'old_flower_colour_measurements_read': False,
        'colour_distribution_tested': False, 'source_alias_mutable': True,
        'frame_species': len(frame), 'frame_sha256': sha256_file(frame_path),
        'current_provider_taxa_sha256': sha256_file(taxa_path),
        'provider_taxon_rows': n_taxa, 'provider_taxon_rows_mapped_to_frame': n_covered,
        'source_taxon_ids_present': len(exact),
        'source_taxon_ids_missing': len(missing),
        'source_name_mismatches': len(name_mismatch),
        'source_rank_mismatches': len(rank_mismatch),
        'species_with_mapped_descendant_taxa': sum(v>0 for v in descendants.values()),
        'source_names_unresolved_sample': missing[:15],
        'source_name_mismatch_sample': name_mismatch[:15],
        'source_rank_mismatch_sample': rank_mismatch[:15],
        'ready_for_full_photo_pixel_measurement': False,
    }
    (output_dir/'taxonomy_source_audit.json').write_text(json.dumps(report,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    print(json.dumps(report,ensure_ascii=False,sort_keys=True))
    return report


def main(argv=None)->int:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--frame',required=True,type=Path)
    p.add_argument('--taxa',required=True,type=Path)
    p.add_argument('--outdir',required=True,type=Path)
    a=p.parse_args(argv)
    audit(a.frame,a.taxa,a.outdir)
    return 0

if __name__=='__main__':
    raise SystemExit(main())
