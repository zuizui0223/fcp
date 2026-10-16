#!/usr/bin/env python3
"""Recompute staged NO-PIXEL workloads from the real immutable FCP source counts.

Counts and identity audit must be independently checked; this tool neither has
source image IDs nor opens original pixels, old flower-colour tables or environment.
Depths are resource-processing frontiers, NOT species/photo inclusion ceilings.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

COUNT_SHA = 'ef26c16252dd785f2d08bde3b18938c7624d73a54386fa3d3974dceb09ed97a9'
AUDIT_SHA = 'ed62aee0f9730ad17347d25180575af1cd101e0fc35c581e470c4a17866ebf56'
ARCHIVE_SHA = 'bfc1a97c436bf5989fccd43f12392953128489b76b92cea24da335ad73da4d8e'
SPECIES = 42111
PHOTOS = 93209351
DEPTHS = (1,10,20,40,64,100,128,256,512,1024,'ALL')
HEADERS = ('inat_taxon_id','species','n_unique_source_photo_ids','snapshot_date')


def sha256(path: Path) -> str:
    hasher = hashlib.sha256()
    with path.open('rb') as f:
        for part in iter(lambda: f.read(8 * 1024 * 1024), b''):
            hasher.update(part)
    return hasher.hexdigest()


def verify(counts: Path, audit: Path, budget: Path, *,
           _test_only_small: bool = False) -> dict:
    hashes = (sha256(counts),sha256(audit))
    if not _test_only_small and hashes != (COUNT_SHA,AUDIT_SHA):
        raise ValueError('Original 42111 source-count or independent audit SHA mismatch')
    a=json.loads(audit.read_text(encoding='utf-8'))
    r=json.loads(budget.read_text(encoding='utf-8'))
    if (a.get('status')!='VERIFIED_COMPLETE_COUNTS_ONLY_NOT_COLOUR_MEASUREMENT'
        or a.get('original_source_taxon_IDs_and_names_independently_matched') is not True
        or a.get('snapshot_date')!='2026-09-27'
        or a.get('all_source_photo_id_manifests_frozen') is not False
        or a.get('all_original_image_pixels_measured') is not False
        or r.get('schema')!='fcp_42111_verified_original_photo_execution_waves_v1'
        or r.get('status')!='VERIFIED_SOURCE_OPPORTUNITY_WAVES_NOT_PHOTO_IDS_OR_PIXELS'
        or r.get('n_original_new_image_pixels_measured_for_full_20260927_census')!=0
        or r.get('all_remaining_photos_still_eligible') is not True
        or r.get('no_color_D_W_H2_or_environment_results_from_counts') is not True
        or r.get('photo_fraction_in_top_one_percent_species') is None
        or r.get('source_photo_opportunity_gini') is None):
        raise ValueError('Count-only FCP scientific boundary was not preserved')
    if not _test_only_small and (
        a.get('source_archive_sha256')!=ARCHIVE_SHA
        or a.get('source_frame_species')!=SPECIES
        or a.get('all_source_photo_opportunity_ids')!=PHOTOS
        or a.get('count_ledger_sha256')!=COUNT_SHA):
        raise ValueError('Real 2026-09-27 provider source audit changed')
    n=[]
    seen=set()
    with counts.open(newline='',encoding='utf-8') as f:
        reader=csv.DictReader(f)
        if tuple(reader.fieldnames or ())!=HEADERS:
            raise ValueError('Original source photo count header changed')
        for x in reader:
            tid=str(x['inat_taxon_id']);count=str(x['n_unique_source_photo_ids'])
            if (not tid.isascii() or not tid.isdecimal() or tid in seen
                or not count.isascii() or not count.isdecimal()
                or not x['species'] or x['snapshot_date']!='2026-09-27'):
                raise ValueError('Duplicate, changed, or malformed original-source taxon')
            seen.add(tid)
            n.append(int(count))
    expected_taxa=a['source_frame_species']
    expected_photos=a['all_source_photo_opportunity_ids']
    if (len(n)!=expected_taxa or sum(n)!=expected_photos or
        r['source_species']!=expected_taxa or
        r['source_photo_opportunities']!=expected_photos or
        r['source_count_ledger_sha256']!=hashes[0] or
        r['independently_verified_audit_sha256']!=hashes[1]):
        raise ValueError('Frozen 42111 full count/source wave denominator drifted')
    if (sum(x>0 for x in n)!=a['taxa_with_at_least_one_photo']
        or n.count(0)!=a['taxa_with_zero_photos']):
        raise ValueError('Original zero-photo and positive-photo taxa changed')
    top=math.ceil(len(n)/100)
    ordered=sorted(n)
    gini=sum((2*i-len(n)-1)*v for i,v in enumerate(ordered,1))/(len(n)*sum(n))
    top_share=sum(ordered[-top:])/sum(n)
    if (abs(gini-r['source_photo_opportunity_gini'])>1e-8 or
        abs(top_share-r['photo_fraction_in_top_one_percent_species'])>1e-10):
        raise ValueError('Independent source opportunity concentration drifted')
    stages=r['stages']
    if len(stages)!=len(DEPTHS):
        raise ValueError('Missing full-study staged workload incl uncapped final pass')
    prev=0
    for step,depth in zip(stages,DEPTHS):
        cumulative=(sum(min(z,depth) for z in n) if depth!='ALL' else sum(n))
        if (step.get('photo_source_depth_work_frontier')!=depth
            or step.get('cumulative_opportunities')!=cumulative
            or step.get('incremental_opportunities')!=cumulative-prev
            or abs(step.get('fraction_of_all_original_opportunities',-1)-cumulative/sum(n))>1e-9
            or step.get('real_pixels_measured')!=0):
            raise ValueError('Frozen source workload stage lost actual original photos')
        prev=cumulative
    if prev!=sum(n) or stages[-1]['incremental_opportunities']<0:
        raise ValueError('Final ALL original-photo stage must include every photo')
    return {'status':'ORIGINAL_42111_REAL_PHOTO_OPPORTUNITY_WORKLOAD_INDEPENDENTLY_REPRODUCED_NO_PIXELS',
            'n_species':len(n),'n_original_source_photo_opportunities':sum(n),
            'first_40_depth_opportunities':stages[3]['cumulative_opportunities'],
            'first_128_depth_opportunities':stages[6]['cumulative_opportunities'],
            'uncapped_final_all_opportunities':stages[-1]['cumulative_opportunities'],
            'biological_measurements_performed':0}


def main(argv=None) -> int:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--counts',type=Path,required=True)
    p.add_argument('--audit',type=Path,required=True)
    p.add_argument('--budget',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    x=p.parse_args(argv)
    checked=verify(x.counts,x.audit,x.budget)
    if x.out.exists():
        raise ValueError('Do not overwrite previous source opportunity verification')
    x.out.parent.mkdir(parents=True,exist_ok=True)
    x.out.write_text(json.dumps(checked,indent=2,sort_keys=True)+'\n',encoding='utf-8')
    print('FCP_REAL_SOURCE_OPPORTUNITY_WAVE_AUDIT',json.dumps(checked,sort_keys=True))
    return 0


if __name__=='__main__':
    raise SystemExit(main())
