#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import urllib.parse
import urllib.request
from collections import defaultdict
from pathlib import Path

API='https://mangal.io/api/v2/'


def fetch_pages(endpoint:str,max_pages:int=200):
    rows=[]
    for page in range(1,max_pages+1):
        q=urllib.parse.urlencode({'count':1000,'page':page})
        req=urllib.request.Request(API+endpoint+'?'+q,headers={
            'User-Agent':'fcp-spatiotemporal-interaction-preflight/0.1',
            'Accept':'application/json'
        })
        with urllib.request.urlopen(req,timeout=90) as h:
            x=json.loads(h.read().decode('utf-8'))
        if isinstance(x,dict):
            x=x.get('data',[x]) if isinstance(x.get('data'),list) else [x]
        if not x:
            break
        rows.extend(x)
        if len(x)<1000:
            break
    return rows


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--out',type=Path,required=True)
    a=ap.parse_args()

    networks=fetch_pages('network')
    nodes=fetch_pages('node')
    interactions=fetch_pages('interaction')

    network_geo={str(n['id']) for n in networks if n.get('geom') is not None}
    network_date={str(n['id']) for n in networks if n.get('date') is not None}
    node_taxon={}
    node_name={}
    for n in nodes:
        nid=str(n['id'])
        tax=n.get('taxonomy')
        if isinstance(tax,dict) and tax.get('id') is not None:
            node_taxon[nid]=str(tax['id'])
            node_name[nid]=str(tax.get('name') or n.get('original_name') or '')
        else:
            node_taxon[nid]=None
            node_name[nid]=str(n.get('original_name') or '')

    # For each interaction type x role, track reference taxon -> networks.
    taxon_networks=defaultdict(lambda:defaultdict(set))
    interaction_counts=defaultdict(int)
    networks_by_type=defaultdict(set)
    missing_taxonomy=0
    for x in interactions:
        typ=str(x.get('type') or 'UNKNOWN')
        net=str(x.get('network_id'))
        if net not in network_geo:
            continue
        interaction_counts[typ]+=1
        networks_by_type[typ].add(net)
        for role,key in [('from','node_from'),('to','node_to')]:
            nid=str(x.get(key))
            tid=node_taxon.get(nid)
            if tid is None:
                missing_taxonomy+=1
                continue
            taxon_networks[(typ,role)][tid].add(net)

    rows=[]
    for (typ,role),mapping in sorted(taxon_networks.items()):
        repeated=[tid for tid,nets in mapping.items() if len(nets)>=3]
        repeated5=[tid for tid,nets in mapping.items() if len(nets)>=5]
        rows.append({
            'interaction_type':typ,
            'role':role,
            'interactions_in_georeferenced_networks':interaction_counts[typ],
            'georeferenced_networks':len(networks_by_type[typ]),
            'reference_taxa_total':len(mapping),
            'reference_taxa_in_ge3_networks':len(repeated),
            'reference_taxa_in_ge5_networks':len(repeated5),
            'coverage_gate':len(repeated)>=20
        })

    eligible=[r for r in rows if r['coverage_gate']]
    ranked=sorted(
        eligible,
        key=lambda r:(-r['reference_taxa_in_ge3_networks'],-r['georeferenced_networks'],r['interaction_type'],r['role'])
    )
    admitted=ranked[:3]
    # Fixed primary families, frozen before any turnover outcome.
    primary_specs=[
        {'interaction_type':'mutualism','focal_role':'either'},
        {'interaction_type':'predation','focal_role':'from'},
        {'interaction_type':'parasitism','focal_role':'from'},
    ]
    primary=[]
    for spec in primary_specs:
        typ=spec['interaction_type']; role=spec['focal_role']
        if role=='either':
            merged=defaultdict(set)
            for side in ('from','to'):
                for tid,nets in taxon_networks.get((typ,side),{}).items():
                    merged[tid].update(nets)
            mapping=merged
        else:
            mapping=taxon_networks.get((typ,role),{})
        repeated=[tid for tid,nets in mapping.items() if len(nets)>=3]
        primary.append({
            'interaction_type':typ,
            'focal_role':role,
            'reference_taxa_total':len(mapping),
            'reference_taxa_in_ge3_networks':len(repeated),
            'coverage_gate':len(repeated)>=20,
        })

    result={
      'version':'v0.1',
      'status':'MANGAL_INTERACTION_ROLE_COVERAGE_PASS' if admitted else 'HOLD_MANGAL_NO_TYPE_ROLE_PASSES_COVERAGE',
      'network_rows':len(networks),
      'georeferenced_networks':len(network_geo),
      'dated_networks':len(network_date),
      'node_rows':len(nodes),
      'interaction_rows':len(interactions),
      'interaction_rows_missing_reference_taxonomy_across_roles':missing_taxonomy,
      'coverage_rule':{
        'minimum_reference_taxa_each_in_at_least_networks':3,
        'minimum_repeated_reference_taxa':20,
        'max_admitted_type_role_families':3
      },
      'type_role_coverage':rows,
      'admitted_candidate_type_role_families':admitted,
      'primary_family_coverage':primary,
      'all_primary_families_pass':all(x['coverage_gate'] for x in primary),
      'partner_similarity_computed':False,
      'rewiring_computed':False,
      'temporal_memory_computed':False,
      'selection_uses_turnover_outcomes':False
    }
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps(result,indent=2,sort_keys=True))


if __name__=='__main__':
    main()
