#!/usr/bin/env python3
"""Paired PHOTO purity diagnostic across identical sampled taxa, not a test of genetic polymorphism."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
from fcp_global_photo_atlas_20261020 import (
    read_original, group_stats, high_purity, geographic_units, STATES, MIN_GLOBAL,
    MIN_OBSERVERS, MIN_PER_UNIT, check,
)


def bootstrap_mean_ci(x, seed=20261022, iterations=999):
    vals=np.asarray(x,float)
    if len(vals)<20:return [None,None]
    rng=np.random.default_rng(seed)
    b=rng.integers(0,len(vals),size=(iterations,len(vals)))
    return [float(a) for a in np.quantile(vals[b].mean(axis=1),[.025,.975])]


def compare(d, pure, keys, minimum, observers, *, region_key=None):
    raw=group_stats(d,keys,minimum,observers)
    good=group_stats(pure,keys,minimum,observers)
    if good.empty:
        return [], {'n':0,'status':'no matched opportunity'}
    keep=[*keys,'n_photo','D_pair','p_white','white_discordance_share','n_states']
    a=raw[keep].rename(columns={x:'all_'+x for x in keep if x not in keys})
    b=good[keep].rename(columns={x:'pure_'+x for x in keep if x not in keys})
    m=a.merge(b,on=keys,how='inner',validate='one_to_one')
    check(len(m)>0,'paired purity no overlapping taxon units')
    if region_key:
        return m, {str(k):stats(g) for k,g in m.groupby(region_key, observed=True)}
    return m,stats(m)


def stats(m):
    diff=m.pure_D_pair-m.all_D_pair
    whitediff=m.pure_p_white-m.all_p_white
    a=m[['all_white_discordance_share','pure_white_discordance_share']].dropna()
    return {
      'n_paired_units':int(len(m)),
      'mean_all_D_on_exact_pure_eligible_taxa':float(m.all_D_pair.mean()),
      'mean_pure_D_on_same_taxa':float(m.pure_D_pair.mean()),
      'pure_minus_all_D_mean':float(diff.mean()),
      'pure_minus_all_D_species_unit_bootstrap_CI95':bootstrap_mean_ci(diff),
      'mean_photo_white_fraction_all':float(m.all_p_white.mean()),
      'mean_photo_white_fraction_pure':float(m.pure_p_white.mean()),
      'pure_minus_all_white_fraction':float(whitediff.mean()),
      'species_with_two_plus_photo_colour_labels_all':int((m.all_n_states>=2).sum()),
      'species_with_two_plus_photo_colour_labels_pure':int((m.pure_n_states>=2).sum()),
      'white_discordance_fraction_n_matched_nonzero':int(len(a)),
      'mean_white_discordance_fraction_all':float(a.all_white_discordance_share.mean()) if len(a) else None,
      'mean_white_discordance_fraction_pure':float(a.pure_white_discordance_share.mean()) if len(a) else None,
    }


def analyse(raw):
    pure=high_purity(raw)
    global_m,global_result=compare(raw,pure,['inat_taxon_id'],MIN_GLOBAL,MIN_OBSERVERS)
    bands,by_band=compare(geographic_units(raw),geographic_units(pure),
                          ['inat_taxon_id','abs_lat_band'],MIN_PER_UNIT,MIN_OBSERVERS,
                          region_key='abs_lat_band')
    counts={s:{'raw':int((raw.morph==s).sum()),'pure':int((pure.morph==s).sum())} for s in STATES}
    for d in counts.values():d['retention_fraction']=d['pure']/d['raw'] if d['raw'] else None
    return {
      'schema':'fcp_real_original_photo_purity_same_taxa_matched_v1',
      'original_source_authenticated_by_upstream_reader':True,
      'all_classifiable_photos':int(len(raw)),
      'high_purity_photos':int(len(pure)),
      'quality_rule':'primary assigned-group ROI fraction >=0.65; all competing groups <=0.20; arbitrary diagnostic, not botanical validation',
      'state_specific_quality_retention':counts,
      'same_species_paired_min40_classifiable_3observers':global_result,
      'same_species_and_lat_band_paired_min10_3observers':by_band,
      'high_purity_filter_cannot_independently_validate_true_genetic_flower_morph':True,
      'species_bootstrap_is_descriptive_not_photo_level_independent':True,
      'needs_independent_flowers_or_spectral_pigment_validation':True,
    }


def selftest():
    rows=[]
    for tid in range(2):
        for i in range(50):
            rows.append({'inat_taxon_id':tid,'photo_id':tid*50+i+1,
                'observer_id':i//2,'morph':'white' if i<25 else 'red_pink',
                'flower_fraction_white':.92 if i<25 else .02,
                'flower_fraction_red':.02 if i<25 else .92,
                'flower_fraction_pink':.0,'flower_fraction_magenta':.0,
                'flower_fraction_blue':.0,'flower_fraction_purple':.0,
                'flower_fraction_yellow':.0,'flower_fraction_orange':.0})
    d=pd.DataFrame(rows)
    r=analyse_counts_only(d)
    check(r['n_paired_units']==2,'paired taxa identities')
    check(abs(r['pure_minus_all_D_mean'])<1e-12,'quality retained identical observations')
    print('FCP_PURITY_PAIRED_SYNTHETIC_PASS')


def analyse_counts_only(d):
    pure=high_purity(d)
    _,res=compare(d,pure,['inat_taxon_id'],MIN_GLOBAL,MIN_OBSERVERS)
    return res


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--photos',type=Path)
    p.add_argument('--out',type=Path)
    p.add_argument('--self-test',action='store_true')
    args=p.parse_args()
    if args.self_test:selftest();return
    check(args.photos is not None and args.out is not None,'paths required')
    selftest()
    out=analyse(read_original(args.photos))
    args.out.mkdir(parents=True,exist_ok=True)
    (args.out/'same_species_purity_diagnostic_result.json').write_text(
      json.dumps(out,sort_keys=True,indent=2,allow_nan=False)+'\n')
    print('FCP_REAL_SAME_SPECIES_PURITY_AUDIT',json.dumps(out,sort_keys=True,allow_nan=False))

if __name__=='__main__':main()
