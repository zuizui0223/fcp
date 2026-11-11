#!/usr/bin/env python3
"""Test whether ~47% photo-white in broad latitudes is unusually even.

Fixed classifiable source photo opportunities. Two separate neutral label shuffles
conditionally preserve all observed genus (or species) white/photo colour totals.
This can reject fine-scale geographic constancy as an exceptional phenotype
pattern but cannot infer adaptation or label accuracy.
"""
from __future__ import annotations
import argparse,json,sys
from pathlib import Path
import numpy as np
import pandas as pd

sys.path.insert(0,str(Path(__file__).resolve().parent))
from audit_fcp_global_photo_rainfall_classifiability_IPW_MNAR_20261010 import load_original,SOURCE_SHA

REPS=999
SEED=202610102257

def uniformity_test(source,*,permutations=REPS):
    s=source.loc[source.classified&source.latitude.between(-90,90)&source.longitude.between(-180,180)].copy()
    lat=np.abs(s.latitude.to_numpy(float))
    band=np.searchsorted([30.,60.],lat,side='right')
    if not (np.bincount(band,minlength=3)>100).all():raise RuntimeError('latitudinal opportunity underidentified')
    w=s.morph.eq('white').to_numpy(np.int8)
    totals=np.bincount(band,minlength=3)
    p_obs=np.bincount(band,weights=w,minlength=3)/totals
    observed=float(p_obs.max()-p_obs.min())
    results={}
    rng=np.random.default_rng(SEED)
    for label,values in [('within_species',s.inat_taxon_id),('within_genus',s.genus)]:
        codes,_=pd.factorize(values,sort=True)
        idx=np.argsort(codes,kind='stable')
        original_counts=pd.DataFrame({'group':codes,'w':w}).groupby('group',sort=True)['w'].sum().to_numpy(int)
        nullspreads=np.empty(permutations,float)
        for i in range(permutations):
            order=np.lexsort((rng.random(len(codes)),codes))
            perm=np.empty_like(w)
            perm[idx]=w[order]
            null_counts=np.bincount(codes,weights=perm,minlength=len(original_counts)).astype(int)
            if not np.array_equal(null_counts,original_counts):
                raise RuntimeError('source colour-composition preservation violated')
            p=np.bincount(band,weights=perm,minlength=3)/totals
            nullspreads[i]=p.max()-p.min()
        results[label]={
           'n_groups':len(original_counts),
           'n_photo_labels_shuffled':len(w),
           'group_white_totals_preserved_exactly':True,
           'null_mean_range_percentage_points':float(100*nullspreads.mean()),
           'null_95pct_range_percentage_points':[float(x) for x in 100*np.quantile(nullspreads,[.025,.975])],
           'one_sided_p_observed_unusually_even_or_lower_range':float((1+int((nullspreads<=observed+1e-12).sum()))/(permutations+1)),
           'one_sided_p_observed_unusually_uneven_or_greater_range':float((1+int((nullspreads>=observed-1e-12).sum()))/(permutations+1)),
           'null_sample_count':permutations,
        }
    return {
       'schema':'fcp_original100543_global_photo_white_broad_latitude_constancy_conditional_null_v1',
       'source_sha256':SOURCE_SHA,
       'n_unique_original_photos_all_states':len(source),
       'n_original_classified_geotagged_photos':len(s),
       'original_latitude_band_n_photos':totals.astype(int).tolist(),
       'photo_white_fractions_by_original_true_site_latitude_band':p_obs.tolist(),
       'observed_range_white_percentage_points':100*observed,
       'latitudinal_group_band_degrees_abs':[[0,30],[30,60],[60,90]],
       'conditional_colour_photo_label_nulls':results,
       'biological_interpretation':'A small raw latitudinal spread is not a biological 50/50 pigment-cost equilibrium; matched colour composition and photo geographic opportunity are fixed in each photo label null.',
       'alternative_photo_colour_unknown_count':int((~source.classified).sum()),
       'source_flower_colour_human_truth_unmeasured':True,
       'latent_adaptive_selection_equilibrium_unidentified':True,
    }

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--source-archive',required=True,type=Path);ap.add_argument('--outdir',required=True,type=Path)
    args=ap.parse_args()
    src=load_original(args.source_archive)
    a=uniformity_test(src)
    args.outdir.mkdir(parents=True,exist_ok=True)
    (args.outdir/'latitude_uniformity_result.json').write_text(json.dumps(a,indent=2)+'\n')
    print(json.dumps(a,indent=2))