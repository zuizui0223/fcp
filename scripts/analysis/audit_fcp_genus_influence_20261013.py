#!/usr/bin/env python3
"""FCP source-exposed posthoc: taxonomic influence of E90-180 white-blue BIO12 slope.

This evaluates sensitivity of the already selected east white-vs-blue/purple
photographic slope difference to source genus prevalence and omission; it is not
new prospective confirmation, genetic polymorphism, pigment or adaptive fitness.
"""
from __future__ import annotations
import argparse, hashlib, io, json, zipfile
from pathlib import Path
import numpy as np
import pandas as pd

SHA='2d6097ba1b0434632e13bbe79dd78180608a287ff5dabebc770b5c0a2cfb52a0'
COLS=('wc_bio12','wc_bio1','wc_bio5','wc_bio15','wc_elevation_m','abs_lat','lon_sin','lon_cos')
HUES=('white','yellow_orange','red_pink','blue_purple')
TARGET='white_minus_blue_purple'


def ingest(zip_path: Path):
    if hashlib.sha256(zip_path.read_bytes()).hexdigest()!=SHA:
        raise ValueError('original FCP photo bytes changed')
    with zipfile.ZipFile(zip_path) as z:
        loc='fcp42111_real_site_abiotic/'
        a=pd.read_csv(io.BytesIO(z.read(loc+'breadth_flower_colour_climate_soil.csv.gz')),compression='gzip',low_memory=False)
        b=pd.read_csv(io.BytesIO(z.read(loc+'taxon_cell_flower_colour_climate_soil.csv.gz')),compression='gzip',low_memory=False)
    if len(a)!=42111 or len(b)!=85337:raise ValueError('panel counts drift')
    f=pd.concat((b,a),ignore_index=True)
    if f.photo_id.duplicated().sum()!=26905:raise ValueError('overlap drift')
    f=f.drop_duplicates('photo_id')
    if len(f)!=100543 or f.inat_taxon_id.nunique()!=42111:raise ValueError('source denominator drift')
    if f.measurement_status.eq('classified_four_state_morph').sum()!=45953:
        raise ValueError('classified source count drift')
    required=['wc_bio1','wc_bio5','wc_bio12','wc_bio15','wc_elevation_m','latitude','longitude']
    f=f.loc[f.measurement_status.eq('classified_four_state_morph') & f.morph.isin(HUES) & f[required].notna().all(axis=1) & f.latitude.between(-90,90) & f.longitude.between(-180,180)].copy()
    f['genus']=f.species.astype(str).str.split().str[0]
    f['abs_lat']=f.latitude.abs()
    f['lon_sin']=np.sin(np.deg2rad(f.longitude))
    f['lon_cos']=np.cos(np.deg2rad(f.longitude))
    for h in HUES:f[h]=f.morph.eq(h).astype(float)
    if len(f)!=45836 or f.inat_taxon_id.nunique()!=25028:raise ValueError('source sample drift')
    mu=f[list(COLS)].mean().to_numpy(float)
    sd=f[list(COLS)].std(ddof=0).to_numpy(float)
    f['east']=f.longitude.between(90,180,inclusive='both')
    return f,mu,sd


def prepare(case,mu,sd):
    cnt=case.groupby(['genus','east']).inat_taxon_id.nunique().unstack(fill_value=0)
    shared=cnt.index[(cnt[True]>=3)&(cnt[False]>=3)]
    f=case.loc[case.genus.isin(shared)].copy()
    if len(shared)<30:raise ValueError('min shared genera gate')
    out={}
    for region in (True,False):
        s=f.loc[f.east==region]
        agg={c:'mean' for c in COLS+HUES}
        agg.update({'genus':'first','photo_id':'size'})
        m=s.groupby('inat_taxon_id',as_index=False).agg(agg)
        m[list(COLS)]=(m[list(COLS)].to_numpy(float)-mu)/sd
        stats={}
        for genus,g in m.groupby('genus',sort=True):
            x=g[list(COLS)].to_numpy(float)
            y=g[list(HUES)].to_numpy(float)
            x-=x.mean(axis=0);y-=y.mean(axis=0)
            stats[genus]=(x.T@x,x.T@y,len(g),int(g.photo_id.sum()))
        out[region]=stats
    assert set(out[True])==set(out[False])==set(shared)
    return f,out


def whole_beta(stat,weight='species_equal',excluded=()):
    excluded=set(excluded)
    selected=sorted(set(stat[True])-excluded)
    if len(selected)<30:return None
    vals=[]
    for region in (True,False):
        xsum=np.zeros((8,8));ysum=np.zeros((8,4))
        for gen in selected:
            xx,xy,n,_=stat[region][gen]
            w=1.0 if weight=='species_equal' else 1.0/n
            xsum+=xx*w;ysum+=xy*w
        vals.append(np.linalg.solve(xsum,ysum)[0,:])
    a,b=vals
    w=HUES.index('white');blu=HUES.index('blue_purple')
    return float((a[w]-a[blu])-(b[w]-b[blu]))


def richness_weight_path(stats):
    """Descriptive bridge from species-equal (alpha=0) to genus-equal (alpha=1).

    Intermediate alpha is a sensitivity family, not a newly selected estimand.
    """
    out={}
    for alpha in (0.,0.25,0.50,0.75,1.):
        b=[]
        for region in (True,False):
            xx=np.zeros((8,8)); yy=np.zeros((8,4))
            for mat, rhs, n, photos in stats[region].values():
                w=n**(-alpha)
                xx+=w*mat;yy+=w*rhs
            b.append(np.linalg.solve(xx,yy)[0])
        e,o=b
        out[str(alpha)]=float(100*((e[0]-e[3])-(o[0]-o[3])))
    return out


def run(path:Path):
    f,mu,sd=ingest(path)
    strict=f.loc[(f.flower_effective_pixels>=f.flower_effective_pixels.median()) & (f.nuisance_pixel_fraction<=f.nuisance_pixel_fraction.median())].copy()
    samples={'all':f,'strict':strict}
    results={}
    ledger=[]
    for sample,frame in samples.items():
        sub,stats=prepare(frame,mu,sd)
        gen=sorted(stats[True]); meta=[]
        for g in gen:
            t=stats[True][g]; b=stats[False][g]
            meta.append({'genus':g,'n_east_species':t[2],'n_other_species':b[2],
                         'n_total_species_units':t[2]+b[2],
                         'n_east_photos':t[3],'n_other_photos':b[3]})
        meta=pd.DataFrame(meta)
        bridge=richness_weight_path(stats)
        for weight in ('species_equal','genus_equal'):
            full=whole_beta(stats,weight)
            vals=[]
            for rec in meta.itertuples(index=False):
                left=whole_beta(stats,weight,(rec.genus,))
                vals.append({'source_sample':sample,'weight':weight,'genus':rec.genus,
                            'n_east_species':int(rec.n_east_species),'n_other_species':int(rec.n_other_species),
                            'n_total_species_units':int(rec.n_total_species_units),
                            'n_east_photos':int(rec.n_east_photos),'n_other_photos':int(rec.n_other_photos),
                            'full_C_pp':100*full,'omit_C_pp':100*left,
                            'influence_delta_pp':100*(full-left)})
            influence=pd.DataFrame(vals)
            ledger.extend(vals)
            sorted_by_size=meta.sort_values(['n_total_species_units','genus'],ascending=[False,True]).genus.tolist()
            sorted_by_influence=influence.assign(absv=lambda d:d.influence_delta_pp.abs()).sort_values(['absv','genus'],ascending=[False,True]).genus.tolist()
            size_removed={str(k):100*whole_beta(stats,weight,sorted_by_size[:k]) for k in (1,5,10,20) if len(gen)>k+30}
            impact_removed={str(k):100*whole_beta(stats,weight,sorted_by_influence[:k]) for k in (1,5,10,20) if len(gen)>k+30}
            # Descriptive ranks / noncausal correlations only.
            top=influence.assign(abs_influence=lambda x:x.influence_delta_pp.abs()).sort_values('abs_influence',ascending=False).head(15)
            rho=float(influence[['n_total_species_units','influence_delta_pp']].corr(method='spearman').iloc[0,1]);
            n=meta.n_total_species_units.to_numpy(dtype=float)
            results[f'{sample}_{weight}']={
              'n_shared_genera':len(gen),'n_shared_source_photos':len(sub),'n_source_species':int(sub.inat_taxon_id.nunique()),
              'full_C_pp_per_sd':100*full,
              'sampled_species_richness_weight_sensitivity_pp':bridge,
              'leave_one_min_pp':float(influence.omit_C_pp.min()),'leave_one_max_pp':float(influence.omit_C_pp.max()),
              'leave_one_nonpositive_n':int((influence.omit_C_pp<=0).sum()),
              'left_positive_fraction':float((influence.omit_C_pp>0).mean()),
              'top_genus_by_absolute_impact':top[['genus','n_east_species','n_other_species','n_total_species_units','influence_delta_pp','omit_C_pp']].to_dict(orient='records'),
              'predefined_remove_largest_species_rich_genera_C_pp':size_removed,
              'posthoc_remove_top_abs_influence_C_pp':impact_removed,
              'impact_pp_min_max':[float(influence.influence_delta_pp.min()),float(influence.influence_delta_pp.max())],
              'spearman_genus_sampled_species_richness_vs_signed_influence':rho,
              'source_species_units_per_genus':{'min':int(n.min()),'median':float(np.median(n)),'p90':float(np.quantile(n,.90)),'max':int(n.max())},
              'n_genus_positive_signed_influence':int((influence.influence_delta_pp>0).sum()),
              'n_genus_negative_signed_influence':int((influence.influence_delta_pp<0).sum()),
              'caution':'Genus omission re-fits pooled coefficient and does not identify that genus own flower-color fitness; outcome-ranked removal is posthoc.'}
    out={'schema':'fcp_genus_influence_on_posthoc_white_blue_region_rain_v1','status':'ORIGINAL_FCP_100543_PHOTO_SOURCE_EXECUTED','source_sha256':SHA,
         'primary_90_180_region_already_selected':True,'full_model_8_covariates':COLS,
         'n_classified_original_photo_env':len(f),
         'benchmark_parent_all_species_equal_C_pp':7.62803895744896,
         'benchmark_parent_all_genus_equal_C_pp':5.309348222172678,
         'results':results,
         'hard_nonclaims':['Genus leave-one-out sensitivity is not proof of genus-local rainfall colour responses','No botanical organ verification, actual pigment cost, adaptation, phylogenetic replication, or independent new data','Preselected eastern sector and white-blue pair remain outcome-exposed','Classification missingness / photo exposure / actual environmental overlap remain unresolved']}
    return out,pd.DataFrame(ledger)


def main():
    p=argparse.ArgumentParser();p.add_argument('--source-archive',type=Path,required=True)
    p.add_argument('--outdir',type=Path,required=True)
    a=p.parse_args();z,ledger=run(a.source_archive)
    a.outdir.mkdir(parents=True,exist_ok=True)
    (a.outdir/'result.json').write_text(json.dumps(z,ensure_ascii=False,indent=2)+'\n')
    ledger.to_csv(a.outdir/'genus_influence.csv',index=False)
    print(json.dumps({k:{v:z['results'][k][v] for v in ('n_shared_genera','full_C_pp_per_sd','leave_one_min_pp','leave_one_max_pp','leave_one_nonpositive_n','predefined_remove_largest_species_rich_genera_C_pp')} for k in z['results']},indent=2))

if __name__=='__main__':main()
