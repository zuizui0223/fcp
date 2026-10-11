#!/usr/bin/env python3
"""Source-locked FCP genus/geography ecological transport of rainfall–flower photo colour.

Sample all 100543 unique original photographed opportunities; only original
four-hue classifiable rows have 4-category labels. Separately use original
nonzero nine-hue floral ROI pixels, including ambiguous palettes. Genre-normal
within-genus regression is NOT true phylogeny or fitness/benefit-cost selection.
"""
from __future__ import annotations
import argparse,hashlib,io,json,zipfile
from pathlib import Path
import numpy as np,pandas as pd
from scipy.stats import norm

SHA="2d6097ba1b0434632e13bbe79dd78180608a287ff5dabebc770b5c0a2cfb52a0"
HUES=("white","yellow_orange","red_pink","blue_purple")
COLS=("wc_bio12","wc_bio1","wc_bio5","wc_bio15","wc_elevation_m","abs_lat","lon_sin","lon_cos")
SEED=20261011

def get_source(p:Path):
    h=hashlib.sha256(p.read_bytes()).hexdigest()
    if h!=SHA:raise RuntimeError("original 100543 photo source ZIP SHA256 mismatch")
    with zipfile.ZipFile(p) as z:
        prefix="fcp42111_real_site_abiotic/"
        def data(x):
            return pd.read_csv(io.BytesIO(z.read(prefix+x+".csv.gz")),compression="gzip",low_memory=False)
        cells=data("taxon_cell_flower_colour_climate_soil")
        broad=data("breadth_flower_colour_climate_soil")
    if len(cells)!=85337 or len(broad)!=42111:raise ValueError("original photo source frame changed")
    u=pd.concat([cells,broad],ignore_index=True)
    if u.photo_id.duplicated().sum()!=26905:raise ValueError("original photo overlap drift")
    original_class=u.loc[u.measurement_status.eq("classified_four_state_morph")].groupby("photo_id").morph.nunique()
    if (original_class>1).any():raise ValueError("photo hue conflict across source panels")
    u=u.drop_duplicates("photo_id")
    u["classifiable"]=u.measurement_status.eq("classified_four_state_morph") & u.morph.isin(HUES)
    if len(u)!=100543 or u.inat_taxon_id.nunique()!=42111 or int(u.classifiable.sum())!=45953:
        raise ValueError("original global 100543/42111/45953 counts drift")
    return u

def eligible(u):
    cols=list(COLS[:5])+["latitude","longitude"]
    v=u[u[cols].notna().all(axis=1)&u.latitude.between(-90,90)&u.longitude.between(-180,180)].copy()
    v["genus"]=v.species.astype(str).str.split().str[0]
    v["abs_lat"]=v.latitude.abs()
    v["lon_sin"]=np.sin(np.deg2rad(v.longitude))
    v["lon_cos"]=np.cos(np.deg2rad(v.longitude))
    v["sector"]=np.minimum(3,np.floor((v.longitude+180)/90).astype(int))
    v["zone"]=np.where(v.abs_lat<23.5,"tropical","extratropical")
    v["region"]=v.sector.astype(str)+"_"+v.zone
    cls=v[v.classifiable].copy()
    if len(cls)!=45836 or cls.inat_taxon_id.nunique()!=25028:raise ValueError("source classifiable photo climate changed")
    for c in HUES:cls[c]=cls.morph.eq(c).astype(float)
    mean=cls[list(COLS)].mean().to_numpy(float)
    scale=cls[list(COLS)].std(ddof=0).to_numpy(float)
    roi_cols=[x for x in u if x.startswith("flower_fraction_")]
    a=v.loc[v[roi_cols].notna().all(axis=1)&v[roi_cols].sum(axis=1).between(.999999,1.000001)].copy()
    a["white"]=a.flower_fraction_white
    a["yellow_orange"]=a[["flower_fraction_yellow","flower_fraction_orange","flower_fraction_bronze"]].sum(axis=1)
    a["red_pink"]=a[["flower_fraction_red","flower_fraction_pink","flower_fraction_magenta"]].sum(axis=1)
    a["blue_purple"]=a[["flower_fraction_blue","flower_fraction_purple"]].sum(axis=1)
    if len(a)!=37575 or int(a.classifiable.sum())!=34291 or not np.allclose(a[list(HUES)].sum(axis=1),1.):
        raise ValueError("nonzero original ROI hue share population changed")
    return v,cls,a,mean,scale

def species_means(rows,mu,sd,by_region=False):
    grouping=["inat_taxon_id"] if not by_region else ["region","inat_taxon_id"]
    agg={x:"mean" for x in list(COLS)+list(HUES)}
    agg.update({"photo_id":"size","genus":"first"})
    x=rows.groupby(grouping,as_index=False).agg(agg).rename(columns={"photo_id":"n_photos"})
    x[list(COLS)]=(x[list(COLS)].to_numpy(float)-mu)/sd
    return x

def fit(sub,equal_genus=False,return_blocks=False):
    # Each source species contributes one mean colour/environment observation.
    sizes=sub.genus.value_counts()
    d=sub[sub.genus.isin(sizes[sizes>=3].index)].copy()
    n=len(d);k=d.genus.nunique()
    if n<300 or k<30:return {"status":"HOLD_SMALL","n_species":n,"n_genera":k}
    genus,ix=np.unique(d.genus.to_numpy(str),return_inverse=True)
    x=d[list(COLS)].to_numpy(float)
    y=d[list(HUES)].to_numpy(float)
    w=np.ones(n) if not equal_genus else 1/sizes.reindex(d.genus).to_numpy(float)
    sz=np.bincount(ix,weights=w)
    gx=np.stack([np.bincount(ix,weights=w*x[:,j])/sz for j in range(8)],axis=1)
    gy=np.stack([np.bincount(ix,weights=w*y[:,j])/sz for j in range(4)],axis=1)
    x-=gx[ix];y-=gy[ix]
    a=x.T@(w[:,None]*x)
    if np.linalg.cond(a)>1e10:return {"status":"HOLD_RANK","n_species":n,"n_genera":k}
    ia=np.linalg.inv(a)
    beta=ia@(x.T@(w[:,None]*y))
    resid=y-x@beta
    hue=[]
    for j,c in enumerate(HUES):
        score=np.zeros((k,8));np.add.at(score,ix,x*(w*resid[:,j])[:,None])
        fac=k/(k-1)*(n-1)/max(1,n-8-k)
        se=float(np.sqrt(max(0,(ia@(score.T@score)@ia*fac)[0,0])))
        v=float(beta[0,j])
        hue.append({"hue":c,"coef":v,"CI95":[v-1.96*se,v+1.96*se],"p_normal":float(2*norm.sf(abs(v/se))) if se>0 else 1.})
    out={"status":"OK","n_species_or_taxon_region_units":n,"n_genera":k,"n_original_photo_records":int(d.n_photos.sum()),"genus_equal":equal_genus,"rainfall_hue":hue}
    if not return_blocks:return out
    blocks={}
    for j,g in enumerate(genus):
        ind=ix==j
        xx=x[ind];yy=y[ind,0];ww=w[ind]
        blocks[g]=(xx.T@(ww[:,None]*xx),xx.T@(ww*yy))
    return out,blocks

def compare_geography(original,mu,sd,boot=999):
    e=species_means(original[original.sector.eq(3)],mu,sd)
    o=species_means(original[original.sector.ne(3)],mu,sd)
    re,A=fit(e,return_blocks=True);ro,B=fit(o,return_blocks=True)
    gs=sorted(set(A)|set(B))
    ids={name:j for j,name in enumerate(gs)}
    n=len(gs)
    aa=np.zeros((n,8,8));bb=np.zeros((n,8,8))
    ay=np.zeros((n,8));by=np.zeros((n,8))
    for genus,(x,y) in A.items():aa[ids[genus]]=x;ay[ids[genus]]=y
    for genus,(x,y) in B.items():bb[ids[genus]]=x;by[ids[genus]]=y
    v=float(np.linalg.solve(aa.sum(axis=0),ay.sum(axis=0))[0]-np.linalg.solve(bb.sum(axis=0),by.sum(axis=0))[0])
    rng=np.random.default_rng(SEED)
    weights=rng.multinomial(n,[1/n]*n,size=boot)
    M1=np.einsum("rn,nij->rij",weights,aa,optimize=True)
    M2=np.einsum("rn,nij->rij",weights,bb,optimize=True)
    b1=np.linalg.solve(M1,(weights@ay)[:,:,None])[:,0,0]
    b2=np.linalg.solve(M2,(weights@by)[:,:,None])[:,0,0]
    diff=b1-b2
    return {"E_90_180_n_species":re["n_species_or_taxon_region_units"],"remaining_n_species":ro["n_species_or_taxon_region_units"],
      "n_global_source_genera_jointly_resampled":n,"bootstrap_repetitions":boot,
      "E_90_180_minus_other_BIO12_photo_white_slope":v,
      "CI95":np.quantile(diff,[.025,.975]).tolist(),
      "one_sided_p_difference_not_positive":float((1+(diff<=0).sum())/(boot+1)),
      "warning":"Longitude sector chosen after viewing original results; this p-value is exploratory, not out-of-sample selection-free."}

def run(path,boot=999):
    orig=get_source(path);v,c,a,mu,sd=eligible(orig)
    qarea=float(c.flower_effective_pixels.median());qnuis=float(c.nuisance_pixel_fraction.median())
    specs={"global_all_4state":c,"E90_180_4state":c[c.sector.eq(3)],
       "outside_E90_180_4state":c[c.sector.ne(3)],
       "global_atleast_two_classified_photos":c,
       "global_source_photo_quality_strict":c[(c.flower_effective_pixels>=qarea)&(c.nuisance_pixel_fraction<=qnuis)],
       "tropical_4state":c[c.zone.eq("tropical")],
       "extratropical_4state":c[c.zone.eq("extratropical")],
       "E90_180_nonzero_ROI":a[a.sector.eq(3)],
       "outside_E90_180_nonzero_ROI":a[a.sector.ne(3)],
       "global_nonzero_ROI":a,
       "E90_180_ROI_quality_strict":a[(a.sector.eq(3))&(a.flower_effective_pixels>=qarea)&(a.nuisance_pixel_fraction<=qnuis)]}
    for b in range(4):specs["drop_longitude_sector_"+str(b)]=c[c.sector.ne(b)]
    models={}
    for k,d in specs.items():
        f=species_means(d,mu,sd)
        if k=="global_atleast_two_classified_photos":f=f[f.n_photos>=2]
        models[k]={"species_equal":fit(f),"genus_equal":fit(f,equal_genus=True),
                   "eligible_named_species_before_genus_support_filter":int(f.inat_taxon_id.nunique())}
    geographically_partitioned={}
    f=species_means(c,mu,sd,by_region=True)
    for k,d in f.groupby("region"):
        geographically_partitioned[k]={"species_equal":fit(d),"genus_equal":fit(d,True),
                                       "n_species_region_units_before_genus_filter":len(d)}
    counts={"all_source_photo_ids":len(orig),"all_source_nominal_species":int(orig.inat_taxon_id.nunique()),
       "all_classifiable_photos":int(orig.classifiable.sum()),"all_unclassifiable_photos":int((~orig.classifiable).sum()),
       "original_climate_site_opportunity_photos":len(v),"four_state_climate_complete_photos":len(c),
       "source_distinct_species_with_four_state_and_climate":int(c.inat_taxon_id.nunique()),
       "nonzero_original_nine_hue_ROI_climate_complete_photos":len(a),
       "nonzero_ROI_ambiguous_coarse_colour_photos":int((~a.classifiable).sum()),
       "photo_source_rainfall_one_SD_mm":float(sd[0]),"photo_quality_area_median":qarea,
       "photo_quality_nuisance_median":qnuis}
    return {"schema":"fcp_full42111_original_photo_rain_hue_realm_generalization_v1",
       "date_jst":"2026-10-11","source_archive_sha256":SHA,"source":counts,
       "source_climate_features":list(COLS),"source_four_hues":list(HUES),
       "four_longitude_sectors":"fixed -180/-90/0/90/180; north and south are not necessarily independent genotypes",
       "geography":geographically_partitioned,"models":models,
       "E90_180_contrast_posthoc":compare_geography(c,mu,sd,boot),
       "nonclaims":["Photographic white/coloured class is not chemically validated anthocyanin pathway or independently confirmed petal.","Genus centering/cluster SE does not adjust actual dated branch lengths and true fine photo site geodesic spatial random field.","Original 54590 unclassifiable photos remain colour unknown; valid nine-hue ROI makes 3284 palette-ambiguous original photographs additionally evaluable only as continuous hues.","Strong-quality subset reduces support for observed rainfall association.","Posthoc geographic winners and P values require independent species/photo confirmation.","No photos supply pigment costs, pollinator/herbivore payoff, fitness selection or 50:50 within natural populations."]}

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source-archive",type=Path,required=True)
    p.add_argument("--outdir",type=Path,required=True)
    p.add_argument("--bootstrap",type=int,default=999)
    o=p.parse_args()
    result=run(o.source_archive,o.bootstrap)
    o.outdir.mkdir(parents=True,exist_ok=True)
    (o.outdir/"result.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps({"source":result["source"],"east_vs_rest":result["E90_180_contrast_posthoc"],
        "global":result["models"]["global_all_4state"]["species_equal"],"regions":{k:v["species_equal"]["rainfall_hue"][0] for k,v in result["geography"].items()}},indent=2))

if __name__=="__main__":main()
