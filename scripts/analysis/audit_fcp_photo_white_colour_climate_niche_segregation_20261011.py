#!/usr/bin/env python3
"""Within-species original photographed white/coloured environmental niche contrasts.

Retrospective original photo labels, exactly as previously classified; four
environment axes, three original-photo-depth gates, 999 per-species colour
shuffles; spatial 10/20 degree conditional shuffles preserve colour in each
species-by-geographic cell. Genetic morphs or fitness are not identified.
"""
from __future__ import annotations
import argparse,hashlib,io,json,zipfile
from pathlib import Path
import numpy as np,pandas as pd

SHA="2d6097ba1b0434632e13bbe79dd78180608a287ff5dabebc770b5c0a2cfb52a0"
SEED=2026101102
REPS=999
AXES=("wc_bio12","wc_bio15","wc_bio5","wc_elevation_m")
HUES=("white","yellow_orange","red_pink","blue_purple")
GATES=((5,2),(10,3),(10,5))

def load_original(path):
    path=Path(path)
    if hashlib.sha256(path.read_bytes()).hexdigest()!=SHA:
        raise RuntimeError("original full-photo source zip SHA256 mismatch")
    with zipfile.ZipFile(path) as archive:
        a=pd.read_csv(io.BytesIO(archive.read(
            "fcp42111_real_site_abiotic/taxon_cell_flower_colour_climate_soil.csv.gz")),
            compression="gzip",low_memory=False)
        b=pd.read_csv(io.BytesIO(archive.read(
            "fcp42111_real_site_abiotic/breadth_flower_colour_climate_soil.csv.gz")),
            compression="gzip",low_memory=False)
    if len(a)!=85337 or len(b)!=42111:raise RuntimeError("Original photo ledger changed")
    z=pd.concat([a,b],ignore_index=True)
    duplicates=z.loc[z.photo_id.duplicated(keep=False)]
    if duplicates.photo_id.nunique()!=26905 or (duplicates.groupby("photo_id").inat_taxon_id.nunique()>1).any():
        raise RuntimeError("original 26905 photo IDs or species correspondence changed")
    mixed=duplicates[duplicates.measurement_status.eq("classified_four_state_morph")].groupby("photo_id").morph.nunique()
    if len(mixed) and mixed.max()>1:raise RuntimeError("original photo colour conflict")
    z=z.drop_duplicates("photo_id")
    yes=z.measurement_status.eq("classified_four_state_morph") & z.morph.isin(HUES)
    if len(z)!=100543 or int(yes.sum())!=45953 or z.inat_taxon_id.nunique()!=42111:
        raise RuntimeError("original 100543/45953/42111 source populations changed")
    valid=yes&z.latitude.between(-90,90)&z.longitude.between(-180,180)&z[list(AXES)].notna().all(axis=1)
    d=z.loc[valid,["inat_taxon_id","species","photo_id","morph","latitude","longitude",*AXES]].copy()
    if len(d)!=45836:raise RuntimeError("original 45836 classified site-climate photos changed")
    d["white"]=d.morph.eq("white").astype(np.int8)
    d["genus"]=d.species.astype(str).str.split().str[0]
    return d.sort_values(["inat_taxon_id","photo_id"],kind="stable").reset_index(drop=True)

def geographic_cell(lat,lon,degrees):
    lat=np.asarray(lat,float);lon=np.asarray(lon,float)
    return (np.floor((lat+90)/degrees)*int(np.ceil(360/degrees)+1)+np.floor((lon+180)/degrees)).astype(int)

def select_species(d):
    n=d.groupby("inat_taxon_id").white.agg(["size","sum"])
    ids=n.index[(n["size"]>=5)&(n["sum"]>=2)&(n["size"]-n["sum"]>=2)]
    v=[]
    for taxon,g in d.loc[d.inat_taxon_id.isin(ids)].groupby("inat_taxon_id",sort=True):
        x=g[list(AXES)].to_numpy(float)
        sd=x.std(axis=0)
        norm=(x-x.mean(axis=0))/np.where(sd>1e-10,sd,1.)
        w=g.white.to_numpy(np.int8)
        v.append({"taxon":int(taxon),"genus":g.genus.iloc[0],
                  "n":len(g),"n_white":int(w.sum()),"labels":w,"z":norm,
                  "cell10":geographic_cell(g.latitude,g.longitude,10),
                  "cell20":geographic_cell(g.latitude,g.longitude,20)})
    if len(v)!=422 or sum(i["n"] for i in v)!=4465:
        raise RuntimeError("original source classifiable mixed species 422/4465 drift")
    return v

def contrast(z,labels):
    labels=np.asarray(labels,np.int8)
    if labels.ndim==1:labels=labels[None,:]
    w=labels.sum(axis=1);c=labels.shape[1]-w
    if np.any(w<1)|np.any(c<1):raise ValueError("both photo hue categories required")
    return (labels@z)/w[:,None] - ((1-labels)@z)/c[:,None]

def shuffle(labels,geo,rng,nrep):
    labels=np.asarray(labels,np.int8)
    x=np.broadcast_to(labels,(nrep,len(labels))).copy()
    group=[np.arange(len(labels))] if geo is None else [np.flatnonzero(geo==k) for k in np.unique(geo)]
    for idx in group:
        if len(idx)<=1:continue
        order=np.argsort(rng.random((nrep,len(idx))),axis=1)
        x[:,idx]=labels[idx][order]
        if not np.array_equal(np.sort(x[:,idx],axis=1),np.tile(np.sort(labels[idx]),(nrep,1))):
            raise RuntimeError("species/geographic block photo hue totals changed")
    if not np.array_equal(x.sum(axis=1),np.full(nrep,labels.sum())):
        raise RuntimeError("species original hue totals changed")
    return x

def source_spatial_exchange(v):
    out={}
    for degrees in (10,20):
        photo=0;blocks=0;species=set()
        for g in v:
            for c in np.unique(g["cell"+str(degrees)]):
                x=g["labels"][g["cell"+str(degrees)]==c]
                if len(x)>1 and 0<x.sum()<len(x):
                    photo+=len(x);blocks+=1;species.add(g["taxon"])
        out[str(degrees)+"degree"]={"n_species_with_both_original_photo_colours_in_one_geocell":len(species),
                                    "informative_source_photo_count":photo,"informative_source_blocks":blocks}
    return out

def evaluate(v,nrep=REPS):
    rng=np.random.default_rng(SEED)
    n=len(v)
    obs=np.empty((n,4),float)
    null={k:np.empty((nrep,n,4),float) for k in ("unrestricted","within_10degree","within_20degree")}
    for i,s in enumerate(v):
        obs[i]=contrast(s["z"],s["labels"])[0]
        for mode,loc in (("unrestricted",None),("within_10degree",s["cell10"]),("within_20degree",s["cell20"])):
            null[mode][:,i,:]=contrast(s["z"],shuffle(s["labels"],loc,rng,nrep))
    results={}
    for min_photo,min_each in GATES:
        ix=np.array([s["n"]>=min_photo and s["n_white"]>=min_each and s["n"]-s["n_white"]>=min_each for s in v])
        if ix.sum()<20:raise RuntimeError("insufficient original source colour-mixed species")
        res={"n_original_species":int(ix.sum()),"n_source_photo_rows":int(sum(s["n"] for s,b in zip(v,ix) if b)),
             "n_distinct_genera":len({s["genus"] for s,b in zip(v,ix) if b}),"environment":{}}
        for j,name in enumerate(AXES):
            x=obs[ix,j];signed=float(np.mean(x));square=float(np.mean(x*x));r={}
            for model,sims in null.items():
                perms=sims[:,ix,j]
                simsigned=np.mean(perms,axis=1)
                simsquares=np.mean(perms*perms,axis=1)
                r[model]={
                  "signed_white_minus_nonwhite_specieswise_environmental_SD":signed,
                  "mean_squared_species_white_nonwhite_environmental_contrast":square,
                  "null_expected_mean_squared_contrast":float(simsquares.mean()),
                  "observed_minus_null_squared_environmental_segregation":float(square-simsquares.mean()),
                  "p_signed_two_sided":float((1+np.count_nonzero(np.abs(simsigned)>=abs(signed)))/(nrep+1)),
                  "p_heterogeneous_species_niche_segregation_upper":float((1+np.count_nonzero(simsquares>=square))/(nrep+1)),
                  "null_95_interval":[float(z) for z in np.quantile(simsquares,[.025,.975])]
                }
            res["environment"][name]=r
        results[f"min{min_photo}_min{min_each}each"]=res
    for v1 in results.values():
        for axis in v1["environment"].values():
            for r in axis.values():
                r["p_heterogeneous_bonferroni_all_12_axis_depth"] = min(1.,12*r["p_heterogeneous_species_niche_segregation_upper"])
                r["p_signed_bonferroni_all_12_axis_depth"] = min(1.,12*r["p_signed_two_sided"])
    return results

def run(path,nrep=REPS):
    d=load_original(path);v=select_species(d);res=evaluate(v,nrep)
    return {"schema":"fcp_original100543_photo_heterogeneous_withinspecies_hue_climatic_niches_v1",
      "date_jst":"2026-10-11","status":"OBSERVATIONAL_RETROSPECTIVE",
      "source_original_archive_sha256":SHA,"n_full_unique_source_photo_ids":100543,
      "n_original_4state_classifiable_source_photos":45953,
      "n_original_classifiable_environment_photo_rows":len(d),
      "n_original_source_species_with_environment":int(d.inat_taxon_id.nunique()),
      "n_species_with_5_photo_2_each":len(v),
      "spatial_source_photo_label_exchangeability":source_spatial_exchange(v),
      "n_original_source_photo_label_permutations":nrep,
      "environment_axes":list(AXES),
      "comparison":"PHOTO white vs other colour environmental centroid within each species; squared contrasts detect opposite niche directions; species-conditional colour-shuffle vs species+10/20-degree-local colour shuffle",
      "results":res,"frozen_original_FCP_H1H2_unchanged":True,
      "hard_nonclaims":["Original photograph morphology classification is not independent botanist hue determination or genotype.",
        "Absolute white-vs-colour environmental niche contrast does not demonstrate causal pigment advantages.",
        "Spatial block-shuffle has low within-species photo exchangeability at 10 degrees (not a full spatial confounding removal).",
        "This is not a population-level white frequency-dependent fitness or payoff analysis.",
        "Original photo unclassifiable records remain unclassified and excluded only from colour response.",
        "All environmental predictors and three thresholds are tested/reported; no winner-only correction.",
        "The source original photo opportunity and reserved prospective taxon cohorts are untouched."]}

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source-archive",type=Path,required=True)
    p.add_argument("--outdir",type=Path,required=True)
    p.add_argument("--permutations",type=int,default=REPS)
    a=p.parse_args();z=run(a.source_archive,a.permutations)
    a.outdir.mkdir(parents=True,exist_ok=True)
    (a.outdir/"result.json").write_text(json.dumps(z,indent=2)+"\n")
    print(json.dumps({"source_rows":z["n_original_classifiable_environment_photo_rows"],"spatial_photo_exchangeability":z["spatial_source_photo_label_exchangeability"],
      "results":z["results"]},indent=2))
if __name__=="__main__":main()
