#!/usr/bin/env python3
"""Original source photo candidate audit for field tests of white-vs-colour co-occurrence.

Do NOT infer genetically segregating floral morphs from this: only original
four-state photo labels/observer/time/coordinates. All geographically compatible
photograph sets are selected subject to deterministic predeclared gates.

No per-photograph colour labels, locations, original photo IDs, or old source
audit-case IDs are published. The species-level ledger is an internal candidate
screen; it is NOT a reviewer-facing 405-photo blind key.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd

SOURCE_SHA={
  "discovery":"ee854126eed2cfe23e333abe2c28d14df24895a5c52cbc389c060a5a4d6f91f4",
  "validation":"0e2ed349122739eecfc725fb2d5e313d284cf30752da91f9a0429cff0eeaa5e6",
  "third":"57630fc9f281bce94a0c40a70aaf7bce879dde93d6154175adcd021e8f5c1186",
}
HIGH_DEPTH={"discovery":369,"validation":363,"third":377}
STATES=("white","yellow_orange","red_pink","blue_purple")
RADII_KM=(10.,25.,50.)
EARTH_KM=6371.0088
MIN_CLASSIFIED=40
YEAR_MIN,YEAR_MAX=1990,2026

def file_sha(path:Path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda:f.read(1<<20),b""):h.update(chunk)
    return h.hexdigest()

def yes(x):
    if pd.api.types.is_bool_dtype(x):return x.fillna(False).astype(bool)
    return x.fillna("").astype(str).str.casefold().str.strip().isin(("true","yes","1"))

def load(path:Path,cohort:str):
    if file_sha(path)!=SOURCE_SHA[cohort]:raise RuntimeError(f"{cohort}: original photo source SHA mismatch")
    d=pd.read_csv(path,low_memory=False)
    req={"inat_taxon_id","species","photo_id","observer_id","observed_on","morph",
         "global_classifiable","latitude","longitude"}
    if not req.issubset(d):raise ValueError("Required original source photo columns missing")
    d=d.loc[yes(d.global_classifiable)&d.morph.isin(STATES),list(req)].copy()
    d["latitude"]=pd.to_numeric(d.latitude,errors="coerce")
    d["longitude"]=pd.to_numeric(d.longitude,errors="coerce")
    d=d.loc[d.latitude.between(-90,90)&d.longitude.between(-180,180)].copy()
    count=d.groupby("inat_taxon_id").size()
    d=d.loc[d.inat_taxon_id.isin(count[count>=MIN_CLASSIFIED].index)].copy()
    if d.inat_taxon_id.nunique()!=HIGH_DEPTH[cohort]:
        raise RuntimeError(f"{cohort}: source high-depth species eligibility has drifted")
    if d.photo_id.duplicated().any():raise RuntimeError("Original source photograph duplicate")
    d["observer"]=d.observer_id.fillna("").astype(str).str.strip()
    d.loc[d.observer.isin(("nan","None","<NA>")),"observer"]=""
    dt=pd.to_datetime(d.observed_on,errors="coerce")
    good=dt.notna()&dt.dt.year.between(YEAR_MIN,YEAR_MAX)
    d["year"]=dt.dt.year.where(good)
    d["month"]=dt.dt.month.where(good)
    d["has_valid_date"]=good
    d["is_white"]=d.morph.eq("white")
    return d.sort_values(["inat_taxon_id","photo_id"],kind="stable").reset_index(drop=True)

def geo_km(lat:np.ndarray,lon:np.ndarray):
    lat=np.deg2rad(np.asarray(lat,float))
    lon=np.deg2rad(np.asarray(lon,float))
    unit=np.stack([np.cos(lat)*np.cos(lon),np.cos(lat)*np.sin(lon),np.sin(lat)],axis=1)
    return np.arccos(np.clip(unit@unit.T,-1,1))*EARTH_KM

def find_mixed_pair(g:pd.DataFrame,dist:np.ndarray,max_diameter:float, mode:str):
    """Exact two-photo witness, distinguish observers, no colour/photo sorting."""
    n=len(g)
    group=np.arange(n)
    date=(g.year.to_numpy(float),g.month.to_numpy(float))
    white=g.is_white.to_numpy(bool)
    observer=g.observer.to_numpy(str)
    available=(observer!="")&g.has_valid_date.to_numpy(bool)
    u,v=np.triu_indices(n,k=1)
    gate=(dist[u,v]<=max_diameter)&(white[u]!=white[v])&available[u]&available[v]&(
        observer[u]!=observer[v])
    if mode=="same_month_year":
        gate&=(date[0][u]==date[0][v])&(date[1][u]==date[1][v])
    elif mode=="same_calendar_month":
        gate&=(date[1][u]==date[1][v])
    else:raise ValueError("Unknown temporal matching")
    idx=np.flatnonzero(gate)
    if not len(idx):return {"witness":False,"min_distance_km":None,"min_time_gap_years":None}
    order=np.lexsort((v[idx],u[idx],dist[u[idx],v[idx]]))
    a,b=u[idx[order[0]]],v[idx[order[0]]]
    return {"witness":True,"min_distance_km":float(dist[a,b]),
            "min_time_gap_years":int(abs(date[0][a]-date[0][b]))}

def four_distinct_observer_anchor(g:pd.DataFrame,dist:np.ndarray,max_diameter:float,mode:str):
    """Conservative single-photo-centred <=radius/2 witness.
    Guarantee ALL 4 source photos within full diameter; may miss off-centre
    complete-link 4-photo solutions. Never report NONE as impossibility.
    """
    white=g.is_white.to_numpy(bool)
    obs=g.observer.to_numpy(str)
    month=g.month.to_numpy(float)
    year=g.year.to_numpy(float)
    avail=g.has_valid_date.to_numpy(bool)&(obs!="")
    n=len(g)
    best=None
    # Best witness uses geometric/temporal keys only, not a source colour-leverage score.
    for anchor in range(n):
        if not avail[anchor]:continue
        near=(dist[anchor]<=max_diameter/2)&avail&(month==month[anchor])
        if mode=="same_month_year":near&=(year==year[anchor])
        elif mode!="same_calendar_month":raise ValueError(mode)
        ix=np.flatnonzero(near)
        if len(ix)<4:continue
        w=[int(i) for i in ix if white[i]]
        c=[int(i) for i in ix if not white[i]]
        if len(w)<2 or len(c)<2:continue
        # Source original same-observer multiple photos do not count as independent.
        for i in w:
            for j in w:
                if j<=i or obs[i]==obs[j]:continue
                for k in c:
                    if obs[k] in (obs[i],obs[j]):continue
                    for l in c:
                        if l<=k or obs[l] in (obs[i],obs[j],obs[k]):continue
                        yspan=float(max(year[[i,j,k,l]])-min(year[[i,j,k,l]]))
                        maxsep=float(dist[np.ix_([i,j,k,l],[i,j,k,l])].max())
                        if maxsep>max_diameter+1e-6:
                            raise RuntimeError("four-photo conservative geographic diameter violated")
                        cand=(yspan,maxsep,anchor,i,j,k,l)
                        if best is None or cand<best:best=cand
    if best is None:return {"witness":False,"maximum_distance_km":None,"span_years":None}
    return {"witness":True,"maximum_distance_km":float(best[1]),
            "span_years":float(best[0])}

def species_screen(g:pd.DataFrame,cohort:str):
    if len(g)<MIN_CLASSIFIED:raise RuntimeError("source classifiability depth")
    nwhite=int(g.is_white.sum()); ncolour=int((~g.is_white).sum())
    row={"cohort":cohort,"inat_taxon_id":int(g.inat_taxon_id.iloc[0]),
         "species":str(g.species.iloc[0]),"genus":str(g.species.iloc[0]).split(" ")[0],
         "n_original_classified_photo_records":int(len(g)),"n_white_photo_labels":nwhite,
         "n_nonwhite_photo_labels":ncolour,"photo_white_fraction":nwhite/len(g),
         "n_distinct_source_observers_with_id":int(g.loc[g.observer.ne(""),"observer"].nunique()),
         "n_source_photos_with_valid_date":int(g.has_valid_date.sum()),
         "n_date_years":int(g.loc[g.has_valid_date,"year"].nunique()),
         "is_mixed_photo_species":bool(nwhite>=1 and ncolour>=1),
         "has_five_each_photo_colours":bool(nwhite>=5 and ncolour>=5)}
    coord=geo_km(g.latitude.to_numpy(float),g.longitude.to_numpy(float))
    for rad in RADII_KM:
        for mode in ("same_month_year","same_calendar_month"):
            prefix=f"{int(rad)}km__{mode}"
            pair=find_mixed_pair(g,coord,rad,mode)
            quad=four_distinct_observer_anchor(g,coord,rad,mode)
            row[f"{prefix}__mixed_different_observer_photo_pair"]=pair["witness"]
            row[f"{prefix}__mixed_photo_pair_min_distance_km"]=pair["min_distance_km"]
            row[f"{prefix}__mixed_photo_pair_year_gap"]=pair["min_time_gap_years"]
            row[f"{prefix}__four_photo_two_each_four_observers_anchor"]=quad["witness"]
            row[f"{prefix}__four_photo_group_max_diameter_km"]=quad["maximum_distance_km"]
            row[f"{prefix}__four_photo_group_year_span"]=quad["span_years"]
    return row

def run(paths):
    cohort_out={}
    frames=[]
    for cohort,path in paths.items():
        photo=load(path,cohort)
        rows=[]
        for _,g in photo.groupby("inat_taxon_id",sort=True):
            rows.append(species_screen(g,cohort))
        data=pd.DataFrame(rows)
        frames.append(data)
        cohort_out[cohort]={"n_source_high_depth_species":len(data),
                 "n_mixed_original_classified_white_nonwhite":int(data.is_mixed_photo_species.sum()),
                 "n_at_least_five_white_and_five_nonwhite":int(data.has_five_each_photo_colours.sum()),
                 "n_original_photo_rows":len(photo)}
    df=pd.concat(frames,ignore_index=True)
    summary={}
    for rad in RADII_KM:
        summary[f"{int(rad)}km"]={}
        for mode in ("same_month_year","same_calendar_month"):
            prefix=f"{int(rad)}km__{mode}"
            pair=df[f"{prefix}__mixed_different_observer_photo_pair"]
            quad=df[f"{prefix}__four_photo_two_each_four_observers_anchor"]
            summary[f"{int(rad)}km"][mode]={
                "n_species_distinct_observer_mixed_source_photo_pair":int(pair.sum()),
                "n_species_four_photo_two_each_four_observers_conservative_anchor":int(quad.sum()),
                "n_genera_four_photo":int(df.loc[quad,"genus"].nunique()),
                "per_cohort":{c:{"n_pair":int(pair[df.cohort.eq(c)].sum()),
                                 "n_four":int(quad[df.cohort.eq(c)].sum())}
                              for c in paths}}
    out={
        "schema":"fcp_source_photo_mixed_white_nonwhite_site_temporal_field_candidates_v1",
        "date_jst":"2026-10-10",
        "role":"Retrospective true-field-candidate feasibility; NOT expert-verified sympatric biological morphs",
        "n_source_species":len(df),
        "cohorts":cohort_out,
        "photo_site_diameters_km":list(RADII_KM),
        "time_windows":["same_month_year","same_calendar_month_across_years"],
        "source_photo_pair_requirement":"one automated-photo-white and one automated-photo-nonwhite, two different nonmissing observers",
        "conservative_four_photo_requirement":"two automatically photo-white and two photo-nonwhite, four different original photographer IDs, all photos within one anchor of <=diameter/2",
        "all_photo_witnesses_selected_from_original_labels_and_coordinates_only":True,
        "no_actual_biological_morph_confirmed":True,
        "no_actual_fitness_or_genetic_data":True,
        "source_sha256":SOURCE_SHA,
        "summary":summary,
        "hard_nonclaims":[
           "A same-month-year photo witness does not prove two colour morphs coexist as independent reproductive genets.",
           "Original photographed white/nonwhite machine classes are NOT chemically/genetically authenticated.",
           "A conservative anchored four-photo witness misses some legitimate full-diameter groups: negative means no witness found with this method, not physical impossibility.",
           "Photos may be of different species or different floral organs; botanical 405 photo review is separate.",
           "Photo records clustered in a 10km region are not guaranteed one freely interbreeding natural population.",
           "Same calendar month across different years does NOT establish same-generation sympatric polymorphism.",
           "Ranking among source machine-classified pictures is not an unbiased estimate of polymorphism prevalence.",
           "Main 1499-source New Phytologist H1/H2, 405 human review status and untouched future source taxa are unchanged.",
        ]}
    return out,df

def main():
    p=argparse.ArgumentParser()
    for c in SOURCE_SHA:p.add_argument("--"+c,type=Path,required=True)
    p.add_argument("--outdir",type=Path,required=True)
    args=p.parse_args()
    output,df=run({k:getattr(args,k) for k in SOURCE_SHA})
    out=args.outdir;out.mkdir(parents=True,exist_ok=True)
    # INTERNAL only: species aggregates but original algorithm labels by species;
    # keep out of any reviewer-facing packet.
    df.to_csv(out/"INTERNAL_source_species_field_candidate_ledger_do_not_show_405_reviewers.csv",index=False)
    (out/"result.json").write_text(json.dumps(output,indent=2)+"\n")
    print(json.dumps(output,indent=2))
if __name__=="__main__":main()
