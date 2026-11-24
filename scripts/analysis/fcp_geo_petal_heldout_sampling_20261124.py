#!/usr/bin/env python3
"""FCP original-447500-photo heldout *design* for open-petal review.

Creates blinded near-vs-far matched *within-species* pair audit forms plus
classified/unclassified original photo review controls. The source photo labels
are used for audit *stratification only*, never as verified petal ground truth.
This does not score the images or refit the frozen spatial estimator.

Input must be the release-authenticated original measured-photo CSV.gz.
Run:
 python scripts/analysis/fcp_geo_petal_heldout_sampling_20261124.py \
   --source /tmp/source/measured/fcp_all4730_original_method_measured_photos.csv.gz \
   --out-dir /tmp/fcp-open-petal-geo-holdout
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[2]
EXPECTED_SHA="8e7f10f0427927208485b2da6fa9c46bb681d6f3611e53b8dfdc606cfc8dabe6"
EXPECTED_N=447500
EXPECTED_TAXA=4475
EXPECTED_CLASS=227107
SEED=20261124
PER_COMBINATION=24
PER_SINGLE_CLASS=64
PAIRS_PATH=ROOT/"manuscripts/climate_colour/Table_S11_pair_audit.csv"
PHOTO_FIG=ROOT/"manuscripts/climate_colour/Table_S1_photographs.csv"
LEDGER=ROOT/"results/fcp_np_all4475_unified_species_ledger_20261021.csv"
LABELS={"white","yellow_orange","red_pink","blue_purple"}
COMBOS=[(0,0),(0,1),(1,0),(1,1)] # near source mismatch / far source mismatch

def guard(ok,msg):
    if not bool(ok):
        raise AssertionError("FCP_HELDOUT_GEO_PETAL_SAMPLE: "+str(msg))

def sha(path):
    h=hashlib.sha256()
    with path.open("rb") as file:
        for b in iter(lambda:file.read(1024*1024),b""):
            h.update(b)
    return h.hexdigest()

def clean_source(path):
    guard(sha(path)==EXPECTED_SHA, "original 447500 source SHA256 differs")
    req=["inat_taxon_id","photo_id","observation_id","observer_id",
         "observed_on","latitude","longitude","positional_accuracy_m",
         "global_classifiable","morph"]
    d=pd.read_csv(path,usecols=req,low_memory=False)
    guard(len(d)==EXPECTED_N and d.photo_id.nunique()==EXPECTED_N,
          "expected 447500 distinct source photo IDs")
    guard(d.inat_taxon_id.nunique()==EXPECTED_TAXA,
          "source registry not 4475 taxa")
    guard(d.groupby("inat_taxon_id").size().eq(100).all(),
          "exactly 100 original photos per taxon expected")
    d["classified"]=d.global_classifiable.astype(str).str.lower().isin(("true","1","yes"))
    guard(int(d.classified.sum())==EXPECTED_CLASS,
          "original classifier denominator changed")
    guard(d.loc[d.classified,"morph"].isin(LABELS).all(),
          "unrecognized source classified photo label")
    older=pd.read_csv(PAIRS_PATH,
                      usecols=["left_photo_id","right_photo_id"])
    past=set(older.left_photo_id.astype(int))|set(older.right_photo_id.astype(int))
    figs=pd.read_csv(PHOTO_FIG,usecols=["photo_id"])
    past.update(figs.photo_id.astype(int))
    d["heldout_eligible_photo"]=~d.photo_id.isin(past)
    # Preserve the date as source-observed *calendar month*, with invalid
    # dates excluded rather than silently assigned to an arbitrary month.
    dt=pd.to_datetime(d.observed_on,utc=True,format="mixed",errors="coerce")
    d["month"]=dt.dt.month
    d["latitude"]=pd.to_numeric(d.latitude,errors="coerce")
    d["longitude"]=pd.to_numeric(d.longitude,errors="coerce")
    d["positional_accuracy_m"]=pd.to_numeric(d.positional_accuracy_m,errors="coerce")
    ledger=pd.read_csv(LEDGER,usecols=["inat_taxon_id","species","source_cohort"],
                       low_memory=False)
    guard(len(ledger)==EXPECTED_TAXA and ledger.inat_taxon_id.is_unique,
          "taxonomic ledger identity changed")
    d=d.merge(ledger,on="inat_taxon_id",validate="many_to_one")
    guard(len(d)==EXPECTED_N,"photo species lookup drift")
    return d,past

def haversine_pairs(g,i,j):
    lat=np.deg2rad(g.latitude.to_numpy(float))
    lon=np.deg2rad(g.longitude.to_numpy(float))
    h=np.sin((lat[i]-lat[j])/2)**2
    h+=np.cos(lat[i])*np.cos(lat[j])*np.sin((lon[i]-lon[j])/2)**2
    return 6371.0088*2*np.arcsin(np.sqrt(np.clip(h,0,1)))

def photo_pair_opportunities(d):
    rng=np.random.default_rng(SEED)
    by_combo={x:[] for x in COMBOS}
    ledger=[]
    for taxon,g0 in d.groupby("inat_taxon_id",sort=True):
        g=g0.loc[g0.heldout_eligible_photo & g0.classified &
                 g0.latitude.between(-90,90) & g0.longitude.between(-180,180) &
                 g0.positional_accuracy_m.between(0,5000) &
                 g0.month.notna() & g0.observer_id.notna()
                ].sort_values("photo_id").reset_index(drop=True)
        n=len(g)
        pairs={}
        counts={"near_same":0,"near_diff":0,"far_same":0,"far_diff":0}
        if n>=2:
            i,j=np.triu_indices(n,1)
            observer=g.observer_id.to_numpy()
            month=g.month.to_numpy()
            acceptable=(observer[i]!=observer[j])&(month[i]==month[j])
            if acceptable.any():
                distance=haversine_pairs(g,i,j)
                morph=g.morph.to_numpy()
                different=(morph[i]!=morph[j])
                for band,selection in (("near",distance<=50),
                                        ("far",distance>=250)):
                    for m in (0,1):
                        mask=acceptable & selection & (different==bool(m))
                        idx=np.flatnonzero(mask)
                        counts[f"{band}_{'diff' if m else 'same'}"]=len(idx)
                        if len(idx):
                            pick=int(rng.choice(idx))
                            left=g.iloc[int(i[pick])]
                            right=g.iloc[int(j[pick])]
                            pairs[(band,m)]={
                                "taxon_id":int(taxon),
                                "photo_A_id":int(left.photo_id),
                                "photo_B_id":int(right.photo_id),
                                "original_colour_A":str(left.morph),
                                "original_colour_B":str(right.morph),
                                "source_photo_month":int(left.month),
                                "distance_km":float(distance[pick]),
                                "source_observer_distinct":
                                    bool(left.observer_id!=right.observer_id),
                                "n_eligible_photo_pairs_in_species_stratum":
                                    int(len(idx)),
                                "n_classified_geovalid_photo":n
                            }
        for near,far in COMBOS:
            if ("near",near) in pairs and ("far",far) in pairs:
                by_combo[(near,far)].append({
                    "source_cohort":str(g0.source_cohort.iloc[0]),
                    "species":str(g0.species.iloc[0]),
                    "taxon_id":int(taxon),
                    "near":pairs[("near",near)],
                    "far":pairs[("far",far)]
                })
        ledger.append({
            "inat_taxon_id":int(taxon),
            "species":str(g0.species.iloc[0]),
            "source_cohort":str(g0.source_cohort.iloc[0]),
            "n_photos_geovalid_classified_for_audit":n,
            **counts,
            "n_available_joint_label_combinations":sum(
                ("near",near) in pairs and ("far",far) in pairs
                for near,far in COMBOS)
        })
    return by_combo,ledger

def select_pairs(by_combo,rand):
    chosen=[]
    counts={}
    for near,far in COMBOS:
        opportunities=by_combo[(near,far)]
        label=f"near_{'different' if near else 'same'}__far_{'different' if far else 'same'}"
        neligible=len(opportunities)
        counts[label]={"n_eligible_species":neligible,"n_selected_species":min(PER_COMBINATION,neligible)}
        guard(neligible>=PER_COMBINATION,
              f"underpowered heldout matched source stratum {label}: {neligible} species")
        for ix in rand.choice(neligible,size=PER_COMBINATION,replace=False):
            s=opportunities[int(ix)]
            for band in ("near","far"):
                r=dict(s[band])
                r.update({
                    "source_stratum":label,
                    "distance_group":band,
                    "joint_near_mismatch":near,
                    "joint_far_mismatch":far,
                    "source_cohort":s["source_cohort"],
                    "species":s["species"],
                    "n_eligible_source_species_in_stratum":neligible,
                    # Within stratum: species and then one compatible pair
                    # were chosen uniformly. Audit pairs' source labels
                    # identify the sampling stratum; no biological results yet.
                    "marginal_draw_inclusion_prob":float(
                        PER_COMBINATION/neligible/
                        r["n_eligible_photo_pairs_in_species_stratum"])
                })
                chosen.append(r)
    guard(len(chosen)==2*len(COMBOS)*PER_COMBINATION,
          "incomplete paired close/far audit")
    return chosen,counts

def single_image_opportunities(d,excluded_photos,rand):
    # This quality-control sample is separate from the near/far paired
    # source-contrast design. It intentionally includes source failures.
    outcomes=[];stats={}
    for classified in (False,True):
        g=d.loc[(d.classified==classified) &
                d.heldout_eligible_photo &
                ~d.photo_id.isin(excluded_photos)]
        bytaxon=[(int(t),x.sort_values("photo_id")) for t,x in
                 g.groupby("inat_taxon_id",sort=True)]
        label="classified" if classified else "unclassified"
        neligible=len(bytaxon)
        guard(neligible>=PER_SINGLE_CLASS, "too few "+label+" taxa")
        eligible_taxon_indices=rand.choice(neligible,size=PER_SINGLE_CLASS,replace=False)
        stats[label]={"n_eligible_taxa":neligible,"n_selected_taxa":PER_SINGLE_CLASS}
        for ix in eligible_taxon_indices:
            t,photos=bytaxon[int(ix)]
            row=photos.iloc[int(rand.integers(len(photos)))]
            outcomes.append({
                "taxon_id":t,
                "species":str(row.species),
                "source_cohort":str(row.source_cohort),
                "photo_id":int(row.photo_id),
                "original_classifiable":bool(classified),
                "original_coarse_photo_label":str(row.morph) if classified else None,
                "n_original_photos_in_species_stratum":int(len(photos)),
                "n_eligible_taxa_in_stratum":neligible,
                "marginal_draw_inclusion_prob":float(PER_SINGLE_CLASS/neligible/len(photos))
            })
    guard(len(outcomes)==2*PER_SINGLE_CLASS,"source classified/unclassified panel incomplete")
    return outcomes,stats

def blind_review_csv(out,rows,kind,rand):
    key=[]
    blind=[]
    for r in rows:
        token=(f"{SEED}:{kind}:{r.get('taxon_id')}:{r.get('source_stratum','')}:"+
               f"{r.get('distance_group','')}:{r.get('photo_A_id',r.get('photo_id'))}")
        ident=hashlib.sha256(token.encode()).hexdigest()[:12]
        r2=dict(r)
        r2["blind_id"]=ident
        if kind=="pairs":
            swap=bool(rand.integers(2))
            if swap:
                r2["photo_A_id"],r2["photo_B_id"]=r2["photo_B_id"],r2["photo_A_id"]
                r2["original_colour_A"],r2["original_colour_B"]=r2["original_colour_B"],r2["original_colour_A"]
            r2["source_A_is_original_order"]=not swap
            blind.append({
                "blind_id":ident,
                "photo_A_page":"https://www.inaturalist.org/photos/"+str(r2["photo_A_id"]),
                "photo_B_page":"https://www.inaturalist.org/photos/"+str(r2["photo_B_id"]),
                "photo_A_open_petal_colour":"",
                "photo_B_open_petal_colour":"",
                "both_photos_comparable_open_corollas":"",
                "same_developmental_stage":"",
                "non_petal_dominance":"",
                "photo_A_visible_open_corolla":"",
                "photo_B_visible_open_corolla":"",
                "photo_quality_sufficient":"",
                "reviewer_confidence":"",
                "reviewer_notes":""
            })
        else:
            blind.append({
                "blind_id":ident,
                "photo_page":"https://www.inaturalist.org/photos/"+str(r2["photo_id"]),
                "open_corolla_visible":"",
                "predominant_open_petal_colour":"",
                "other_organ_dominant":"",
                "background_or_nonflower_contribution":"",
                "flower_bud_or_fruit":"",
                "photo_quality_sufficient":"",
                "reviewer_confidence":"",
                "reviewer_notes":""
            })
        key.append(r2)
    guard(len({r["blind_id"] for r in blind})==len(blind),
          "blind source identifier collision")
    for who in ("A","B"):
        seq=list(blind)
        rand.shuffle(seq)
        f=out/f"reviewer_{who}_geographic_{kind}_blank.csv"
        pd.DataFrame([{"reviewer":who,**x} for x in seq]).to_csv(f,index=False)
    pd.DataFrame(key).to_csv(out/f"UNBLIND_AFTER_HUMAN_REVIEWS_{kind}_key.csv",index=False)
    return len(blind)

def main(source,out):
    d,past=clean_source(source)
    out.mkdir(parents=True,exist_ok=True)
    pool,ledger=photo_pair_opportunities(d)
    rand=np.random.default_rng(SEED+91)
    selected,counts=select_pairs(pool,rand)
    taken={r["photo_A_id"] for r in selected}|{r["photo_B_id"] for r in selected}
    single,stats=single_image_opportunities(d,taken,rand)
    selectedids=taken|{r["photo_id"] for r in single}
    guard(not (selectedids & past),"heldout audit reuses the earlier visual pilot or illustrated photos")
    guard(all(r["source_observer_distinct"] for r in selected),
          "same original observer in validation pair")
    guard(all((r["distance_km"]<=50 if r["distance_group"]=="near"
               else r["distance_km"]>=250) for r in selected),
          "paired distance band drift")
    pd.DataFrame(ledger).to_csv(out/"all4475_geography_open_petal_source_opportunity_ledger.csv",
                                index=False)
    npair=blind_review_csv(out,selected,"pairs",rand)
    nsingle=blind_review_csv(out,single,"singles",rand)
    receipt={
        "schema":"fcp447500_geographic_petal_validation_holdout_source_design_v1",
        "status":"BLANK_INDEPENDENT_HUMAN_REVIEW_QUEUES_NOT_PETAL_RESULTS",
        "source_sha256":EXPECTED_SHA,
        "n_original_photo_rows":len(d),
        "n_source_taxa":int(d.inat_taxon_id.nunique()),
        "n_source_classifiable":int(d.classified.sum()),
        "sampling_seed":SEED,
        "design":"Equal-probability species selection within near/far original photographic same/different 2x2 cells; one near and far source pair per sampled species-cell. Pair probabilities conditional on source label strata are recorded in private unblinding key. Both pairs have different original observers, same calendar month, source accuracy <=5km, near<=50km versus far>=250km. Independent classified/unclassified source photo control panel sampled one photo per species per original status.",
        "pair_strata":counts,
        "single_image_strata":stats,
        "selected_heldout_matched_taxon_draws":PER_COMBINATION*len(COMBOS),
        "selected_pairs":npair,
        "selected_single_source_photos":nsingle,
        "selected_unique_photo_ids":len(selectedids),
        "n_photo_ids_excluded_from_earlier_pilot_and_figures":len(past),
        "n_photo_id_overlap_with_prior_audit":len(selectedids & past),
        "minimum_photo_quality":"source original georeferenced classifiable pairs with observer id, same calendar month, and <=5000m coordinate uncertainty",
        "restricted_claims":[
            "No visual human scores are provided; these are BLANK review sheets only",
            "Stratification by original photo labels is deliberately outcome-enriched and NOT a random population estimate",
            "The dual-opportunity taxon frame and >=250km far comparator are different from the original all-range P2 estimand",
            "No claims of validated petals, genetic morph persistence, temperature adaptation or global classifier accuracy",
            "Pairs may share source photographs and species across strata; analyse and bootstrap by species and genus, not raw pair count",
            "Source photo pages provide links, not reproduced original images; reviewers must check licence, access and tissue",
            "Form scoring, actual third reader adjudication and source-weighted re-estimation are separate future steps"
        ]}
    (out/"fcp_geographic_petal_heldout_design_receipt.json").write_text(
        json.dumps(receipt,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    print(json.dumps(receipt,indent=2,ensure_ascii=False))

if __name__=="__main__":
    import argparse
    p=argparse.ArgumentParser()
    p.add_argument("--source",type=Path,required=True)
    p.add_argument("--out-dir",type=Path,required=True)
    a=p.parse_args()
    main(a.source,a.out_dir)
