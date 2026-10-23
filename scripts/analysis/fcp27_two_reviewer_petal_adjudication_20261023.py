#!/usr/bin/env python3
"""Two-independent-reviewer FCP petal colour adjudication, with a strict null gate.

Never fabricate human labels. On blank review sheets, return a PENDING manifest
and no biological colour estimate. After two independently completed sheets,
freeze their SHA256 and the source mapping, evaluate petal-stage validity
before comparing low vs high original photo-PC1 strata, and keep morphology
distinct from inherited/genetic flower-colour polymorphism.
"""
from __future__ import annotations
import argparse, hashlib, json
from pathlib import Path
import numpy as np
import pandas as pd

N_IMAGES = 594
N_SPECIES = 54
REQUIRED = [
    "review_id", "target_flower_visible", "petal_region_unambiguous",
    "observed_dominant_colour","flower_stage","lighting_white_balance",
    "multiple_flower_or_organ_confusion","cultivation_or_taxon_concern",
    "review_confidence"
]
EXTRA = ["target_organ", "petal_pattern"]
STATUS = {"yes","no","uncertain","not_available"}
COLORS = {"white","cream","yellow","orange","red","pink","magenta",
          "purple","blue","green","mixed","other","unassessable"}
STAGE = {"open_fresh","bud","senescent","fruit_postfloral","mixed",
         "unknown","not_available"}
ORGAN = {"corolla_petal","petaloid_sepal","sepal_bract","stamen",
         "fruit","bud","vegetative","mixed","uncertain","not_available"}
PATTERN = {"uniform","basal_blotch","veined","bicolour_sector",
           "irregular","unknown","not_available"}


def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for part in iter(lambda:f.read(1<<20),b""):
            h.update(part)
    return h.hexdigest()


def req(v,why):
    if not v:raise RuntimeError("FCP27_PETAL_TWO_REVIEWER: "+why)


def prepare_blank(frame):
    d=frame.copy()
    for col in EXTRA:
        if col not in d.columns:d[col]=""
    return d


def read_reviewer(path,tag):
    d=prepare_blank(pd.read_csv(path,dtype=str,keep_default_na=False))
    req(len(d)==N_IMAGES and d.review_id.is_unique,
        tag+" must have exactly 594 unique source blind ids")
    req(set(REQUIRED).issubset(d.columns),tag+" missing mandatory column")
    return d


def completed(row):
    cols=["target_flower_visible","petal_region_unambiguous",
          "observed_dominant_colour","flower_stage",
          "lighting_white_balance","multiple_flower_or_organ_confusion",
          "cultivation_or_taxon_concern","review_confidence",
          "target_organ","petal_pattern"]
    return all(str(row[c]).strip()!="" for c in cols)


def score_one(row):
    return (row["target_flower_visible"]=="yes"
            and row["petal_region_unambiguous"]=="yes"
            and row["flower_stage"]=="open_fresh"
            and row["target_organ"] in {"corolla_petal","petaloid_sepal"}
            and row["lighting_white_balance"] in {"no","acceptable"}
            and row["multiple_flower_or_organ_confusion"]=="no"
            and row["review_confidence"] in {"high","medium"}
            and row["observed_dominant_colour"] not in {"unassessable","other","mixed"}
            and row["cultivation_or_taxon_concern"] in {"no","uncertain"})


def most_common(x):
    vc=x.value_counts()
    if vc.empty:return None,0
    return str(vc.index[0]),int(vc.iloc[0])


def agreement(a,b):
    if len(a)==0:return None
    return float(np.mean(a.to_numpy()==b.to_numpy()))


def calc(r1,r2,mapping,sha_record):
    r1=prepare_blank(r1)
    r2=prepare_blank(r2)
    req(len(r1)==N_IMAGES and len(r2)==N_IMAGES,
        "incomplete matched review universe")
    req(r1.review_id.is_unique and r2.review_id.is_unique,
        "duplicate reviewer identifiers")
    req(set(r1.review_id)==set(r2.review_id),
        "source image identities differ between independent reviewers")
    req(len(mapping)==N_IMAGES and mapping.review_id.is_unique,
        "sealed PI source identity mapping incomplete")
    req(mapping.inat_taxon_id.nunique()==N_SPECIES,
        "source review has not exactly 54 taxa")
    pending={
        "schema":"fcp27_two_independent_botanical_review_adjudication_v1",
        "source_n_selected_images":N_IMAGES,
        "source_n_taxa":N_SPECIES,
        "reviewer1_n_fully_labelled":int(sum(completed(v) for _,v in r1.iterrows())),
        "reviewer2_n_fully_labelled":int(sum(completed(v) for _,v in r2.iterrows())),
        "input_sha256":sha_record,
        "no_photo_colour_or_genetic_fcp_imputation":True,
        "status":None,
        "estimates":None,
    }
    if (pending["reviewer1_n_fully_labelled"]!=N_IMAGES
        or pending["reviewer2_n_fully_labelled"]!=N_IMAGES):
        pending["status"]="human_review_incomplete_not_estimable"
        return pending,pd.DataFrame()
    # Only after all human review rows have been annotated do identities get joined.
    o=r1.merge(r2,on="review_id",suffixes=("_r1","_r2"),
               validate="one_to_one").merge(
        mapping[["review_id","species","inat_taxon_id","review_arm",
                 "source_stratum","pair_no"]],
        on="review_id",validate="one_to_one")
    req(len(o)==N_IMAGES,"human review records changed after unblinding")
    for col,choices in (("target_flower_visible",STATUS),
                        ("petal_region_unambiguous",STATUS),
                        ("observed_dominant_colour",COLORS),
                        ("flower_stage",STAGE),("target_organ",ORGAN),
                        ("petal_pattern",PATTERN)):
        for side in ("r1","r2"):
            req(o[f"{col}_{side}"].isin(choices).all(),
                f"invalid {col} values for {side}")
    a=o[[f"target_flower_visible_{s}" for s in ("r1","r2")]]
    val1=pd.Series([score_one(v) for _,v in r1.iterrows()],index=r1.review_id)
    val2=pd.Series([score_one(v) for _,v in r2.iterrows()],index=r2.review_id)
    o["review1_valid"]=o.review_id.map(val1).astype(bool)
    o["review2_valid"]=o.review_id.map(val2).astype(bool)
    o["both_valid"]=(o.review1_valid&o.review2_valid)
    o["agreed_colour"]=(o.both_valid
                       &(o.observed_dominant_colour_r1==o.observed_dominant_colour_r2))
    o["agreed_pattern"]=(o.both_valid
                        &(o.petal_pattern_r1==o.petal_pattern_r2))
    species=[]
    for (tid,arm),g in o.groupby(["inat_taxon_id","review_arm"],sort=True):
        low=g.loc[g.source_stratum=="low_pc1"]
        high=g.loc[g.source_stratum=="high_pc1"]
        color={}
        patt={}
        eligible=True
        for name,subset in (("low_pc1",low),("high_pc1",high)):
            yes=subset.loc[subset.agreed_colour]
            hue,count=most_common(yes.observed_dominant_colour_r1)
            pat,npat=most_common(subset.loc[subset.agreed_pattern].petal_pattern_r1)
            color[name]={"n_consensus_valid":len(yes),
                         "dominant_hue":hue,"dominant_n":count}
            patt[name]={"dominant_pattern":pat,"dominant_n":npat}
            if len(yes)<2 or count<2:eligible=False
        colour_contrast=(eligible
                         and color["low_pc1"]["dominant_hue"]!=color["high_pc1"]["dominant_hue"])
        pattern_contrast=(eligible and patt["low_pc1"]["dominant_n"]>=2
                          and patt["high_pc1"]["dominant_n"]>=2
                          and patt["low_pc1"]["dominant_pattern"]!=patt["high_pc1"]["dominant_pattern"])
        species.append({
            "inat_taxon_id":int(tid),"species":str(g.species.iloc[0]),
            "review_arm":str(arm),"pair_no":int(g.pair_no.iloc[0]),
            "n_total_photographs":len(g),
            "n_both_reviewers_agree_valid":int(g.both_valid.sum()),
            "n_consensus_petal_colour":int(g.agreed_colour.sum()),
            "n_consensus_basal_pattern":int(g.agreed_pattern.sum()),
            "low_consensus":color["low_pc1"],"high_consensus":color["high_pc1"],
            "photo_corolla_hue_contrast_evaluable":bool(eligible),
            "photo_corolla_hue_contrast_supported_if_evaluable":bool(colour_contrast) if eligible else None,
            "petal_pattern_contrast_if_evaluable":bool(pattern_contrast) if eligible else None,
            "not_inherited_genetic_morphs":True
        })
    out=pd.DataFrame(species)
    req(len(out)==N_SPECIES,"not one row per taxon after review")
    tally={}
    for arm,g in out.groupby("review_arm"):
        known=g.loc[g.photo_corolla_hue_contrast_evaluable]
        tally[arm]={
            "n_taxa":len(g),
            "n_evaluable_after_both_reviewers":len(known),
            "n_consensus_photo_petals_hue_contrast":int(
                known.photo_corolla_hue_contrast_supported_if_evaluable.sum()),
            "n_unknown_not_equivalent_to_monochrome":len(g)-len(known)
        }
    pending.update({
        "status":"complete_two_reviewers_photo_based_not_genetic",
        "estimates":{
            "photo_validity_agreement":agreement(o.review1_valid,o.review2_valid),
            "raw_observed_hue_agreement":agreement(
                o.observed_dominant_colour_r1,o.observed_dominant_colour_r2),
            "photo_hue_evaluable_by_arm":tally,
            "requires_taxon_level_adjudication_for_divergent_reviewer_labels":True,
            "natural_population_repeated_genotype_or_fitness_evidence":False,
        }
    })
    return pending,out


def synthetic():
    a=pd.DataFrame([{"review_id":f"x{i}",
      "target_flower_visible":"","petal_region_unambiguous":"",
      "observed_dominant_colour":"","flower_stage":"",
      "lighting_white_balance":"","multiple_flower_or_organ_confusion":"",
      "cultivation_or_taxon_concern":"","review_confidence":"",
      "target_organ":"","petal_pattern":""} for i in range(N_IMAGES)])
    mapping=pd.DataFrame({"review_id":a.review_id,
             "inat_taxon_id":np.repeat(np.arange(N_SPECIES),11),
             "species":[f"Taxon {j}" for j in np.repeat(np.arange(N_SPECIES),11)],
             "review_arm":np.tile(np.repeat(["candidate","matched_control"],11),27),
             "source_stratum":np.tile(["low_pc1"]*3+["central_pc1"]*3+
                                      ["high_pc1"]*3+["unclassified_photo"]*2,N_SPECIES),
             "pair_no":np.repeat(np.arange(N_SPECIES)%27+1,11)})
    status,rows=calc(a,a,mapping,{})
    req(status["status"]=="human_review_incomplete_not_estimable"
        and status["estimates"] is None and rows.empty,
        "unreviewed photographs got invented human result")
    print("FCP27_TWO_REVIEWER_PENDING_GUARD_SYNTHETIC_PASS")


def main():
    p=argparse.ArgumentParser()
    for col in ("reviewer1","reviewer2","mapping","out"):
        p.add_argument("--"+col,type=Path)
    p.add_argument("--self-test",action="store_true")
    a=p.parse_args()
    synthetic()
    if a.self_test:return
    req(a.reviewer1 and a.reviewer2 and a.mapping and a.out,
        "both independent reviews and sealed PI map must be supplied")
    r1=read_reviewer(a.reviewer1,"reviewer1")
    r2=read_reviewer(a.reviewer2,"reviewer2")
    mapping=pd.read_csv(a.mapping,dtype={"review_id":str})
    receipt,per_taxon=calc(r1,r2,mapping,{
        "reviewer1":sha(a.reviewer1),"reviewer2":sha(a.reviewer2),
        "original_source_mapping":sha(a.mapping)
    })
    a.out.mkdir(parents=True,exist_ok=True)
    (a.out/"human_review_decision.json").write_text(
        json.dumps(receipt,indent=2,allow_nan=False)+"\n")
    if not per_taxon.empty:
        per_taxon.to_json(a.out/"per_taxon_human_photo_outcomes.json",
                          orient="records",indent=2)
    print("FCP27_TWO_REVIEWER_ADJUDICATION",
          json.dumps({k:v for k,v in receipt.items() if k!="input_sha256"}))


if __name__=="__main__":
    main()
