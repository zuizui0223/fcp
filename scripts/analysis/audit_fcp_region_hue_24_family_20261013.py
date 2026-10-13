#!/usr/bin/env python3
"""FCP original-source all-four-geographies x six-hue-pairs retrospective family.

Scientific context: selected eastern sector and white-versus-blue/purple hue pair
were chosen after historical flower-colour results were inspected. This analysis
makes their multiplicity explicit by evaluating all four original longitude
sectors x all six four-hue probability coefficient contrasts in matched genera.
Within each source/quality population a GLOBAL genus bootstrap is SHARED across
all four sector comparisons. This estimates a descriptive simultaneous envelope,
not confirmatory ecological p-values or fitness effects.
"""
from __future__ import annotations
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import numpy as np

from audit_fcp_same_genus_east_elsewhere_rain_tradeoff_20261011 import (
    get_source, eligible, sufficient_statistics, HUES, COLS, SHA,
)

SEED = 2026101308
REPLICATES = 999
MIN_NAMED_SPECIES_PER_GENUS_PER_SIDE = 3
MIN_GENERA = 30
SECTORS = ((-180.,-90.),(-90.,0.),(0.,90.),(90.,180.))
PAIRS = tuple(itertools.combinations(range(4),2))
CASE_NAMES = ("all_classified_original_source_photos", "original_photo_quality_strict")
SCHEMA = "fcp_four_sectors_six_hue_source_matched_genus_posthoc_v1"


def source_sectors(c):
    # The 90–180° E final sector is closed on its eastern endpoint, matching PR164.
    k = np.floor((c.longitude.to_numpy(float) + 180.)/90.).astype(int)
    return np.minimum(3,np.maximum(0,k))


def sector_genus_statistics(c, mu, sd, equal_genus:bool):
    """Gather unmodified parent model sufficient statistics for all 4 sectors."""
    result = {}
    sectors = source_sectors(c)
    if not np.all((sectors>=0)&(sectors<4)):
        raise RuntimeError("invalid original photo sectors")
    src=c.copy()
    src["fcp_sector"]=sectors
    for sector in range(4):
        src["fcp_inside"]=src.fcp_sector.eq(sector)
        counts=src.groupby(["genus","fcp_inside"]).inat_taxon_id.nunique().unstack(fill_value=0)
        if not {True,False}.issubset(set(counts.columns)):
            raise ValueError("missing original photographed region support")
        keep=counts.index[(counts[True]>=MIN_NAMED_SPECIES_PER_GENUS_PER_SIDE)&
                          (counts[False]>=MIN_NAMED_SPECIES_PER_GENUS_PER_SIDE)]
        sel=src.loc[src.genus.isin(keep)]
        meta={"sector":sector,"longitude_degrees":SECTORS[sector],
              "n_genera":len(keep),"n_photos":len(sel),
              "n_unique_named_species":int(sel.inat_taxon_id.nunique()),
              "n_region_named_species":int(sel[sel.fcp_inside].inat_taxon_id.nunique()),
              "n_other_named_species":int(sel[~sel.fcp_inside].inat_taxon_id.nunique())}
        if len(keep)<MIN_GENERA:
            result[sector]={"status":"HOLD_FEWER_THAN_30_SHARED_GENERA", "meta":meta}
            continue
        ins=sufficient_statistics(sel[sel.fcp_inside],mu,sd,
                                  MIN_NAMED_SPECIES_PER_GENUS_PER_SIDE,equal_genus)
        out=sufficient_statistics(sel[~sel.fcp_inside],mu,sd,
                                  MIN_NAMED_SPECIES_PER_GENUS_PER_SIDE,equal_genus)
        if set(ins)!=set(out) or set(ins)!=set(keep):
            raise RuntimeError("matched source genera changed between sides")
        result[sector]={"status":"OK","meta":meta,"inside":ins,"outside":out}
    return result


def simultaneous_source_bootstrap(regions, reps:int=REPLICATES,seed:int=SEED):
    """A global matched-genus multinomial redraw is shared by 4 sectors.

    This is NOT the parent marginal bootstrap: previously each region alone
    sampled only the qualifying genus pool. Here the universe is the UNION
    of eligible source genera across all four chosen sectors. Shared draws
    keep cross-sector, cross-colour covariance for the familywise maxT.
    """
    all_genera=sorted({g for reg in regions.values() if reg["status"]=="OK"
                        for g in reg["inside"]})
    G=len(all_genera)
    if G<MIN_GENERA:
        return {"status":"HOLD_FEWER_THAN_30_UNION_GENERA"}
    gen_index={g:i for i,g in enumerate(all_genera)}
    rng=np.random.default_rng(seed)
    weights=rng.multinomial(G,np.repeat(1/G,G),size=reps)
    points=np.full((4,len(PAIRS)),np.nan,dtype=float)
    draws=np.full((reps,4,len(PAIRS)),np.nan,dtype=float)
    metas=[]
    for sector in range(4):
        z=regions[sector]
        if z["status"]!="OK":
            metas.append({**z["meta"],"status":z["status"]})
            continue
        gen=sorted(z["inside"])
        ids=np.array([gen_index[g] for g in gen])
        w=weights[:,ids]
        estimates=[]
        simulated=[]
        for side in ("inside","outside"):
            xx=np.stack([z[side][g][0] for g in gen])
            xy=np.stack([z[side][g][1] for g in gen])
            orig_xx=xx.sum(axis=0)
            if np.linalg.cond(orig_xx)>1e10:
                raise RuntimeError(f"ill-conditioned original fitted design: {sector}/{side}")
            estimates.append(np.linalg.solve(orig_xx,xy.sum(axis=0))[0,:])
            bx=np.einsum("rg,gij->rij",w,xx,optimize=True)
            by=np.einsum("rg,gij->rij",w,xy,optimize=True)
            try:
                simulated.append(np.linalg.solve(bx,by)[:,0,:])
            except np.linalg.LinAlgError as e:
                raise RuntimeError(f"bootstrap singular: {sector}/{side}") from e
        point=estimates[0]-estimates[1]
        sim=simulated[0]-simulated[1]
        for k,(i,j) in enumerate(PAIRS):
            points[sector,k]=point[i]-point[j]
            draws[:,sector,k]=sim[:,i]-sim[:,j]
        metas.append({**z["meta"],"status":"OK",
                      "four_hue_inside_slopes_per_sd":dict(zip(HUES,estimates[0].tolist())),
                      "four_hue_outside_slopes_per_sd":dict(zip(HUES,estimates[1].tolist()))})
    included=np.isfinite(points)
    ntests=int(included.sum())
    if ntests!=24:
        return {"status":"HOLD_INCOMPLETE_24_DIRECTIONAL_TEST_FAMILY",
                "tested":ntests,"regions":metas}
    flat_draw=draws.reshape(reps,-1)
    flat_point=points.flatten()
    std=flat_draw.std(axis=0,ddof=1)
    if (std<=1e-12).any():
        return {"status":"HOLD_DEGENERATE_DIRECTIONAL_PAIR",
                "n_degenerate":int((std<=1e-12).sum())}
    max_t=np.max(np.abs((flat_draw-flat_point)/std),axis=1)
    critical=float(np.quantile(max_t,0.95))
    all_rows=[]
    for sector in range(4):
        for ip,(i,j) in enumerate(PAIRS):
            index=sector*len(PAIRS)+ip
            est=flat_point[index]
            ci=np.quantile(flat_draw[:,index],[0.025,0.975])
            sim_ci=(est+np.array([-1.,1.])*critical*std[index])
            all_rows.append({"sector":sector,"hue_A":HUES[i],"hue_B":HUES[j],
                             "point_pp_per_BIO12_SD":float(100*est),
                             "nominal_bootstrap_CI95_pp":(100*ci).tolist(),
                             "maxT_24_pair_simultaneous_approx_CI95_pp":(100*sim_ci).tolist(),
                             "maxT_excludes_zero":bool(sim_ci[0]>0 or sim_ci[1]<0),
                             "bootstrap_positive_fraction":float(np.mean(flat_draw[:,index]>0))})
    selected=[x for x in all_rows if x["sector"]==3 and x["hue_A"]=="white" and x["hue_B"]=="blue_purple"]
    if len(selected)!=1:
        raise RuntimeError("selected historical white-vs-blue pair missing")
    return {"status":"OK","n_union_genera":G,"bootstrap_replicates":reps,
            "shared_genus_draw_seed":seed,"n_test_family":ntests,
            "maxT_critical_95":critical,"sector_metadata":metas,
            "all_directional_contrasts":all_rows,
            "historical_white_blue_east":selected[0],
            "maxT_scope":"24=4 original longitude sectors x all 6 four-hue pairs within one source sample and taxonomic weight; posthoc photo/type/selection not controlled",
            "selection_adjustment_all_geography_hypotheses":False,
            "not_an_independent_validation":True}


def run(source_file:Path,repetitions:int=REPLICATES):
    if repetitions<199:
        raise ValueError("at least 199 paired bootstraps are required")
    if hashlib.sha256(source_file.read_bytes()).hexdigest()!=SHA:
        raise RuntimeError("FCP archived photo source SHA mismatch")
    rows=get_source(source_file)
    _,classifiable,_,mu,sd=eligible(rows)
    if len(classifiable)!=45836 or classifiable.inat_taxon_id.nunique()!=25028:
        raise RuntimeError("original 45,836 source photo climate support changed")
    qarea=float(classifiable.flower_effective_pixels.median())
    qnuis=float(classifiable.nuisance_pixel_fraction.median())
    variants={"all_classified_original_source_photos":classifiable,
              "original_photo_quality_strict":classifiable.loc[
                  (classifiable.flower_effective_pixels>=qarea)&
                  (classifiable.nuisance_pixel_fraction<=qnuis)].copy()}
    result={}
    for case,photo_frame in variants.items():
        result[case]={"original_photo_n":len(photo_frame),
                      "original_named_species_n":int(photo_frame.inat_taxon_id.nunique()),
                      "weighting":{}}
        for label,equal_genus in (("equal_species",False),("equal_genera",True)):
            st=sector_genus_statistics(photo_frame,mu,sd,equal_genus)
            obj=simultaneous_source_bootstrap(st,repetitions,
                seed=SEED+(1000 if equal_genus else 0)+(10000 if case=="original_photo_quality_strict" else 0))
            result[case]["weighting"][label]=obj
    return {"schema":SCHEMA,"role":"posthoc_evidence_and_falsification_not_confirmatory",
            "date":"2026-10-13","parent_pr":188,
            "old_selected_sector_and_colour_pair":"90–180 E; white minus blue/purple",
            "source_SHA256":SHA,"parent_source_colour_photos":45836,
            "parent_source_species":25028,"source_photo_classifiable_n":45953,
            "source_photo_opportunity_n":100543,"source_photo_unclassified_n":54590,
            "source_photo_quality_thresholds":{"minimum_flower_effective_pixels":qarea,
                                              "maximum_nuisance_pixel_fraction":qnuis},
            "cases":result,
            "claim_ceiling":["Four-sector/six-hue correction still conditions on retrospective hypothesis/design and the same old iNaturalist source.",
                             "No colour-altering rescue model, posthoc significance selection, biological adaptation, genetic pigment or fitness inference.",
                             "Photo quality changes taxon opportunity; source image exposure and unclassifiable colour MNAR remain.",
                             "Same genus and photo geography controls are not a dated phylogeny nor a complete photo-site spatial autocorrelation model."]}


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source-archive",type=Path,required=True)
    p.add_argument("--outdir",type=Path,required=True)
    p.add_argument("--bootstrap",type=int,default=REPLICATES)
    a=p.parse_args()
    obj=run(a.source_archive,repetitions=a.bootstrap)
    a.outdir.mkdir(parents=True,exist_ok=True)
    (a.outdir/"result.json").write_text(json.dumps(obj,indent=2,ensure_ascii=False)+"\n")
    print(json.dumps({
       "source_sha":obj["source_SHA256"],
       "family_summary":{k:{w:{"n_union_genera":v["n_union_genera"],
          "maxT_critical_95":v["maxT_critical_95"],
          "historical_white_blue_east":v["historical_white_blue_east"]}
          for w,v in z["weighting"].items() if v["status"]=="OK"}
          for k,z in obj["cases"].items()}},indent=2))

if __name__=="__main__": main()
