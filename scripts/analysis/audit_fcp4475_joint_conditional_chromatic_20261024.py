#!/usr/bin/env python3
"""Independent chromatic-only maxT audit on SHA-locked FCP joint climate results.

Reads original-photo seven-climate source NPZ arrays, NOT reconstructed
coefficients from a manuscript description. Confidence intervals are
whole-genus 999 redraws with maxT over the actual number of eligible
climate-by-eight-chromatic axes for EACH source panel and PHOTO component.
Do not call post-hoc chromatic hue choices preregistered causal effects.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

SOURCE_SHA={
 "summary":"cd29072f6d9f31966cc8fc7696866822000f1987d753704b7d34fad059306174",
 "original__thermal_moisture_solar_core_npz":"6fb31a6ea0369f773af2f9ecf52341d28a4cffa8887a7a269229d13a75217baf",
 "original__seven_climates_joint_conditional_npz":"971434a797933cfadfd035a9337c85e1c6ddc139e6b365aa5c223404c187dd44",
 "observer_thinned__thermal_moisture_solar_core_npz":"fb6d64447191a8d1fefdb6035d44da6a3bc1253f1cd1ec206078defeaf8f89db",
 "observer_thinned__seven_climates_joint_conditional_npz":"ffa2de03dc70e2b65147cb5d3f9efbda100fb46b741ee3ca13f473c80ba29aee",
 "original__thermal_moisture_solar_core_csv":"8f8bd18e6fde4c09730d4f79e1e13676aa6e9c0130772d489a0c73e1ef42bb68",
 "original__seven_climates_joint_conditional_csv":"9e47f13472493d4583fea751991a2e4aaae3dacb711cf99ab926f607291ebcb7",
 "observer_thinned__thermal_moisture_solar_core_csv":"9aee32b6a7b6abeeef8e58241b5fc3085880aa736ea5cb5bfe9a434db6abdb89",
 "observer_thinned__seven_climates_joint_conditional_csv":"a6acdc73d3bfbef643c88b4695e2c85d98dd098634fbe0db326c94e22146f28b",
}
MODELS=("thermal_moisture_solar_core","seven_climates_joint_conditional")
PANELS=("original","observer_thinned")
CLIM=("BIO1","BIO4","BIO5","BIO12","BIO14","BIO15","SRADmean")
HUE=("white","yellow","orange","red","pink","magenta","purple","blue","bronze")
COMP=("total","between_original_photo_labels","within_original_photo_labels")
Q=np.array([8.]+[-1.]*8)/np.sqrt(72.)
REPS=999
SEED=20261024


def check(c,reason):
    if not c: raise RuntimeError("FCP4475_CONDITIONAL_CHROMATIC_AUDIT: "+reason)


def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for part in iter(lambda:f.read(1<<20),b""):h.update(part)
    return h.hexdigest()


def mean_boot(taxa,coeff,reps=REPS,seed=SEED):
    group,genus=pd.factorize(taxa["genus"],sort=True)
    G=len(genus)
    valid=np.isfinite(coeff)
    S=np.zeros((G,7,3,9),float)
    N=np.zeros_like(S)
    np.add.at(S,group,np.where(valid,coeff,0.))
    np.add.at(N,group,valid.astype(float))
    sm=S.sum(axis=0)
    nn=N.sum(axis=0)
    mean=np.divide(sm,nn,out=np.full_like(sm,np.nan),where=nn>0)
    draws=np.random.default_rng(seed).integers(G,size=(reps,G))
    weights=np.zeros((reps,G),float)
    np.add.at(weights,(np.repeat(np.arange(reps),G),draws.ravel()),1)
    BS=weights@S.reshape(G,-1)
    BN=weights@N.reshape(G,-1)
    boot=np.divide(BS,BN,out=np.full_like(BS,np.nan),where=BN>0)
    return mean,nn,boot.reshape((reps,7,3,9)),G


def analyse_component(mean,denom,boot,comp):
    # The null-removed white colour contrast is orthogonalized exactly.
    point=mean[:,comp,:]
    samples=boot[:,:,comp,:]
    achro=np.einsum("ah,h->a",point,Q)
    b_achro=np.einsum("bah,h->ba",samples,Q)
    chrom=point-achro[:,None]*Q
    b_chrom=samples-b_achro[:,:,None]*Q
    check(np.allclose(chrom[np.isfinite(chrom).all(axis=1),0],0.,atol=1e-10),
          "chromatic residual inadvertently contains white source photo coordinate")
    test_axes=np.flatnonzero(np.isfinite(chrom).all(axis=1) & (denom[:,comp,0]>0))
    M=int(len(test_axes)*8)
    if not M:
        return {
            "status":"not_estimable","chromatic_evaluable_8hue_tests":0,
            "per_axis":{}
        }
    v=chrom[test_axes,1:].reshape(M)
    draws=b_chrom[:,test_axes,1:].reshape((REPS,M))
    check(np.isfinite(v).all() and np.isfinite(draws).all(),
          "bootstrap lost original identifiable chromatic source taxon")
    sd=draws.std(axis=0,ddof=1)
    check(np.all(sd>1e-13),"genus chromatic bootstrap zero variance")
    crit=float(np.quantile(np.max(np.abs((draws-v)/sd),axis=1),.95))
    lo=v-crit*sd
    hi=v+crit*sd
    per_axis={}
    for ai,axisid in enumerate(test_axes):
        axis=CLIM[axisid]
        results={}
        for j,hue in enumerate(HUE[1:]):
            k=ai*8+j
            results[hue]={
                "conditional_chromatic_photo_palette_mean":float(v[k]),
                "simultaneous95_chromatic_climate_family":[float(lo[k]),float(hi[k])],
                "genus_bootstrap_marginal95":[
                    float(x) for x in np.quantile(draws[:,k],[.025,.975])],
                "family_excludes_zero":bool(lo[k]>0 or hi[k]<0),
            }
        per_axis[axis]={
            "n_original_source_taxa":int(denom[axisid,comp,0]),
            "n_evaluable_chromatic_components":8,
            "chromatic_photo_hue_components":results,
            "n_simultaneous_family_effects":M,
        }
    return {
        "status":"genus_bootstrap_chromatic_8hue_conditional_photo_source",
        "chromatic_evaluable_8hue_tests":M,
        "maxT_genus95_critical":crit,
        "per_axis":per_axis,
    }


def audit(root):
    sp=root/"joint_climate_conditional_hue_summary.json"
    check(sha(sp)==SOURCE_SHA["summary"],"source result JSON SHA changed")
    canonical=json.loads(sp.read_text())
    check(canonical["original_photo_count"]==447500 and
          canonical["original_taxon_count"]==4475 and
          canonical["original_photo_classified_count"]==227107,
          "source original 447,500 photo/4,475 taxon denominators changed")
    results={}
    for idx,panel in enumerate(PANELS):
        for mdl in MODELS:
            key=panel+"__"+mdl
            p=root/f"original4475_{key}_nine_hue_conditional_source.npz"
            c=root/f"original4475_{key}_conditional_climate_support.csv"
            check(sha(p)==SOURCE_SHA[key+"_npz"] and sha(c)==SOURCE_SHA[key+"_csv"],
                  "authenticated original photo taxa or tensor bytes changed "+key)
            frame=pd.read_csv(c,low_memory=False)
            check(len(frame)==4475 and frame.inat_taxon_id.is_unique and
                  int(frame.n_original_classifiable.sum())==227107,
                  "source 4475 species photo opportunity lost")
            with np.load(p,allow_pickle=False) as z:
                coef=z["conditional_photo_palette_slope"]
                ids=z["inat_taxon_id"]
                check(np.array_equal(ids,frame.inat_taxon_id.to_numpy()),
                      "taxonomic source order differs from climate tensor")
                check(coef.shape==(4475,7,3,9),"joint-climate photo colour tensor corrupted")
                check(list(z["climate_axes"])==list(CLIM) and list(z["photo_hues"])==list(HUE),
                      "original seven climate/nine hue source order drift")
            mean,n,bootstrap,G=mean_boot(frame,coef,seed=SEED+idx*191+len(mdl)*13)
            full={}
            for j,comp in enumerate(COMP):
                full[comp]=analyse_component(mean,n,bootstrap,j)
            results[key]={
                "n_genus_source_clusters":int(G),
                "source_n_taxa":4475,
                "n_with_any_joint_climate":int(frame.status.eq("source_joint_climate_identifiable").sum()),
                "chromatic_by_original_photo_component":full
            }
    return {
        "schema":"fcp4475_joint_climate_source_nine_hue_chromatic_independent_audit_v1",
        "source_photo_run_id":40061576683,
        "source_photo_artifact_id":12547264579,
        "source_photo_artifact_zip_sha256":"ce145c38144d833be48b1a00a21414f6fb9e6d234e7101c5183e0f5d9e5dede9",
        "source_array_sha256":SOURCE_SHA,
        "n_original_source_photos":447500,
        "n_original_taxa":4475,
        "n_source_classifiable_photos":227107,
        "bootstrap_whole_genera":REPS,
        "multiple_testing":"simultaneous maxT within each of 4 source panels x 3 original photographic label component families across eligible climate axes x 8 chromatic hues, not jointly corrected over these 12 exploratory families",
        "white_axis_projection_removed_exactly":True,
        "original_photo_colour_categories_not_heritable_morphs":True,
        "chromatic_hue_posthoc_not_preregistered":True,
        "no_background_corrected_independently_measured_petals":True,
        "climate_conditional_coefficients_not_causal":True,
        "results":results,
    }


def selftest():
    example=np.array([.2,-.05,-.02,-.03,-.04,-.01,-.02,-.02,-.01])
    example-=example.mean()
    check(np.isfinite(example).all(),"synthetic source hue not finite")
    projection=np.dot(example,Q)
    residual=example-projection*Q
    check(abs(residual[0])<1e-12,
          "a purely chromatic source photo contrast retains white")
    print("FCP4475_JOINT_CHROMATIC_INDEPENDENT_SELFTEST_PASS")


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source",type=Path)
    p.add_argument("--out",type=Path)
    p.add_argument("--self-test",action="store_true")
    a=p.parse_args()
    selftest()
    if a.self_test:return
    check(a.source and a.out,"original source archive output and audit folder required")
    output=audit(a.source)
    a.out.mkdir(parents=True,exist_ok=True)
    (a.out/"joint_conditional_chromatic_56hue_genus_audit.json").write_text(
        json.dumps(output,indent=2,sort_keys=True,allow_nan=False)+"\n")
    print("FCP4475_JOINT_CHROMATIC_INDEPENDENT_AUDIT_COMPLETE",
          json.dumps({name:{
            component:{
                "n_effects":v["chromatic_evaluable_8hue_tests"],
                "BIO12":v["per_axis"].get("BIO12",{}),
                "BIO5":v["per_axis"].get("BIO5",{}),
            }
            for component,v in blocks["chromatic_by_original_photo_component"].items()}
            for name,blocks in output["results"].items()
          },sort_keys=True,allow_nan=False))

if __name__=="__main__":
    main()
