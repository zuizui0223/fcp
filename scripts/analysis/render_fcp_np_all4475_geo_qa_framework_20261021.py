#!/usr/bin/env python3
"""Additional source-locked New Phytologist figures: all-4475 completeness,
50-km geographical allocation, and H1–H5 explanatory-stage boundaries.

No new biological measurements. Nature of missingness, nulls, and each
hypothesis' empirical versus theoretical status must stay visible.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import xml.etree.ElementTree as ET

from render_fcp_np_all4475_figures_20261021 import (
    BLUE, GREY, INK, PALE, ORANGE,
    check, txt, line, rect, dot, wrap, val,
)

COHORTS=(
    ("discovery","Discovery"),
    ("validation","Validation"),
    ("prospective","Prospective"),
    ("wave4","Wave 4"),
    ("new_U100","New U100"),
)
EXPECTED_SHA="8e7f10f0427927208485b2da6fa9c46bb681d6f3611e53b8dfdc606cfc8dabe6"


def need_sources(core,missing,D,white):
    r=core["result"]
    check(core["matched_photo_csv_sha256"]==EXPECTED_SHA,"photo source mismatch")
    check(core["result"]["all_taxa"]==4475 and core["result"]["source_classified_photo_rows"]==227107,
          "entire-photo source census lost taxa or changed classifiability")
    check(missing["n_taxa"]==4475 and missing["source_actual_447500_photo_sha256"]==EXPECTED_SHA,
          "missingness source differs from 4475 authentic photos")
    check(D["n_taxa"]==4475 and D["n_colour_usable"]==227107,
          "D-classifiability source numbers changed")
    check(white["result"]["H1_W"]["all_4475"]["n_species"]==4475,
          "whole-frame W composition null omitted source taxa")
    check(r["n_h2_any_local50_pair"]==3600
          and r["n_h3_two_colours_with_complete_BIO5"]==4112,
          "all-4475 H2/H3 observed denominator changed")
    check(r["n_h4_genetic_coexistence_verified"]==0
          and r["n_H5_individual_matched_molecular_phenotypes"]==0,
          "new genetic or chemical evidence now exists: redesign inferred ceiling")


def render_geo(core):
    z=core["result"]
    p=z["H2_pair_opportunity_sensitivity"]
    check(p["at_least_1_local_photo_pair"]["n_species"]==3600,
          "at-least-one 50km source availability changed")
    cohort=z["by_species_disjoint_cohort"]
    W,H=1400,715
    b=[txt(55,51,"Geographic colour allocation is recurrent, not universal",size=28,weight="bold"),
       txt(55,84,"4,475 taxa retained; 3,600 have at least one same-species ≤50-km photo pair",size=18,color=GREY)]
    b.append(line(700,114,700,595,color="#d6e1e9",w=2))
    b.append(txt(56,130,"A  Five independent species cohorts",size=20,weight="bold"))
    b.append(txt(740,130,"B  Photo-pair opportunity thresholds",size=20,weight="bold"))
    leftX=lambda r:300+(r+.01)/.05*295
    b.append(line(leftX(0),160,leftX(0),515,color="#9fabb2",dash="5,5",w=2))
    for i,(name,label) in enumerate(COHORTS):
        c=cohort[name];n=int(c["n_H2_local"]);value=float(c["H2_mean"])
        check(0<n<=c["n_species"] and abs(value)<.2,"corrupt signed photo mismatch")
        y=185+70*i
        b.append(txt(56,y+4,label,size=19))
        b.append(line(leftX(0),y,leftX(value),y,w=4,color=BLUE))
        b.append(dot(leftX(value),y,r=9,fill=BLUE))
        b.append(txt(618,y+4,f"{value:+.4f}",size=17,anchor="end",weight="bold"))
        b.append(txt(56,y+26,f"{n:,} eligible taxa",size=14,color=GREY))
    rightX=lambda r:1050+(r+.005)/.035*170
    b.append(line(rightX(0),160,rightX(0),515,color="#9fabb2",dash="5,5",w=2))
    for i,(threshold,label) in enumerate(((1,"≥1 pair"),(5,"≥5 pairs"),(30,"≥30 pairs"),(100,"≥100 pairs"))):
        name=f"at_least_{threshold}_local_photo_pair"+("s" if threshold!=1 else "")
        if name not in p:
            name=f"at_least_{threshold}_local_photo_pairs"
        check(name in p,"pair-opportunity support missing for "+label)
        x=p[name];n=int(x["n_species"]);value=float(x["observed_mean_depletion"])
        check(x["upper_p_199"]==.005,"threshold null permutation evidence changed")
        y=190+i*93
        b.append(txt(738,y-7,label,size=17,weight="bold"))
        b.append(txt(738,y+15,f"{n:,} species",size=14,color=GREY))
        b.append(line(rightX(0),y,rightX(value),y,color=BLUE,w=4))
        b.append(dot(rightX(value),y,r=9,fill=BLUE))
        b.append(txt(1260,y+6,f"{value:+.4f}",size=17,weight="bold"))
    b.append(txt(55,576,"Source-conditional mismatch depletion (range-wide minus local; species-equal).",size=17,color=INK))
    b.append(txt(55,604,"199 same-species composition-and-coordinate-preserving null worlds; pooled P = 0.005.",size=16,color=GREY))
    b.append(txt(55,629,"875 taxa have no local photo-pair opportunity and stay NA, not zero.",size=16,color=GREY))
    b.append(txt(55,654,"Neither local inherited morph coexistence nor historical demography is inferred.",size=16,color=GREY))
    return wrap(W,H,b)


def render_missingness(core,missing,D):
    z=core["result"]; W,H=1450,710
    b=[txt(55,52,"All 4,475 taxa retained: observational support is not biological prevalence",size=27,weight="bold"),
       txt(55,85,"447,500 original photos • 227,107 classifiable labels • 220,393 unresolved labels",
           size=18,color=GREY)]
    b.append(txt(55,144,"A  Distinct empirically estimable measures",size=21,weight="bold"))
    b.append(txt(750,144,"B  The effect of unclassified photographs on D",size=21,weight="bold"))
    b.append(line(715,157,715,560,color="#d6e1e9",w=2))
    items=[
        ("H1: colour-axis geometry",int(z["n_h1_continuous_axis_evaluable"])),
        ("H2: local ≤50-km colour pairs",int(z["n_h2_any_local50_pair"])),
        ("H3: same-species white/BIO5",int(z["n_h3_two_colours_with_complete_BIO5"])),
    ]
    for i,(label,n) in enumerate(items):
        y=202+i*113
        b.append(txt(55,y-17,label,size=20))
        b.append(rect(55,y,525,22,fill="#e1e9ed",rad=3))
        b.append(rect(55,y,525*n/4475,22,fill=BLUE,rad=3))
        b.append(txt(605,y+20,f"{n:,}/4,475",size=19,color=BLUE,weight="bold"))
    b.append(txt(55,556,"H4/H5: 0 individuals with linked inherited-morph",size=15,color=ORANGE))
    b.append(txt(55,576,"and tissue-genotype-fitness evidence; NOT zero morphs.",size=15,color=ORANGE))
    rmin=float(missing["mean_D_hypothetical_all100_sharp_lower"])
    rmax=float(missing["mean_D_hypothetical_all100_sharp_upper"])
    actual=float(missing["mean_D_conditional_on_classified"])
    check(abs(actual-D["mean_D_all"])<1e-10,"sharp-bound D and observed D disagree")
    check(0<=rmin<actual<rmax<=1,"source missingness bounds invalid")
    x=lambda a:765+float(a)/.8*560
    y=295
    b.append(rect(x(rmin),y-32,x(rmax)-x(rmin),66,fill="#e1e9ed",rad=4))
    b.append(line(x(actual),y-43,x(actual),y+47,w=4,color=BLUE))
    b.append(dot(x(actual),y,r=11,fill=BLUE))
    b.append(txt(x(rmin),230,f"{rmin:.3f}",size=20,color=GREY,weight="bold",anchor="middle"))
    b.append(txt(x(actual),395,f"{actual:.3f}",size=21,color=BLUE,weight="bold",anchor="middle"))
    b.append(txt(x(rmax),230,f"{rmax:.3f}",size=20,color=GREY,weight="bold",anchor="middle"))
    b.append(txt(763,430,"Grey: sharp bounds on hypothetical four-colour assignment",size=17,color=GREY))
    b.append(txt(763,458,"Blue: conditional D from the classified photo subset",size=17,color=BLUE))
    b.append(txt(763,502,"They are NOT statistical CIs or bounds on genetic morph diversity.",size=16,color=ORANGE))
    b.append(txt(55,631,"Every taxon has an original 100-photo row; NA means a contrast cannot be estimated.",size=17,color=GREY))
    b.append(txt(55,658,"Unclassified images may not contain an identifiable flower organ, so all-100 bounds are hypothetical.",size=15,color=GREY))
    return wrap(W,H,b)


def render_framework(core,white):
    w=float(white["result"]["H1_W"]["all_4475"]["observed_mean_W"])
    null=float(white["result"]["H1_W"]["all_4475"]["null_median_W"])
    check(abs(w-.478871911)<1e-4,"source W changed and requires an updated concept")
    W,H=1450,800;b=[]
    b.extend([txt(55,52,"A common phenotype direction, three separable evolutionary stages",size=28,weight="bold"),
              txt(55,86,"H1–H5 distinguish accessibility, geographic allocation, and inherited persistence",size=19,color=GREY)])
    boxes=[
      (55,155,415,278,"Stage 1 · Phenotypic accessibility",
       "H1: recurrent colour direction",
       f"W = {w:.4f} vs null {null:.4f}",
       "4,475 photographed source taxa",
       "Positive source-conditioned evidence"),
      (518,155,415,278,"Stage 2 · Geographic allocation",
       "H2: historical isolation (unknown)",
       "H3: environmental response (limited)",
       "3,600/4,112 estimable photo taxa",
       "Source correlations ≠ causal history"),
      (981,155,415,278,"Stage 3 · Inherited maintenance",
       "H4: local genetic coexistence (unknown)",
       "H5: tissue pigment costs (unknown)",
       "No matched genotype/fitness data",
       "Published PAL/WAL evidence is external"),
    ]
    for i,(x,y,wid,h,title,a,c,d,e) in enumerate(boxes):
        b.append(rect(x,y,wid,h,fill="#ecf3f6" if i==0 else "#f3f6f7",rad=12))
        b.append(rect(x,y,8,h,fill=BLUE if i==0 else GREY,rad=2))
        b.append(txt(x+27,y+47,title,size=20,weight="bold"))
        b.append(txt(x+27,y+103,a,size=17))
        b.append(txt(x+27,y+147,c,size=17))
        b.append(txt(x+27,y+195,d,size=16,color=GREY))
        b.append(txt(x+27,y+244,e,size=15,color=ORANGE if i>0 else BLUE))
    for xx in (475,938):
        b.append(line(xx,294,xx+36,294,color=GREY,w=2,dash="6,6"))
    b.append(txt(55,492,"Distinguishing predictions for independent future tests",size=23,weight="bold"))
    rows=[
      ("H2 | HISTORY","Same present environment, different neutral ancestry → different colour frequencies?"),
      ("H3 | SELECTION","Same neutral ancestry, contrasting environment → genotype × treatment fitness difference?"),
      ("H4–H5 | MAINTENANCE","Same population, multi-year true morphs → frequency-dependent fitness and petal/leaf pleiotropy?"),
    ]
    for i,(head,msg) in enumerate(rows):
        yy=537+i*66
        b.append(rect(55,yy-27,1332,52,fill="#f4f7f9",rad=4))
        b.append(txt(76,yy+4,head,size=18,color=BLUE,weight="bold"))
        b.append(txt(322,yy+4,msg,size=16))
    b.append(txt(55,765,"Dashed links are hypotheses, not fitted evolutionary causal arrows; 447,500 photo images contain no neutral SNP or seed-fitness assays.",size=15,color=GREY))
    return wrap(W,H,b)


def synthetic():
    x=render_framework({"result":{}},{"result":{"H1_W":{"all_4475":{"observed_mean_W":.478871911,"null_median_W":.432039123}}}})
    ET.fromstring(x)
    check("not fitted evolutionary" in x,"conceptual boundary missing")
    print("FCP_NP_PHENOTYPE_GEOGRAPHY_MISSINGNESS_VISUAL_SYNTHETIC_PASS")


def main():
    p=argparse.ArgumentParser()
    for k in ("core","white","missing","D","out"):
        p.add_argument("--"+k,type=Path)
    p.add_argument("--self-test",action="store_true")
    a=p.parse_args()
    synthetic()
    if a.self_test:return
    check(a.core is not None and a.white is not None and a.missing is not None
          and a.D is not None and a.out is not None,"figure input missing")
    read=lambda path:json.loads(path.read_text())
    core,white,missing,D=(read(x) for x in (a.core,a.white,a.missing,a.D))
    need_sources(core,missing,D,white)
    graphics={
       "fcp_np_all4475_figure2_geo_50km_20261021.svg":render_geo(core),
       "fcp_np_all4475_supp_photo_missingness_20261021.svg":render_missingness(core,missing,D),
       "fcp_np_all4475_figure5_mechanism_schema_20261021.svg":render_framework(core,white),
    }
    a.out.mkdir(parents=True,exist_ok=True)
    for fname,s in graphics.items():
        check(len(s)>3000 and s.endswith("</svg>\n"),"bad SVG source length")
        root=ET.fromstring(s)
        check(root.tag.endswith("svg"),"invalid SVG root")
        (a.out/fname).write_text(s,encoding="utf-8")
        print("FCP_NP_ALL4475_FIGURE_SOURCE_QA_PASS",fname,len(s))
    print("FCP_NP_ALL4475_GEO_MISSINGNESS_FRAMEWORK_FIGURES_COMPLETE")


if __name__=="__main__":
    main()
