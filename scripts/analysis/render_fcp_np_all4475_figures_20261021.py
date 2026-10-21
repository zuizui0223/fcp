#!/usr/bin/env python3
"""Publication-grade source-locked, coordinate-free SVG summaries for FCP NP.

All visible points and uncertainty intervals are rendered from verified actual
447500-photo analysis receipts committed to the same Git analysis branch.
No fabricated significance, colour morph genotype, phylogeny or geography.
Uses no font file and no opaque chart package state.
"""
from __future__ import annotations
import argparse
import json
import math
from pathlib import Path
from xml.sax.saxutils import escape

BLUE="#226376"
GREY="#7a8491"
INK="#233548"
PALE="#e7eef2"
ORANGE="#a3603f"
BG="#ffffff"
COHORTS=[
 ("all_4475","All taxa • 4,475"),
 ("cohort_discovery","Discovery • 500"),
 ("cohort_validation","Validation • 500"),
 ("cohort_prospective","Prospective • 499"),
 ("cohort_wave4","Wave 4 • 498"),
 ("cohort_new_U100","New U100 • 2,478"),
]
AXES=[
 ("BIO1","Annual mean temperature"),
 ("BIO4","Temperature seasonality"),
 ("BIO5","Warmest-month maximum"),
 ("BIO12","Annual precipitation"),
 ("BIO14","Driest-month precipitation"),
 ("BIO15","Precipitation seasonality"),
 ("SRADmean","Mean solar radiation"),
]
PAIRS=[
 ("W_geometry_vs_geographic_sorting","W axis × locality sorting"),
 ("W_geometry_vs_photo_BIO5_slope","W axis × BIO5 response"),
 ("W_geometry_vs_photo_temperature_seasonality","W axis × temperature seasonality"),
 ("geographic_sorting_vs_photo_BIO5_slope","Locality sorting × BIO5 response"),
]


def check(ok,msg):
    if not ok:
        raise RuntimeError("FCP_NP_FIGURE_SOURCE_INVALID: "+msg)


def val(x):
    return f"{float(x):.4g}"


def txt(x,y,phrase,*,size=18,color=INK,weight=None,anchor="start",italic=False):
    attr=f' font-weight="{weight}"' if weight else ""
    sty=' font-style="italic"' if italic else ""
    return f'<text x="{x:.1f}" y="{y:.1f}" font-family="Arial,Helvetica,sans-serif" font-size="{size}" fill="{color}" text-anchor="{anchor}"{attr}{sty}>{escape(str(phrase))}</text>'


def line(x1,y1,x2,y2,*,color=INK,w=2,dash=None):
    d=f' stroke-dasharray="{dash}"' if dash else ""
    return f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" stroke="{color}" stroke-width="{w}"{d}/>'


def rect(x,y,w,h,*,fill=PALE,rad=0):
    return f'<rect x="{x:.1f}" y="{y:.1f}" width="{w:.1f}" height="{h:.1f}" fill="{fill}" rx="{rad}"/>'


def dot(x,y,*,r=9,fill=BLUE,stroke=None):
    outer=f' stroke="{stroke}" stroke-width="2"' if stroke else ""
    return f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{r}" fill="{fill}"{outer}/>'


def wrap(w,h,body):
    return f'<svg xmlns="http://www.w3.org/2000/svg" role="img" width="{w}" height="{h}" viewBox="0 0 {w} {h}">'+rect(0,0,w,h,fill=BG)+''.join(body)+"</svg>\n"


def axes_ticks(parts,zero,line_y,mapping,ticks):
    for z in ticks:
        x=mapping(z)
        parts.append(line(x,line_y-4,x,line_y+5,color="#8d9ca8",w=1))
        parts.append(txt(x,line_y+29,f"{z:+.2f}",size=15,color=GREY,anchor="middle"))
    parts.append(line(mapping(zero),132,mapping(zero),line_y-5,color="#91a7b5",w=2,dash="5,5"))


def render_fate(root):
    z=root["results"]
    check(root["all_species_in_common_parent_roster"]==4475,"fate source not complete 4475")
    W,H=1350,650;pieces=[]
    pieces+=[txt(60,54,"Repeated flower-colour geometry does not predict a uniform fate",
                 size=28,weight="bold"),txt(60,85,"Source-conditional rank associations; original 4,475 species",
                 size=17,color=GREY)]
    L,R=620,1100
    x=lambda a:L+(float(a)+.12)/.24*(R-L)
    axes_ticks(pieces,0,532,x,[-.10,-.05,0,.05,.10])
    for i,(key,label) in enumerate(PAIRS):
        item=z[key]
        check(item["n_species"]>100 and item["partial_rank_r_controlling_photo_opportunity_and_source"] is not None,
              "fate association missing "+key)
        ci=item["genus_bootstrap_95CI_conditional_fixed_residuals"]
        check(ci and ci[0]<=ci[1],"genus CI missing")
        r=item["partial_rank_r_controlling_photo_opportunity_and_source"]
        y=175+i*90
        if i==3:
            pieces.append(line(60,y-45,1270,y-45,color="#d8e3e9",w=2))
        pieces.append(txt(60,y-1,label,size=21,weight="bold" if i<3 else None))
        pieces.append(txt(60,y+24,f"{item['n_species']:,} source taxa with estimable association",
                          size=15,color=GREY))
        pieces.append(line(x(ci[0]),y,x(ci[1]),y,color=BLUE if i<3 else ORANGE,w=5))
        pieces.append(line(x(ci[0]),y-8,x(ci[0]),y+8,color=BLUE if i<3 else ORANGE,w=2))
        pieces.append(line(x(ci[1]),y-8,x(ci[1]),y+8,color=BLUE if i<3 else ORANGE,w=2))
        pieces.append(dot(x(r),y,fill=BLUE if i<3 else ORANGE,r=9))
        pieces.append(txt(1150,y+6,f"{r:+.3f}",size=21,weight="bold",color=BLUE if i<3 else ORANGE))
    pieces.append(txt(620,582,"Adjusted partial rank association",size=17,anchor="middle"))
    pieces.append(txt(60,614,"95% genus-cluster intervals from fixed nuisance residuals; exploratory and source-overlapping.",size=15,color=GREY))
    pieces.append(txt(60,636,"No test of causal history, inherited morph maintenance, fitness or molecular origins.",size=15,color=GREY))
    return wrap(W,H,pieces)


def render_white(data):
    z=data["result"]["H1_W"] if "result" in data else data["H1_W"]
    check(z["all_4475"]["n_species"]==4475,"white axis complete source mismatch")
    W,H=1350,650;parts=[]
    parts+=[txt(60,54,"White–chromatic phenotype alignment recurs across all source cohorts",size=28,weight="bold"),
            txt(60,83,"Observed W versus the species/colour-state-preserving, 199-world null",size=17,color=GREY)]
    L,R=540,1140
    x=lambda a:L+(float(a)-.37)/(.56-.37)*(R-L)
    axes_ticks(parts,.45,555,x,[.38,.42,.46,.50,.54])
    for i,(key,name) in enumerate(COHORTS):
        v=z[key]
        lo,hi=v["null_q025_q975"]
        med=v["null_median_W"];obs=v["observed_mean_W"]
        check(lo<=med<=hi and v["n_null_worlds"]==199,"null cannot be plotted for "+key)
        y=155+i*67
        if i==1:parts.append(line(60,y-30,1280,y-30,color="#dae3e8",w=2))
        parts.append(txt(60,y+5,name,size=20))
        parts.append(rect(x(lo),y-12,x(hi)-x(lo),25,fill="#dae2e7",rad=3))
        parts.append(line(x(med),y-16,x(med),y+16,color=GREY,w=2))
        parts.append(dot(x(obs),y,r=8,fill=BLUE))
        parts.append(txt(1170,y+5,val(obs),size=18,weight="bold"))
    parts.append(txt(540,595,"W: mean squared alignment with the fixed nine-colour white axis",size=16))
    parts.append(txt(60,621,"Blue: observed  •  Grey band/line: 95% null range / null median. All six upper p = 0.005.",size=15,color=GREY))
    parts.append(txt(60,643,"Exploratory shared-provider result; the separate frozen prospective 158-species test remains unchanged.",size=14,color=GREY))
    return wrap(W,H,parts)


def render_seven(data):
    z=data.get("axes") or data["result"]["axes"]
    order={j["axis"]:j for j in z}
    check(len(z)==7 and all(k in order for k,_ in AXES),"seven-axis source test incomplete")
    W,H=1350,675;parts=[]
    parts+=[txt(60,54,"Seven abiotic gradients: only two source-photo associations remain",size=28,weight="bold"),
            txt(60,83,"White-minus-nonwhite location contrasts; genus-adjusted simultaneous 95% intervals",size=17,color=GREY)]
    L,R=590,1130
    x=lambda z:L+(float(z)+.08)/.17*(R-L)
    axes_ticks(parts,0,565,x,[-.06,-.03,0,.03,.06,.09])
    for i,(key,name) in enumerate(AXES):
        z=order[key]
        lo,hi=z["genus_bootstrap_family7_maxT_95"]
        r=z["mean_species_equal_white_minus_nonwhite_SD"]
        check(z["n_species_estimable"]==4112,"seven-axis availability changed")
        y=150+i*57
        parts.append(txt(62,y+7,f"{key}  {name}",size=18,weight="bold" if z["family7_excludes_zero"] else None))
        c=BLUE if z["family7_excludes_zero"] else GREY
        parts.append(line(x(lo),y,x(hi),y,color=c,w=5 if z["family7_excludes_zero"] else 3))
        parts.append(line(x(lo),y-9,x(lo),y+9,color=c,w=2))
        parts.append(line(x(hi),y-9,x(hi),y+9,color=c,w=2))
        parts.append(dot(x(r),y,fill=c,r=8))
        parts.append(txt(1180,y+6,f"{r:+.3f}",size=19,weight="bold" if z["family7_excludes_zero"] else None,color=c))
    parts.append(txt(590,608,"Within-species climate SD units (white − nonwhite photo localities)",size=17,anchor="middle"))
    parts.append(txt(60,637,"4,112 estimable of 4,475 source taxa; simultaneous seven-axis max-|t| genus bootstrap.",size=15,color=GREY))
    parts.append(txt(60,659,"Post hoc WorldClim normals; within-year BIO4 is not interannual selection or an inherited pigment effect.",size=14,color=GREY))
    return wrap(W,H,parts)


def synthetic():
    x=wrap(100,100,[txt(15,15,"A&B",size=14),line(10,10,40,40)])
    check("&amp;" in x and x.endswith("</svg>\n"),"SVG escaping broken")
    print("FCP_NP_ALL4475_SVG_AUTHORING_SYNTHETIC_PASS")


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--joint",type=Path)
    ap.add_argument("--white",type=Path)
    ap.add_argument("--seven",type=Path)
    ap.add_argument("--out",type=Path)
    ap.add_argument("--self-test",action="store_true")
    a=ap.parse_args()
    synthetic()
    if a.self_test:return
    check(a.joint and a.white and a.seven and a.out,"missing source graphic materials")
    read=lambda p:json.loads(p.read_text())
    figures={
      "fcp_np_all4475_geometry_vs_fate_20261021.svg":render_fate(read(a.joint)),
      "fcp_np_all4475_white_axis_five_cohorts_20261021.svg":render_white(read(a.white)),
      "fcp_np_all4475_seven_abiotic_axes_20261021.svg":render_seven(read(a.seven)),
    }
    a.out.mkdir(parents=True,exist_ok=True)
    for name,value in figures.items():
        check(len(value)>1800 and "<svg" in value and value.strip().endswith("</svg>"),
              "figure construction invalid "+name)
        (a.out/name).write_text(value)
        print("FCP_NP_REPRODUCIBLE_SVG_READY",name,len(value))
    print("FCP_NP_ALL4475_REPRODUCIBLE_THREE_FIGURES_COMPLETE")

if __name__=="__main__":main()
