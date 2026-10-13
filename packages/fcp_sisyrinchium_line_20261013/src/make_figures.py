#!/usr/bin/env python3
"""Render source-locked site-year observational effects; no significance by plot itself."""
from pathlib import Path
import argparse,json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def make(result,output):
    r=json.loads(Path(result).read_text())
    out=Path(output);out.mkdir(exist_ok=True,parents=True)
    main=r['full']['results']
    site=main['exposure']['names'];color=[main['exposure']['siteyear_estimates'][s]*100 for s in site]
    adjusted=[main['distance']['siteyear_estimates'][s]*100 for s in site]
    fig,ax=plt.subplots(figsize=(9.2,4.9))
    yy=list(range(len(site)))
    ax.scatter(color,yy,label='Observed congener-share contrast (within line)',marker='o',s=72)
    ax.scatter(adjusted,yy,label='After removing linear water-distance trend',marker='s',s=55)
    ax.axvline(0,color='.6',linewidth=1)
    ax.set_yticks(yy,labels=site);ax.invert_yaxis()
    ax.set_xlabel('Other white species fraction: violet M minus white M (percentage points)')
    ax.set_title('Same-site, same-line association and distance-sensitive contrast')
    ax.legend(fontsize=8,loc='lower right')
    fig.tight_layout();fig.savefig(out/'sisyrinchium_line_congener_exposure.png',dpi=170);plt.close(fig)
    sites=main['clumping']['names'];vals=[main['clumping']['siteyear_estimates'][s]*100 for s in sites]
    fig,ax=plt.subplots(figsize=(9.2,5.0))
    ax.barh(sites,vals)
    ax.invert_yaxis();ax.axvline(0,color='.5',linewidth=1)
    ax.set_xlabel('Observed - within-recorded-line fixed-inventory expected colour diversity (%)')
    ax.set_title('Residual M-morph quadrat segregation inside recorded field lines')
    fig.tight_layout();fig.savefig(out/'sisyrinchium_line_residual_clumping.png',dpi=170);plt.close(fig)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--result',default='results/result.json');parser.add_argument('--out',default='docs');a=parser.parse_args()
    make(a.result,a.out)
