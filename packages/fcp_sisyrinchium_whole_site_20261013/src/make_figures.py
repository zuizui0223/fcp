#!/usr/bin/env python3
"""Plots using only source-verified frozen analysis summary, no simulated photographs."""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def make(root):
    r=json.loads((Path(root)/'results/result.json').read_text())
    d=Path(root)/'docs';d.mkdir(parents=True,exist_ok=True)
    labels=['BASE','DIST','CONGEN','BOTH'];ind=np.arange(4);width=.34
    fig,ax=plt.subplots(figsize=(8.2,5.2))
    local=[r['within_site_reference'][k] for k in labels]
    out=[r['main']['losses'][k] for k in labels]
    ax.bar(ind-width/2,local,width,label='Leave one line out (same site)')
    ax.bar(ind+width/2,out,width,label='Leave one entire site out')
    for x,v in zip(ind-width/2,local):ax.text(x,v+0.008,f'{v:.3f}',ha='center',fontsize=9)
    for x,v in zip(ind+width/2,out):ax.text(x,v+0.008,f'{v:.3f}',ha='center',fontsize=9)
    ax.set_xticks(ind,labels)
    ax.set_ylabel('Site-equal held-out binomial log loss (lower better)')
    ax.set_ylim(0,max(max(local),max(out))+0.14)
    ax.legend(frameon=False,loc='upper left')
    ax.spines[['top','right']].set_visible(False)
    fig.tight_layout();fig.savefig(d/'within_vs_whole_site_loss.png',dpi=180);plt.close(fig)
    # Dist vs BASE performance for each unseen original site.
    sites=r['site_scores'];names=[s['population'] for s in sites]
    gains=[s['gain_distance'] for s in sites]
    fig,ax=plt.subplots(figsize=(8.3,4.9))
    ax.barh(names,gains)
    ax.axvline(0,linestyle='--',lw=1)
    for i,(s,x) in enumerate(zip(sites,gains)):
        ax.text(x+(0.008 if x>=0 else -0.008),i,f'{x:+.3f}',va='center',
                ha='left' if x>=0 else 'right',fontsize=9)
    ax.set_xlim(min(gains)-.13,max(gains)+.07)
    ax.set_xlabel('BASE - DIST held-out log loss  (positive = distance helps)')
    ax.set_title('Water-distance transport into excluded original sites')
    ax.spines[['top','right']].set_visible(False)
    fig.tight_layout();fig.savefig(d/'whole_site_distance_gain_each_site.png',dpi=180);plt.close(fig)
    # The stronger out-of-range issue in one source site.
    fig,ax=plt.subplots(figsize=(8.2,4.8))
    p=np.array([s['gain_distance'] for s in sites])
    q=np.array([s['loss']['BASE']-s['loss']['CONGEN'] for s in sites])
    y=np.arange(len(sites));h=.33
    ax.barh(y-h/2,p,h,label='DIST gain')
    ax.barh(y+h/2,q,h,label='CONGEN gain')
    ax.set_yticks(y,names);ax.axvline(0,linestyle='--',lw=1)
    ax.set_xlabel('Source-site-only mean - model held-out log loss')
    ax.legend(frameon=False)
    ax.spines[['top','right']].set_visible(False)
    fig.tight_layout();fig.savefig(d/'whole_site_distance_vs_congener_site_gains.png',dpi=180);plt.close(fig)
    return [str(x) for x in sorted(d.glob('*.png'))]

if __name__=='__main__':
    import argparse
    a=argparse.ArgumentParser();a.add_argument('--root',default='.');x=a.parse_args();print('\n'.join(make(x.root)))
