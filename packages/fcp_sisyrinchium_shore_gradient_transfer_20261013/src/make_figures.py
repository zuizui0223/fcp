#!/usr/bin/env python3
"""Source-backed illustrations from deterministic within-site heldout results."""
import argparse,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def plot(root):
    root=Path(root);out=root/'docs';out.mkdir(exist_ok=True)
    result=json.loads((root/'results/result.json').read_text())
    extra=json.loads((root/'results/posthoc_distance_scale.json').read_text())
    sites=sorted(result['sites'],key=lambda x:x['siteyear'])
    labels=[s['siteyear'].split('|')[1]+' ('+s['siteyear'].split('|')[0]+')' for s in sites]
    x=np.arange(len(sites))
    fig,ax=plt.subplots(figsize=(10.5,5.4),dpi=165)
    slope=np.array([s['full_site_slope_raw'] for s in sites]);ax.barh(x,slope)
    ax.set_yticks(x,labels);ax.axvline(0,linestyle='--',linewidth=1)
    ax.invert_yaxis()
    ax.set_xlabel('Full-site logit gradient per recorded water-distance unit')
    ax.set_title('Shoreline-gradient direction varies among source populations')
    for i,v in enumerate(slope):ax.text(v+(0.008 if v>=0 else -0.008),i,f'{v:+.3f}',ha='left' if v>=0 else 'right',va='center',fontsize=9)
    fig.tight_layout();p=out/'full_site_distance_slopes.png';fig.savefig(p);plt.close(fig)
    fig,ax=plt.subplots(figsize=(11.5,5.6),dpi=165)
    width=.22
    post={s['siteyear']:s for s in extra['sites']}
    values=[
        [s['logloss_BASE'] for s in sites],
        [s['logloss_LOCAL_DIST'] for s in sites],
        [s['logloss_BORROW_MEDIAN'] for s in sites],
        [post[s['siteyear']]['score']['BORROW_Z_MEDIAN'] for s in sites],
    ]
    for j,(vals,lbl) in enumerate(zip(values,['Baseline site mean','Local distance','Borrowed raw-distance slope','Borrowed standardized slope (posthoc)'])):
        ax.bar(x+(j-1.5)*width,vals,width,label=lbl)
    ax.set_xticks(x,labels,rotation=20,ha='right')
    ax.set_ylabel('Heldout-line binomial log loss / observed M plant (lower better)')
    ax.set_title('Same 32 source lines; portability depends on distance scaling')
    ax.legend(frameon=False,fontsize=8)
    fig.tight_layout();p2=out/'heldout_cross_site_gradient_transport.png';fig.savefig(p2);plt.close(fig)
    return [str(p),str(p2)]

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--root',default='.');p=a.parse_args();print(plot(p.root))
