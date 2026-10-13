#!/usr/bin/env python3
"""Source-backed charts from frozen site-year morph allocations. No invented coordinates."""
import json, argparse
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

def plot(root):
    root=Path(root);r=json.loads((root/'results/result.json').read_text())
    main=r['median'];names=main['names']
    base=json.loads((root/'source/parent_line_result.json').read_text())
    old=base['full']['results']['clumping']['siteyear_estimates']
    a=[old[x]*100 for x in names];b=[main['siteyear_estimates'][x]*100 for x in names]
    fig,ax=plt.subplots(figsize=(10,5.4));y=np.arange(len(names))
    ax.barh(y-0.21,a,height=.4,label='Original line-conditioned excess')
    ax.barh(y+0.21,b,height=.4,label='Near/far distance-conditioned excess')
    ax.set_yticks(y,names);ax.invert_yaxis();ax.set_xlabel('Excess among-quadrat phenotypic diversity (percentage points)')
    ax.axvline(0,linewidth=.7);ax.legend(loc='lower right',fontsize=8)
    fig.tight_layout();fig.savefig(root/'docs/within_line_distance_heterogeneity_by_site.png',dpi=170);plt.close(fig)
    labels=['Original observed B','Null expectation (fixed distance bins)','Residual excess']
    values=[main['mean_observed_B']*100,main['mean_exact_null_B']*100,main['mean_delta']*100]
    fig,ax=plt.subplots(figsize=(7.5,4.2));bars=ax.barh(labels,values);ax.invert_yaxis()
    for bar,v in zip(bars,values):ax.text(v+0.5,bar.get_y()+bar.get_height()/2,f'{v:.2f}%',va='center')
    ax.set_xlim(0,max(values)*1.23);ax.set_xlabel('Between-quadrat phenotypic diversity (%)')
    fig.tight_layout();fig.savefig(root/'docs/observed_vs_distance_matched_null.png',dpi=170);plt.close(fig)
    return len(names)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',default='.');a=p.parse_args();print('PLOTS',plot(a.root))
