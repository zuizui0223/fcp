#!/usr/bin/env python3
"""Source-backed figures; chart styling without externally fetched images."""
import argparse
import json
from pathlib import Path
import matplotlib.pyplot as plt


def figures(root):
    root=Path(root);data=json.loads((root/'results'/'result.json').read_text())
    out=root/'docs';out.mkdir(exist_ok=True)
    keys=['BASE','CONGEN','BOTH','DIST']
    names=['Site mean only','Other white species','Distance + other','Distance only']
    y=[data['main']['mean']['logloss_'+k] for k in keys]
    fig,ax=plt.subplots(figsize=(9,4.5))
    ax.barh(names,y)
    ax.invert_yaxis()
    ax.set_xlim(0,0.85)
    ax.set_xlabel('Held-out M-individual binomial log loss (lower is better)')
    ax.set_title('Predicting purple/white flowers in unseen field lines: six sites')
    for i,s in enumerate(y):ax.text(s+0.008,i,f'{s:.4f}',va='center')
    fig.tight_layout()
    fig.savefig(out/'heldout_line_model_comparison.png',dpi=170)
    plt.close(fig)
    sites=data['selected_latest']
    names=[s['siteyear'] for s in sites]
    g=[s['extra_congen_beyond_distance'] for s in sites]
    fig,ax=plt.subplots(figsize=(9,4.8))
    ax.barh(names,g)
    ax.invert_yaxis()
    ax.axvline(0,color='black',lw=1)
    ax.set_xlabel('Incremental gain of other-species share beyond distance\nDIST log loss - BOTH log loss (positive = improvement)')
    ax.set_title('Incremental predictor utility is not consistent among six sites')
    fig.tight_layout()
    fig.savefig(out/'heldout_line_congener_increment_by_site.png',dpi=170)
    plt.close(fig)
    return [str(out/'heldout_line_model_comparison.png'),str(out/'heldout_line_congener_increment_by_site.png')]

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--root',default='.');print(figures(a.parse_args().root))
