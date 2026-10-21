"""Create reproducible preview figures from independently source-replayed FCP PR #193/#194 CSVs.

These are *manuscript preview* figures, not biological validations or final
submission figure files. All displayed values are read from original result rows.

Example:
  python make_prototype_figures.py \
    --macro-ledger /path/to/PR193/direction_ledger.csv \
    --matched-ledger /path/to/PR194/matched_spatial_support_differences.csv \
    --outdir /tmp/figures
"""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

AXES = [('wc_bio12','Annual precipitation (BIO12)'),
        ('wc_vapr_annual_mean','Atmospheric vapour pressure'),
        ('wc_bio5','Warmest-month temperature (BIO5)'),
        ('wc_srad_annual_mean','Annual solar radiation')]
HUES = [('white','white'),('yellow_orange','yellow/orange'),
        ('red_pink','red/pink'),('blue_purple','blue/purple')]


def source_row(table, **filters):
    q = table.copy()
    for key,value in filters.items():
        q=q.loc[q[key]==value]
    if len(q)!=1:
        raise ValueError(f'Expected one source row, got {len(q)} for {filters}')
    return q.iloc[0]


def fig2(macro, outdir):
    dat=[]
    for anchor,axis_label in AXES:
        for hue,hue_label in HUES:
            row=source_row(macro,mode='climate',quality='all',region='world',
                           adjustment='genus_geo',weight='genus_equal',anchor=anchor,hue=hue)
            dat.append((f'{axis_label}  |  {hue_label}',
                        float(row['direction_pp_per_env_global_sd']),
                        float(row['maxT_ci_lo_pp']),float(row['maxT_ci_hi_pp'])))
    labels,centres,los,his = zip(*dat)
    centres=np.array(centres);los=np.array(los);his=np.array(his)
    positions=np.arange(len(dat))[::-1]
    fig,ax=plt.subplots(figsize=(10.1,8.9))
    for i in range(len(dat)):
        err=np.array([[centres[i]-los[i]],[his[i]-centres[i]]])
        ax.errorbar(centres[i],positions[i],xerr=err,fmt='o',capsize=2.9,markersize=4)
    ax.axvline(0,ls='--',lw=0.9)
    ax.set_yticks(positions,labels,fontsize=8.1)
    ax.set_xlabel('Change in photographed hue probability (pp per original global environmental SD)')
    ax.set_title('Global genus-conditioned flower-colour directions',loc='left')
    ax.tick_params(axis='x',labelsize=9)
    ax.grid(axis='x',alpha=.15)
    fig.text(.07,.019,
       'Within-cohort maxT 95% intervals across six climate axes x four hues; these four hue probabilities sum to one.\n'
       'Retrospective photographed outcomes; taxonomic and measurement limits apply. Source: FCP #193.',fontsize=8)
    fig.subplots_adjust(left=.43,right=.98,top=.93,bottom=.12)
    for ext in ('png','svg'):
        fig.savefig(outdir/f'FIG2_macro_four_hue_direction_preview.{ext}',dpi=250)
    plt.close(fig)


def fig3(macro,matched,outdir):
    items=[]
    for region,region_name in [('world','World'),('E90_E180','East 90-180 E')]:
        full=source_row(macro,mode='climate',quality='all',region=region,adjustment='genus_geo',
                        weight='genus_equal',anchor='wc_bio12',hue='white')
        strict=source_row(macro,mode='climate',quality='strict',region=region,adjustment='genus_geo',
                          weight='genus_equal',anchor='wc_bio12',hue='white')
        b=source_row(matched,region=region,weight='genus_equal',type='stage',
                     comparison='same_species_with_broad_geography',anchor='wc_bio12',hue='white')
        c=source_row(matched,region=region,weight='genus_equal',type='stage',
                     comparison='same_species_with_genus_cell',anchor='wc_bio12',hue='white')
        for kind,x in [('All classifiable',full),('Strict image quality',strict)]:
            items.append((f'{region_name}  |  {kind}',
                float(x.direction_pp_per_env_global_sd),float(x.cluster_ci_lo_pp),float(x.cluster_ci_hi_pp)))
        for kind,x in [('Fixed species, broad geography (B)',b),('Same species, cell conditioning (C)',c)]:
            items.append((f'{region_name}  |  {kind}',
                float(x.point_pp),float(x.CI95_lo_pp),float(x.CI95_hi_pp)))
    lab,mid,lo,hi=zip(*items)
    mid=np.array(mid);lo=np.array(lo);hi=np.array(hi)
    y=np.arange(len(items))[::-1]
    fig,ax=plt.subplots(figsize=(10.2,6.2))
    for i in range(len(items)):
        ax.errorbar(mid[i],y[i],xerr=np.array([[mid[i]-lo[i]],[hi[i]-mid[i]]]),
                    fmt='o',capsize=3,markersize=4.4)
    ax.axvline(0,ls='--',lw=.9)
    ax.set_yticks(y,lab,fontsize=9)
    ax.set_xlabel('Annual precipitation -> photographed white probability (pp per original source BIO12 SD)')
    ax.set_title('Geographic and photographic support changes the estimated direction',loc='left')
    ax.grid(axis='x',alpha=.15)
    fig.text(.055,.028,
      'Shown: nominal genus-bootstrap 95% intervals for individual stage slopes, not 24-comparison intervals.\n'
      'Only B vs C is species-identical. B-C 24-way maxT intervals span zero in World and East. Sources: FCP #193/#194.',fontsize=8)
    fig.subplots_adjust(left=.39,right=.98,top=.90,bottom=.18)
    for ext in ('png','svg'):
        fig.savefig(outdir/f'FIG3_BIO12_spatial_support_preview.{ext}',dpi=250)
    plt.close(fig)


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--macro-ledger',type=Path,required=True)
    p.add_argument('--matched-ledger',type=Path,required=True)
    p.add_argument('--outdir',type=Path,required=True)
    a=p.parse_args()
    a.outdir.mkdir(parents=True,exist_ok=True)
    macro=pd.read_csv(a.macro_ledger)
    matched=pd.read_csv(a.matched_ledger)
    assert macro.shape[0]==1280 and matched.shape[0]==576, (macro.shape,matched.shape)
    fig2(macro,a.outdir)
    fig3(macro,matched,a.outdir)
    print('PASS: two publication-panel prototypes derived from source-verified 1,280/576-row ledgers')

if __name__=='__main__':main()
