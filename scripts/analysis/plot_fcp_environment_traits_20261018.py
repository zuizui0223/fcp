#!/usr/bin/env python3
"""Source-backed publication figure for FCP environmental and GIFT analyses."""
from pathlib import Path
import argparse
import hashlib
import json
from datetime import datetime, timezone

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT=Path(__file__).resolve().parents[2]
COLORS={'BIO5':'#205C83','BIO14':'#967123','SRADmean':'#85567D'}
VARIABLES=['BIO5','BIO14','SRADmean']
MODELS=['primary','observer_thinned','spatial_thinned','geographic_trend','primary_ge10']
TRAITS=['self_fertilization','lifecycle','woodiness']

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def select_one(frame,**keys):
    use=np.ones(len(frame),bool)
    for key,value in keys.items():use&=frame[key].eq(value).to_numpy()
    selected=frame.loc[use]
    if len(selected)!=1:raise ValueError(f'Expected exactly one source row: {keys}')
    return selected.iloc[0]

def style(ax):
    ax.axvline(0,color='#69747C',lw=1.1,ls=(0,(3,3)),zorder=0)
    ax.spines[['top','right','left']].set_visible(False)
    ax.spines['bottom'].set_color('#7D868E')
    ax.tick_params(axis='y',length=0,pad=12)
    ax.tick_params(axis='x',colors='#3C4852',labelsize=10)
    ax.xaxis.grid(True,color='#E8EBEE',linewidth=.7)
    ax.set_axisbelow(True)

def interval(ax,row,y,color,estimate='reml_mean',low='reml_mean_ci_low',high='reml_mean_ci_high',lw=2.5):
    m,l,h=[float(row[k]) for k in (estimate,low,high)]
    if not np.isfinite([m,l,h]).all() or not l<=m<=h:raise ValueError('Invalid source interval')
    ax.hlines(y,l,h,color=color,lw=lw,zorder=2)
    ax.plot(m,y,'o',color=color,ms=6.2,mec='white',mew=.7,zorder=4)
    return [m,l,h]

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--outdir',type=Path,default=ROOT/'results/fcp_environment_trait_synthesis_20261018')
    args=parser.parse_args();args.outdir.mkdir(parents=True,exist_ok=True)
    envpath=ROOT/'results/fcp_real_447500_environment_model_20261018/meta_summary.csv'
    giftpath=ROOT/'results/fcp_gift_modifiers_20261018/slope_modifier_contrasts.csv'
    env=pd.read_csv(envpath);gift=pd.read_csv(giftpath)
    primary=[select_one(env,model='primary',group='pooled',variable=v) for v in VARIABLES]
    sensitivity=[select_one(env,model=m,group='pooled',variable='BIO5') for m in MODELS]
    modifiers=[select_one(gift,trait=t,variable=v) for t in TRAITS for v in VARIABLES]
    assert len(modifiers)==9 and all(r.status=='estimated' for r in modifiers)
    assert len({int(r.n_species) for r in primary})==1
    assert all(int(r.bootstrap_repetitions)==5000 for r in modifiers)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.labelcolor':'#243340',
                         'text.color':'#243340','svg.fonttype':'none','pdf.fonttype':42})
    fig=plt.figure(figsize=(14,9),facecolor='white')
    a=fig.add_axes([.215,.655,.275,.185])
    b=fig.add_axes([.215,.25,.275,.26])
    c=fig.add_axes([.705,.25,.275,.59])
    for ax in (a,b,c):style(ax)
    fig.text(.045,.955,'Within-species climate associations and trait contrasts',fontsize=19,weight='bold')
    fig.text(.045,.921,'White versus chromatic photographs • Firth logistic slopes, with climate predictors fitted jointly',fontsize=11,color='#52616C')
    fig.text(.045,.872,'A',fontsize=15,weight='bold')
    fig.text(.073,.873,'Primary environmental slopes',fontsize=13,weight='bold')
    fig.text(.073,.847,f'REML means and modified-HKSJ 95% CIs • {int(primary[0].n_species):,} species',fontsize=10,color='#52616C')
    plotted=[]
    for y,row in zip([2,1,0],primary):
        vals=interval(a,row,y,COLORS[row.variable]);plotted.append({'panel':'A','source':'environment','variable':row.variable,'values':vals})
    a.set_yticks([2,1,0],['BIO5\nWarmest-month maximum\ntemperature','BIO14\nDriest-month precipitation','SRADmean\nMean solar radiation'],fontsize=10)
    a.set_ylim(-.55,2.55);a.set_xlim(-.05,.12);a.set_xticks([-.04,0,.04,.08,.12])
    a.set_xlabel('log odds per within-species climate SD',fontsize=10,labelpad=9)
    fig.text(.045,.566,'B',fontsize=15,weight='bold')
    fig.text(.073,.567,'BIO5 across analysis specifications',fontsize=13,weight='bold')
    fig.text(.073,.541,'REML means and modified-HKSJ 95% CIs; species panels differ',fontsize=10,color='#52616C')
    names=['Primary','One photo per observer','One photo per 1° cell','Geographic-trend adjustment','At least 10 photos per color']
    labels=[]
    for y,name,row in zip(range(4,-1,-1),names,sensitivity):
        vals=interval(b,row,y,COLORS['BIO5']);plotted.append({'panel':'B','source':'environment','model':row.model,'variable':'BIO5','values':vals})
        labels.append(f'{name}\nN = {int(row.n_species):,} species')
    b.set_yticks(range(4,-1,-1),labels,fontsize=10);b.set_ylim(-.6,4.6)
    b.set_xlim(-.05,.12);b.set_xticks([-.04,0,.04,.08,.12])
    b.set_xlabel('log odds per within-species climate SD',fontsize=10,labelpad=9)
    fig.text(.548,.872,'C',fontsize=15,weight='bold')
    fig.text(.576,.873,'GIFT trait contrasts in mean slopes',fontsize=13,weight='bold')
    fig.text(.576,.847,'Equal-species means: first group minus reference group',fontsize=10,color='#52616C')
    ys=[9,8,7,5,4,3,1,0,-1]
    for row,y in zip(modifiers,ys):
        color=COLORS[row.variable]
        lo,hi=float(row.family9_percentile_ci_low),float(row.family9_percentile_ci_high)
        assert lo<=row.percentile_ci_low<=row.percentile_ci_high<=hi
        c.hlines(y,lo,hi,color=color,lw=1.1,zorder=1)
        c.plot([lo,hi],[y,y],'|',color=color,ms=7,mew=1.1,zorder=2)
        vals=interval(c,row,y,color,estimate='difference_first_minus_reference',low='percentile_ci_low',high='percentile_ci_high',lw=4.2)
        plotted.append({'panel':'C','source':'GIFT','trait':row.trait,'variable':row.variable,'values':vals,'family9_bounds':[lo,hi]})
    headings=['Self-fertilization: present − absent','Life cycle: annual − perennial','Woodiness: woody − non-woody']
    for j,(heading,y) in enumerate(zip(headings,[10.05,6.05,2.05])):
        row=modifiers[j*3]
        c.text(-.45,y,heading,transform=c.get_yaxis_transform(),fontsize=10.5,weight='bold',ha='left',va='center')
        c.text(1.0,y,f'N = {int(row.n_first):,} / {int(row.n_reference):,}',transform=c.get_yaxis_transform(),fontsize=9.5,ha='right',va='center',color='#52616C')
    c.set_yticks(ys,VARIABLES*3,fontsize=10)
    c.set_ylim(-1.7,10.65);c.set_xlim(-.28,.28);c.set_xticks([-.2,-.1,0,.1,.2])
    c.set_xlabel('difference in log odds per within-species climate SD',fontsize=10,labelpad=9)
    handles=[Line2D([0],[0],color='#52616C',lw=4.2,label='Nominal 95% genus-bootstrap CI'),
             Line2D([0],[0],color='#52616C',lw=1.1,marker='|',markersize=7,label='Approx. Bonferroni-family9 interval')]
    c.legend(handles=handles,loc='upper left',bbox_to_anchor=(-.4,-.11),frameon=False,fontsize=10,handlelength=2.5,labelspacing=.8)
    fig.text(.045,.112,'A–B: approximate two-stage Firth + REML analysis. Intervals concern mean slopes; no prediction intervals are shown.',fontsize=10,color='#52616C')
    fig.text(.045,.087,'B: thinning reapplies eligibility; geographic adjustment has its own retained panel. These are not paired-effect estimates.',fontsize=10,color='#52616C')
    fig.text(.045,.062,'C: 5,000 genus-cluster resamples; nine reported contrasts. Trait coverage and species counts vary by comparison.',fontsize=10,color='#52616C')
    fig.text(.045,.037,'Within-species SDs differ in physical size. Climate associations and trait contrasts are observational, not causal effects.',fontsize=10,color='#52616C')
    png=args.outdir/'fcp_environment_traits.png';svg=args.outdir/'fcp_environment_traits.svg'
    fig.savefig(png,dpi=220,facecolor='white');fig.savefig(svg,facecolor='white',metadata={'Date':None})
    svg.write_text('\n'.join(line.rstrip() for line in svg.read_text().splitlines())+'\n')
    plt.close(fig)
    provenance={'schema':'fcp_environment_trait_synthesis_figure_v1','created_utc':datetime.now(timezone.utc).isoformat(),
       'script_sha256':sha(Path(__file__)),'matplotlib_version':matplotlib.__version__,
       'source_sha256':{str(p.relative_to(ROOT)):sha(p) for p in (envpath,giftpath)},
       'output_sha256':{p.name:sha(p) for p in (png,svg)},
       'plotted_rows':plotted,'verification':{'primary_rows':3,'sensitivity_rows':5,'modifier_rows':9,
       'exactly_one_source_row_per_mark':True,'all_bounds_finite_and_contain_estimate':True,
       'family9_encloses_nominal':True},
       'notes':['A and B read REML hypermeans, not equal-species arithmetic means.',
                'B shows different retained species panels, not paired changes.',
                'C reads difference_first_minus_reference and existing genus-bootstrap interval columns.',
                'Family9 intervals are approximate Bonferroni-familywise intervals, not nine simultaneous exact intervals.',
                'No model refitting or fabricated estimates; all plotted statistics read from the two source CSVs.']}
    (args.outdir/'figure_provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
    print(png);print(svg);print('SOURCE_ROWS_AND_BOUNDS_VERIFIED')

if __name__=='__main__':main()
