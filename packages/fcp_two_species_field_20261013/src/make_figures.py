#!/usr/bin/env python3
"""Static plots from frozen observed site-year outcome ledger, no model refitting."""
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
docs=ROOT/'docs';docs.mkdir(parents=True,exist_ok=True)
d=pd.read_csv(ROOT/'results/siteyear_metrics.csv')
for sp,display in [('sisyrinchium','Sisyrinchium M-species'),('boechera','Boechera stricta')]:
    x=d[(d.species==sp)&(d.status=='ELIGIBLE')].copy()
    x=x.sort_values('delta',ascending=True)
    labels=[i.replace('single_census_unreported_year|','').replace('Y11|','2011 · ').replace('Y13|','2013 · ').replace('_',' ') for i in x.siteyear]
    fig,ax=plt.subplots(figsize=(9,0.43*len(x)+2))
    ax.barh(labels,100*x.delta.to_numpy())
    ax.axvline(0,linestyle='--',linewidth=.8)
    ax.set_xlabel('Excess between-quadrat colour diversity vs fixed-count null (percentage points)')
    ax.set_title(display + ' · within-field-site morph allocation')
    ax.grid(axis='x',alpha=.2)
    ax.set_axisbelow(True)
    fig.tight_layout()
    path=docs/(sp+'_quadrat_excess_by_site.png');fig.savefig(path,dpi=220)
    plt.close(fig)

x=pd.read_csv(ROOT/'results/posthoc_leave_one_siteyear_out.csv')
f=x[x.species=='boechera'].copy().sort_values('after_leave_one_out_delta_percentage_points')
fig,ax=plt.subplots(figsize=(9,5))
labels=[s.split('|')[-1].replace('_',' ') for s in f.excluded_siteyear]
ax.barh(labels,f.after_leave_one_out_delta_percentage_points)
ax.axvline(f.full_mean_delta_percentage_points.iloc[0],linestyle='--',label='All six sites',linewidth=.9)
ax.axvline(0,linewidth=.7)
ax.set_xlabel('Mean excess after removing named site (percentage points)')
ax.set_ylabel('Removed site')
ax.set_title('Boechera stricta · posthoc leave-one-site-out fragility')
ax.legend(loc='lower right')
ax.grid(axis='x',alpha=.2)
ax.set_axisbelow(True);fig.tight_layout()
fig.savefig(docs/'boechera_leave_one_site_out.png',dpi=220)
plt.close(fig)
print('PLOTS_CREATED', [p.name for p in sorted(docs.glob('*.png'))])
