#!/usr/bin/env python3
"""Optional illustrative plots of the audited observational outputs."""
from pathlib import Path
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
plt.rcParams['font.family']='Noto Sans CJK JP'
plt.rcParams['font.size']=10
r=json.loads((ROOT/'results/result.json').read_text(encoding='utf8'))
rows=r['primary']['annual_side_contrasts']
years=[x['year'] for x in rows]
west=[100*x['west']['equal_transect_white_fraction'] for x in rows]
east=[100*x['east']['equal_transect_white_fraction'] for x in rows]
fig,ax=plt.subplots(figsize=(9,4.6))
ax.plot(years,west,marker='o',linewidth=1.75,label='西側（南北2トランセクト等重み）')
ax.plot(years,east,marker='s',linewidth=1.75,label='東側（南北2トランセクト等重み）')
ax.set_xlabel('調査年');ax.set_ylabel('白花頻度 (%)');ax.set_ylim(-3,103)
ax.set_xticks(list(range(1992,1999)));ax.grid(axis='y',alpha=0.23)
for yr,lbl in [(1994,'1994 / 保留: 9個体'),(1997,'1997 / 植物0')]:
 ax.axvline(yr,linewidth=0.8,linestyle=':',color='gray',alpha=0.7)
 ax.text(yr,68,lbl,ha='center',va='center',fontsize=9,color='gray')
ax.legend(loc='center right',frameon=False)
ax.set_title('同じ固定区の花色型頻度は年ごとにどれほど変わるか')
fig.tight_layout()
fig.savefig(ROOT/'docs/linanthus_observed_annual_ravine_sides.png',dpi=170)
plt.close(fig)

rows=r['primary']['heldout_years']
fig,ax=plt.subplots(figsize=(9,4.6))
for name,key,marker in [('全域平均のみ','mse_global','o'),('谷の東西','mse_side','s'),('4区別','mse_section','^')]:
 ax.plot([r['year'] for r in rows],[r[key] for r in rows],label=name,marker=marker,linewidth=1.6)
ax.set_yscale('log');ax.set_xticks(years);ax.set_xlabel('評価年（その年は学習から全除外）')
ax.set_ylabel('調査区間の白花頻度の二乗予測誤差（対数軸）')
ax.set_title('新しい調査年の花色頻度：全域平均 vs 地理情報')
ax.legend(frameon=False,loc='center right')
ax.grid(axis='y',alpha=0.23)
fig.tight_layout();fig.savefig(ROOT/'docs/linanthus_leftout_year_prediction.png',dpi=170)
plt.close(fig)
print('FIGURE_OUTPUT_COMPLETE')
