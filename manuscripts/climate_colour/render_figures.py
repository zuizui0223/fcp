"""Rebuild the manuscript figures from the accompanying photo and species tables."""
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image

root = Path(__file__).resolve().parent
hues = 'white yellow orange red pink magenta purple blue bronze'.split()
colours = ['#f8f8f5','#eac52f','#ef8c2e','#b32f32','#ef9fbc','#ac418e','#7355a4','#4a78b0','#94653e']
thermal, rain = '#a44730', '#267479'
plt.rcParams.update({'font.size':10, 'axes.spines.right':False, 'axes.spines.top':False})
photos = pd.read_csv(root/'Table_S1_photographs.csv')
fig, axes = plt.subplots(2,2,figsize=(7.6,7.6),layout='constrained')
for ax, row, letter, note in zip(axes.flat, photos.itertuples(), 'ABCD',
        ['Purple and yellow flowers','Yellow flower','Pale petals and yellow centre','Flip instability']):
    ax.axis('off')
    ax.imshow(Image.open(root/f'photos/{row.photo_id}.jpg'))
    ax.set_title(f'{letter}   {row.species}',loc='left',fontstyle='italic',fontsize=11)
    ax.set_xlabel('')
    ax.text(0,-.035,note,va='top',transform=ax.transAxes)
    bar=ax.inset_axes([0,-.20,1,.075])
    bar.set_xlim(0,1)
    bar.set_ylim(0,1)
    bar.axis('off')
    if row.global_classifiable:
        left=0
        for h,c in zip(hues,colours):
            width=getattr(row,'flower_fraction_'+h)
            bar.barh(.5,width,left=left,height=.9,color=c,edgecolor='#555',linewidth=.35)
            left+=width
        assert abs(left-1)<1e-5
    else:
        bar.text(0,.5,'Excluded from colour analysis',va='center',fontweight='bold')
from matplotlib.patches import Patch
fig.legend(handles=[Patch(facecolor=c,edgecolor='#999',label=h) for h,c in zip(hues,colours)],
           loc='lower center',ncol=5,frameon=False,fontsize=9,bbox_to_anchor=(.5,.003))
fig.get_layout_engine().set(rect=(0,.10,1,.90),h_pad=.20)
fig.savefig(root/'Figure_1_measurement.png',dpi=225,bbox_inches='tight')
plt.close(fig)

s=pd.read_csv(root/'Table_S2_species_coefficients.csv.gz')
domains=['thermal_level','precip_amount']
cols=[[f'{d}_{h}' for h in hues] for d in domains]
paired=s.dropna(subset=cols[0]+cols[1])
v=np.stack([paired[c].mean().to_numpy() for c in cols])
q=np.array([8]+[-1]*8)/np.sqrt(72)
align=(v@q)**2/(v*v).sum(axis=1)
assert len(paired)==2689 and np.allclose(align,[.837,.243],atol=.00055)
fig,ax=plt.subplots(3,1,figsize=(8.1,8.15),layout='constrained',height_ratios=[1.55,.70,1.72])
x=np.arange(9)
for i,(name,color) in enumerate([('Temperature',thermal),('Precipitation',rain)]):
    ax[0].bar(x+(i-.5)*.38,v[i]*100,.37,color=color,label=name)
ax[0].axhline(0,color='#444',lw=.7)
ax[0].set_xticks(x,hues,rotation=25)
ax[0].set_ylabel('Mean partial coefficient\n(percentage points per index unit)')
ax[0].legend(frameon=False)
ax[0].set_title('A   Mean colour change on the same 2,689 species',loc='left')
ax[1].barh([1,0],align*100,color=[thermal,rain],height=.55)
for y,val in zip([1,0],align*100):
    ax[1].text(val+1,y,f'{val:.1f}%',va='center')
ax[1].set_xlim(0,100)
ax[1].set_yticks([1,0],['Temperature','Precipitation'])
ax[1].set_xlabel('Squared alignment with the white axis (%)')
ax[1].set_title('B   Orientation of each mean vector',loc='left')
contrasts = pd.read_csv(root/'Table_S3_cohort_contrasts.csv')
y = np.arange(len(contrasts))
delta = contrasts.difference_pp.to_numpy()
ax[2].errorbar(delta, y, xerr=[delta-contrasts.lower95_pp, contrasts.upper95_pp-delta],
               fmt='o', color='#275465', capsize=3, markersize=5)
ax[2].axvline(0, color='#555', lw=.8, ls='--')
ax[2].set_yticks(y, contrasts.panel)
ax[2].invert_yaxis()
ax[2].set_xlim(-40, 100)
ax[2].set_xlabel('Temperature minus precipitation alignment (percentage points, 95% interval)')
ax[2].set_title('C   Difference in alignment after averaging species vectors', loc='left')
fig.savefig(root/'Figure_2_climate.png',dpi=225)
plt.close(fig)

fig,ax=plt.subplots(1,2,figsize=(8,4),layout='constrained')
b=s.thermal_level_white.dropna()*100
ax[0].hist(b,bins=np.linspace(-55,55,48),color='#798f99')
ax[0].axvline(0,color='#333',ls=':',label='Zero')
ax[0].axvline(b.mean(),color=thermal,label=f'Mean {b.mean():+.2f}')
ax[0].set_xlabel('White coefficient (percentage points)')
ax[0].set_ylabel('Species')
ax[0].set_title('A   All 2,837 thermal estimates',loc='left')
ax[0].legend(frameon=False)
f=s[['forward','reverse']].dropna()*100
assert len(f)==1919 and round(b.gt(0).mean()*100,1)==54.2
ax[1].scatter(f.forward,f.reverse,s=6,c=rain,alpha=.35,edgecolors='none',rasterized=True)
ax[1].axhline(0,color='#555',lw=.6)
ax[1].axvline(0,color='#555',lw=.6)
lim=np.ceil(f.abs().max().max()/10)*10
ax[1].set(xlim=(-lim,lim),ylim=(-lim,lim),aspect='equal',xlabel='Forward fold (percentage points)',ylabel='Reverse fold (percentage points)')
ax[1].set_title('B   Paired observer estimates',loc='left')
ax[1].text(.04,.97,'1,919 species\nr = 0.018',va='top',transform=ax[1].transAxes)
fig.savefig(root/'Figure_3_species.png',dpi=225)
plt.close(fig)
