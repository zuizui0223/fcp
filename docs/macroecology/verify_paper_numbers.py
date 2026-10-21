"""Check manuscript evidence against frozen original analytical CSVs.

This script does NOT estimate independent biological effects; see PR #193/194/
#191 source-exact GitHub Actions for the archived numerical computations.
"""
from __future__ import annotations
import argparse,json,hashlib
from datetime import date
from pathlib import Path
import pandas as pd

def one(data, **kw):
    s=data
    for k,v in kw.items():s=s.loc[s[k]==v]
    if len(s)!=1:raise AssertionError((kw,len(s)))
    return s.iloc[0]

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--macro-ledger',type=Path,required=True)
    p.add_argument('--matched-ledger',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--photo-source',type=Path)
    p.add_argument('--abiotic-source',type=Path)
    a=p.parse_args()
    m=pd.read_csv(a.macro_ledger)
    t=pd.read_csv(a.matched_ledger)
    assert m.shape==(1280,17) and t.shape==(576,11)
    checks=[]
    def compare(k,obs,expected,tol=.004):
        assert abs(obs-expected)<=tol,(k,obs,expected)
        checks.append({'key':k,'observed':round(float(obs),9),'manuscript_rounding':expected,'abs_delta':round(abs(float(obs)-expected),9),'status':'PASS'})
    for anchor,point,lo,hi in [
       ('wc_bio12',2.255,.319,4.190),
       ('wc_vapr_annual_mean',3.952,1.199,6.704),
       ('wc_bio5',1.257,-.648,3.161),
       ('wc_srad_annual_mean',-.533,-2.944,1.878)]:
        v=one(m,mode='climate',quality='all',region='world',adjustment='genus_geo',weight='genus_equal',anchor=anchor,hue='white')
        compare('world '+anchor+' point',v.direction_pp_per_env_global_sd,point)
        compare('world '+anchor+' maxT low',v.maxT_ci_lo_pp,lo)
        compare('world '+anchor+' maxT high',v.maxT_ci_hi_pp,hi)
    for reg,qual,adj,point,lo,hi in [
       ('world','strict','genus_geo',1.229,-2.046,4.504),
       ('world','all','genus_geocell',.808,-3.556,5.172),
       ('E90_E180','all','genus_geo',4.716,1.565,7.866),
       ('E90_E180','all','genus_geocell',-2.076,-8.794,4.643)]:
        v=one(m,mode='climate',quality=qual,region=reg,adjustment=adj,weight='genus_equal',anchor='wc_bio12',hue='white')
        compare(reg+' '+qual+' '+adj+' point',v.direction_pp_per_env_global_sd,point)
        compare(reg+' '+qual+' '+adj+' ci low',v.maxT_ci_lo_pp,lo)
        compare(reg+' '+qual+' '+adj+' ci hi',v.maxT_ci_hi_pp,hi)
    for reg,point,lo,hi in [('world',.518,-1.253,2.289),('E90_E180',3.178,-.258,6.613)]:
        v=one(t,region=reg,weight='genus_equal',type='paired_difference',comparison='spatial_adjustment_on_same_species',anchor='wc_bio12',hue='white')
        compare(reg+' B minus C point',v.point_pp,point)
        compare(reg+' B minus C maxT low',v.maxT95_lo_pp,lo)
        compare(reg+' B minus C maxT high',v.maxT95_hi_pp,hi)
    hashes={}
    for key,path,expected in [
        ('photo_source_sha256',a.photo_source,'2d6097ba1b0434632e13bbe79dd78180608a287ff5dabebc770b5c0a2cfb52a0'),
        ('environment_source_sha256',a.abiotic_source,'6e287294b3a05eaa383be3de644551825915bfecca9c4bde6950b4890f5e16ca')]:
        if path is not None:
            digest=hashlib.file_digest(path.open('rb'),'sha256').hexdigest()
            assert digest==expected,(key,digest)
            hashes[key]=digest
    out={'schema':'fcp_macro_paper_numbers_reconciliation_v0_1','date':'2026-10-14',
         'status':'PASS','n_checked_rounded_manuscript_values':len(checks),
         'macro_source_rows':len(m),'matched_source_rows':len(t),
         'source_hashes':hashes,'checks':checks,
         'note':'Manuscript rounding/extraction audit only. Original statistical code and source replay are under PR #193/#194 and hosted Actions #38502660130; independent-source biological confirmation is not claimed.'}
    a.output.write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({k:out[k] for k in ['status','n_checked_rounded_manuscript_values','macro_source_rows','matched_source_rows','source_hashes']},indent=2))

if __name__=='__main__':main()
