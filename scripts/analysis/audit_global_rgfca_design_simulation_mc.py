#!/usr/bin/env python3
"""Read-only MC uncertainty audit of the completed S1-S3 artifact, not a rerun."""
from __future__ import annotations
import argparse, csv, hashlib, io, json, math, zipfile
from pathlib import Path

EXPECTED_RESULT_BLOB = '87962d88f66f86bae534604d15e237a91b98cf6e'
RESULT = 'global_rgfca_design_power_identifiability_simulation_result_v1.json'
FILES = [RESULT, 'global_rgfca_design_power_s1_v1.csv', 'global_rgfca_outer_stability_s2_v1.csv', 'global_rgfca_environmental_holm_power_s3_v1.csv']

def blob(b: bytes) -> str:
    return hashlib.sha1(b'blob '+str(len(b)).encode()+b'\0'+b).hexdigest()

def rate(p: float, n: int) -> dict:
    k = round(p*n)
    if n <= 0 or not 0 <= p <= 1 or abs(k-p*n) > 1e-6:
        raise ValueError('probability does not correspond to a complete event count')
    z=1.959963984540054; den=1+z*z/n
    center=(p+z*z/(2*n))/den
    half=z*math.sqrt(p*(1-p)/n+z*z/(4*n*n))/den
    return {'events':k,'replicates':n,'rate':p,'mc_se':math.sqrt(p*(1-p)/n),
            'mc_wilson95_low':max(0.,center-half),'mc_wilson95_high':min(1.,center+half)}

def audit(archive: Path, output: Path) -> dict:
    with zipfile.ZipFile(archive) as z:
        content={name:z.read(name) for name in FILES}
    if blob(content[RESULT]) != EXPECTED_RESULT_BLOB:
        raise ValueError('artifact result differs from verified repository result blob')
    r=json.loads(content[RESULT]); s1,s2,s3=[list(csv.DictReader(io.StringIO(content[name].decode()))) for name in FILES[1:]]
    if (len(s1),len(s2),len(s3)) != (38,5,105): raise ValueError('incomplete scenario tables')
    calibration=[{'frame':x['frame'],**rate(x['null']['evaluation_type1_fraction'],x['null']['evaluation_replicates'])} for x in r['S1']['frames']]
    audited_s1=[{'frame':x['frame'],'amplitude':float(x['signal_amplitude']),'shared_fraction':float(x['shared_species_fraction']),**rate(float(x['power_ge_threshold']),int(x['simulation_replicates']))} for x in s1]
    stability=[]
    for x,ref in zip(s2,r['S2']['rows']):
        for key in x:
            if not math.isclose(float(x[key]),float(ref[key]),rel_tol=1e-12,abs_tol=1e-14):
                raise ValueError('S2 JSON/CSV disagreement: '+key)
        stability.append({'outer_count':int(x['outer_count']),
                          'pair_median_r':float(x['independent_pair_median_r']),
                          'to_frozen_median_r':float(x['to_frozen_200_median_r']),
                          'pair_ge_09':rate(float(x['independent_pair_prob_r_ge_0_9']),int(x['bootstrap_replicates'])),
                          'to_frozen_ge_09':rate(float(x['to_frozen_200_prob_r_ge_0_9']),int(x['bootstrap_replicates']))})
    mde=[]; all_s3=[]
    for target in sorted({x['target_block'] for x in s3}):
        for corr in sorted({float(x['correlation']) for x in s3}):
            rows=[x for x in s3 if x['target_block']==target and float(x['correlation'])==corr]
            if len(rows)!=7: raise ValueError('incomplete block/correlation effect grid')
            above=[float(x['true_target_effect']) for x in rows if float(x['target_holm_power'])>=.8]
            mde.append({'target':target,'correlation':corr,'first_grid_effect_power_ge_08':min(above) if above else None})
            for x in rows:
                all_s3.append({'target':target,'correlation':corr,'effect':float(x['true_target_effect']),
                               'target_holm':rate(float(x['target_holm_power']),int(x['replicates'])),
                               'family_any':rate(float(x['family_any_discovery_probability']),int(x['replicates']))})
    out={'status':'complete_read_only_mc_uncertainty_audit','artifact_sha256':hashlib.sha256(archive.read_bytes()).hexdigest(),
         'source_git_blobs':{key:blob(value) for key,value in content.items()},
         'calibration':calibration,'s1':audited_s1,'s2':stability,'s3_mde_grid':mde,'s3':all_s3,
         'observed_decisions_changed':False,
         'limits':['Intervals are Monte Carlo uncertainty conditional on the stored calibration/generator, not biological confidence intervals.',
                   'S1 absolute latent differences do not validate tied palette JSD.',
                   'Species sharing cannot be inferred from S1 rejection against a spatially unstructured null.',
                   'The two S1 geometry frames change species, locations and schedules; their contrast is not an isolated causal missingness effect.',
                   'S2 bootstrap maps resample the same fixed balanced outer maps, not independent biological data.',
                   'S3 is a Gaussian shifted-null approximation, not evidence that the observed effects are true.']}
    if output.exists(): raise FileExistsError('refusing to overwrite audit')
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    return out

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--archive',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();r=audit(a.archive,a.output)
    print(json.dumps({'status':r['status'],'calibration':r['calibration'],'s3_mde_grid':r['s3_mde_grid']}))
