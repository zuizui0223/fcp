#!/usr/bin/env python3
"""Read the actual downloaded GitHub Actions ZIP and independently compare outputs.

No network access, no new ecological hypothesis. Source bytes exact; numerical results
may vary at <=1e-8 with pinned floating-point libraries across hardware/runtime.
"""
import argparse,csv,hashlib,io,json,math,zipfile
from pathlib import Path

SHARE='fcp-sisyrinchium-shore-gradient-transport-20261013'
ART_SHA='2ee07238f26a7cc1f92cb719f074d5819d21f081b7549cd82c0ac802530175c6'


def digest(b):return hashlib.sha256(b).hexdigest()

def nested(a,b,checks,label,maxdiff):
    checks[0]+=1
    if isinstance(a,dict):
        assert isinstance(b,dict) and set(a)==set(b),label
        for k in a:nested(a[k],b[k],checks,label+'.'+str(k),maxdiff)
    elif isinstance(a,list):
        assert isinstance(b,list) and len(a)==len(b),label
        for i,(x,y) in enumerate(zip(a,b)):nested(x,y,checks,f'{label}[{i}]',maxdiff)
    elif type(a) in (int,float) and type(b) in (int,float):
        assert math.isfinite(a) and math.isfinite(b)
        delta=abs(a-b);maxdiff[0]=max(maxdiff[0],delta)
        assert math.isclose(a,b,abs_tol=1e-8,rel_tol=1e-8),(label,a,b)
    else:assert a==b,(label,a,b)


def audit(root,artifact,receipt_path):
    root=Path(root);artifact=Path(artifact)
    b=artifact.read_bytes();assert digest(b)==ART_SHA
    checks=[0];maxdiff=[0.];files=[]
    with zipfile.ZipFile(io.BytesIO(b)) as z:
        assert z.testzip() is None
        names=set(z.namelist())
        for p in ['Sisyrinchium_data.zip','parent_result.json','parent_predictions.csv']:
            key='source/'+p
            assert key in names
            a=(root/key).read_bytes();remote=z.read(key)
            assert a==remote,(key,digest(a),digest(remote))
            files.append({'file':key,'bytes':len(remote),'source_byte_identical':True,'sha256':digest(remote)})
        for fn in ['result.json','posthoc_distance_scale.json','independent_verification.json']:
            key='results/'+fn
            original=json.loads((root/key).read_text());rem=json.loads(z.read(key).decode())
            nested(original,rem,checks,key,maxdiff)
            files.append({'file':key,'exact':(root/key).read_bytes()==z.read(key),
                          'numerically_verified':True,'source_sha256':digest((root/key).read_bytes()),'remote_sha256':digest(z.read(key))})
        for fn in ['named_site_transfer_scores.csv','borrowed_line_predictions.csv','outer_site_and_inner_line_fit_ledger.csv','site_deletion_transfer_sensitivity.csv','posthoc_distance_scale_predictions.csv']:
            key='results/'+fn
            a=list(csv.DictReader((root/key).open(newline='')))
            remote=z.read(key).decode('utf-8-sig')
            r=list(csv.DictReader(io.StringIO(remote)))
            assert len(a)==len(r)
            for j,(x,y) in enumerate(zip(a,r)):
                assert set(x)==set(y)
                for k in x:
                    checks[0]+=1
                    if x[k]==y[k]:continue
                    try:v=float(x[k]);w=float(y[k])
                    except ValueError:raise AssertionError(('NONNUMERIC',fn,j,k,x[k],y[k]))
                    maxdiff[0]=max(maxdiff[0],abs(v-w))
                    assert math.isclose(v,w,rel_tol=1e-8,abs_tol=1e-8),(fn,j,k,v,w)
            files.append({'file':key,'exact':(root/key).read_bytes()==z.read(key),'numerically_verified':True,
                          'source_sha256':digest((root/key).read_bytes()),'remote_sha256':digest(z.read(key))})
        for fn in ['full_site_distance_slopes.png','heldout_cross_site_gradient_transport.png']:
            key='docs/'+fn
            assert (root/key).read_bytes()==z.read(key),(key,'FIGURE_BYTES_CHANGED')
            files.append({'file':key,'exact':True,'sha256':digest(z.read(key))})
    rec={'schema':'fcp_sisyrinchium_shoreline_gradient_transport_remote_execution_audit_v1',
         'date':'2026-10-13','scientific_commit':'0656a88c6ee13eed4d8d71056d431d39746f2429',
         'github_pr':186,'github_push_run':38367758188,'github_pr_run':38367980019,
         'github_push_job':115196889358,'status':'PASS',
         'github_actions_source_hashes_exact':True,'remote_artifact_sha256':digest(b),
         'remote_artifact_size_bytes':len(b),'remote_artifact_id':11797493503,
         'focused_tests_passed':16,'separate_independent_newton_checks_passed':1112,
         'local_vs_remote_recursive_or_csv_checks':checks[0],
         'local_vs_remote_max_absolute_numeric_difference':maxdiff[0],
         'source_original_2016_field_study_only':True,'all_six_sites_same_original_2016_study':True,
         'posthoc_standardized_diagnostic_not_primary':True,
         'not_exact_byte_output_same':any(x.get('exact') is False for x in files),
         'whole_repo_pytest_not_run':True,
         'files':files,
         'hard_limits':['No new or independent field sampling','Six selected named sites are not random evolutionary replicates',
                        'No causal attribution to water, congener interference, selection or genetic differentiation',
                        'Transferred slope excludes target site, but target intercept trained from target other lines',
                        'Floating-point final results numerically equivalent but not byte-exact across execution environments']}
    Path(receipt_path).write_text(json.dumps(rec,ensure_ascii=False,sort_keys=True,indent=2)+'\n')
    print(json.dumps({'status':'PASS','files':len(files),'checks':checks[0],'maxdiff':maxdiff[0],'artifact_sha256':digest(b)},ensure_ascii=False))
    return rec

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',required=True);p.add_argument('--artifact',required=True);p.add_argument('--out',required=True)
    a=p.parse_args();audit(a.root,a.artifact,a.out)
