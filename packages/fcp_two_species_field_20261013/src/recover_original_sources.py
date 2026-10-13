#!/usr/bin/env python3
"""Recover public-author study originals via DOI-matched Zenodo mirror with hash gates.

Source acquisition, not a biological analysis. Never uses Dryad challenge endpoints.
"""
from __future__ import annotations
import argparse,hashlib,json,urllib.request,zipfile,io
from pathlib import Path

RECORDS={
  5003124:{'doi':'10.5061/dryad.j8q34','needed':{
        'data.zip':('Sisyrinchium_data.zip','008ec196aa5acb769dbfd36e835e32ba',127502)}},
  4993384:{'doi':'10.5061/dryad.q0032','needed':{
        'natural_populations_individual.txt':('Boechera_natural_populations_individual.txt','85990c8a30a8b9ef3a2c045bb6ffb72a',22333),
        'natural_populations_plot.txt':('Boechera_natural_populations_plot.txt','1c284e1543450aaacf161cb7400d0a5a',2224),
        'README_for_natural_populations_individual.txt':('Boechera_README_individual.txt',None,900),
        'README_for_natural_populations_plot.txt':('Boechera_README_plot.txt',None,872)
  }}
}

def remote_get(url:str,limit=2000000):
    if not url.startswith('https://zenodo.org/'):
        raise ValueError('UNAPPROVED_HOST')
    req=urllib.request.Request(url,headers={'User-Agent':'FCP-open-field-data-source-audit/2026-10-13'})
    with urllib.request.urlopen(req,timeout=40) as h:
        content=h.read(limit+1)
    if len(content)>limit:raise ValueError('RAW_ARCHIVE_TOO_LARGE')
    return content

def audit_path(path:Path,expect_md5:str|None,expected_bytes:int):
    data=path.read_bytes()
    assert len(data)==expected_bytes,(str(path),len(data),expected_bytes)
    md5=hashlib.md5(data).hexdigest()
    if expect_md5 is not None: assert md5==expect_md5,(str(path),md5,expect_md5)
    return {'bytes':len(data),'md5':md5,'sha256':hashlib.sha256(data).hexdigest()}

def retrieve(outdir:Path,reuse_existing=False):
    outdir.mkdir(parents=True,exist_ok=True)
    records=[]
    for record_id,desc in RECORDS.items():
        entry={'record_id':record_id,'original_dryad_doi':desc['doi'],'files':[]}
        records.append(entry)
        if reuse_existing:
            metadata=None
        else:
            metadata=json.loads(remote_get(f'https://zenodo.org/api/records/{record_id}'))
            assert metadata.get('id')==record_id and metadata.get('doi')==desc['doi'],'SOURCE_RECORD_DOI_MISMATCH'
            hosted={f['key']:f for f in metadata.get('files',[])}
        for hosted_name,(local_name,md5_expect,len_expect) in desc['needed'].items():
            out=outdir/local_name
            item={'source_name':hosted_name,'local_name':local_name,'status':'UNCONFIRMED'}
            entry['files'].append(item)
            if not reuse_existing:
                if hosted_name not in hosted: raise RuntimeError('MISSING_EXPECTED_FILE '+hosted_name)
                h=hosted[hosted_name]
                assert h['size']==len_expect
                url=h['links']['self']
                assert url==f'https://zenodo.org/api/records/{record_id}/files/{hosted_name}/content'
                reported=h['checksum'].split(':',1)[-1]
                if md5_expect is not None: assert reported==md5_expect
                data=remote_get(url)
                if len(data)!=len_expect or hashlib.md5(data).hexdigest()!=reported:
                    raise RuntimeError('PUBLISHER_BYTES_HASH_MISMATCH '+hosted_name)
                out.write_bytes(data)
                item['source_url']=url
                item['publisher_md5']=reported
            item.update(audit_path(out,md5_expect,len_expect))
            item['status']='SOURCE_BYTES_VERIFIED'
            if local_name.endswith('.zip'):
                with zipfile.ZipFile(out) as z:
                    assert z.testzip() is None,'CORRUPT_SOURCE_ZIP'
                    item['verified_frequency_txt']='data/frequency.txt' in z.namelist()
                    item['apple_double_skipped']=sum(name.startswith('__MACOSX/') for name in z.namelist())
                    assert item['verified_frequency_txt']
    receipt={'schema':'fcp_two_original_dryad_doi_zenodo_mirror_recovery_v1','status':'SOURCE_BYTES_VERIFIED','runs_without_causal_inference':True,'sources':records}
    (outdir/'source_receipt.json').write_text(json.dumps(receipt,indent=2,ensure_ascii=False)+'\n')
    return receipt

def main():
    p=argparse.ArgumentParser();p.add_argument('--source-dir',type=Path,default=Path(__file__).resolve().parents[1]/'source');p.add_argument('--reuse-existing',action='store_true')
    args=p.parse_args()
    result=retrieve(args.source_dir,args.reuse_existing)
    print(json.dumps({'status':result['status'],'files':sum(len(z['files']) for z in result['sources'])},ensure_ascii=False))

if __name__=='__main__':main()
