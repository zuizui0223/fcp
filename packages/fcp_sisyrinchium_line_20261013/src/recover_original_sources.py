#!/usr/bin/env python3
"""Strict source recovery from DOI-preserving public Zenodo original mirror."""
from pathlib import Path
import urllib.request,hashlib,zipfile,io,json,argparse
URL='https://zenodo.org/api/records/5003124/files/data.zip/content'
ZIP_MD5='008ec196aa5acb769dbfd36e835e32ba'
ZIP_LEN=127502
FREQ_SHA='a12cdf46483a6d8b78bd00110c196a25f3c8849967ae8a2c37d6aee3d29f516e'

def recover(outdir):
    out=Path(outdir);out.mkdir(parents=True,exist_ok=True)
    urls=[URL,'https://zenodo.org/records/5003124/files/data.zip?download=1']
    attempts=[]
    raw=None
    successful_url=None
    for candidate in urls:
        try:
            req=urllib.request.Request(candidate,headers={
                'User-Agent':'FCP-original-source-audit/2026-10-13 (noncommercial ecological science)',
                'Accept':'application/json, text/plain, application/octet-stream, */*'})
            with urllib.request.urlopen(req,timeout=45) as conn:
                data=conn.read(200000)
                if conn.status!=200:raise ValueError('NOT_HTTP_200')
            md5=hashlib.md5(data).hexdigest()
            if len(data)!=ZIP_LEN or md5!=ZIP_MD5:
                raise ValueError(f'SOURCE_ZIP_AUTHENTICATION_FAILED:{len(data)}:{md5}')
            raw=data
            successful_url=candidate
            break
        except Exception as ex:
            attempts.append({'url':candidate,'error_type':type(ex).__name__,'reason':str(ex)[:180]})
    if raw is None:
        raise RuntimeError('ORIGINAL_PUBLIC_ZENODO_SOURCE_UNAVAILABLE: '+json.dumps(attempts))
    md5=hashlib.md5(raw).hexdigest()
    with zipfile.ZipFile(io.BytesIO(raw)) as z:
        if z.testzip() is not None:raise ValueError('SOURCE_ZIP_CRC_FAILED')
        freq=z.read('data/frequency.txt')
    sha=hashlib.sha256(freq).hexdigest()
    if sha!=FREQ_SHA:raise ValueError('SOURCE_FREQUENCY_MEMBER_SHA256_FAILED')
    (out/'Sisyrinchium_data.zip').write_bytes(raw)
    (out/'frequency.txt').write_bytes(freq)
    rec={'status':'SOURCE_BYTES_VERIFIED','original_doi':'10.5061/dryad.j8q34',
          'zenodo_public_record':5003124,'original_url':successful_url,'attempts':attempts,'original_zip_bytes':len(raw),
          'original_zip_md5':md5,'original_zip_sha256':hashlib.sha256(raw).hexdigest(),
          'original_frequency_member':'data/frequency.txt','frequency_sha256':sha,
          'not_independent_new_field_plant_observations':True,
          'other_species_l_s_not_white_M_morph':True}
    (out/'source_receipt.json').write_text(json.dumps(rec,indent=2,sort_keys=True)+'\n')
    print(json.dumps(rec))
    return rec

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--outdir',default='source');a=p.parse_args()
    recover(a.outdir)
