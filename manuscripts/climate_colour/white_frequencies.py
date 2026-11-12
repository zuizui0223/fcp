"""Describe observed white-label frequencies. Usage: python white_frequencies.py MEASURED_CSV"""
from pathlib import Path
import hashlib
import sys
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist

root = Path(__file__).resolve().parent
source = Path(sys.argv[1])
assert hashlib.sha256(source.read_bytes()).hexdigest() == "8e7f10f0427927208485b2da6fa9c46bb681d6f3611e53b8dfdc606cfc8dabe6"
ledger = pd.read_csv(root.parents[1]/"results/fcp_np_all4475_unified_species_ledger_20261021.csv").set_index("inat_taxon_id")
counts = ledger.filter(regex="^miss_coarse_photo_")
repeated = counts.ge(5).sum(axis=1).ge(2)
assert len(ledger)==4475 and repeated.sum()==2400
s = pd.DataFrame({"repeated":repeated, "n":counts.sum(axis=1), "white":counts.miss_coarse_photo_white})
cols = "inat_taxon_id photo_id observer_id observed_on latitude longitude positional_accuracy_m global_classifiable morph".split()
photos = pd.read_csv(source,usecols=cols,low_memory=False)
photos = photos[photos.global_classifiable & photos.inat_taxon_id.isin(ledger.index)].copy()
assert len(photos)==int(s.n.sum()) and photos.photo_id.is_unique
photos=photos.sort_values("photo_id")
photos["white"]=photos.morph.eq("white")
assert photos.groupby("inat_taxon_id").white.sum().reindex(s.index,fill_value=0).eq(s.white).all()
one=photos.dropna(subset=["observer_id"]).drop_duplicates(["inat_taxon_id","observer_id"]).groupby("inat_taxon_id").white.agg(["size","sum"])
s["observer_n"],s["observer_white"]=one["size"],one["sum"]
# Choose a single photo-centred neighbourhood per species. Selection never reads colour.
# Distinct observers, same year-month, max 2.5 km from the centre (max 5 km pairwise).
q=photos.dropna(subset=["observer_id","latitude","longitude"]).copy()
q=q[q.positional_accuracy_m.between(0,5000)]
q["ym"]=pd.to_datetime(q.observed_on,format="mixed",utc=True,errors="coerce").dt.strftime("%Y-%m")
records={}
for (tax,ym),g in q.groupby(["inat_taxon_id","ym"]):
    if len(g)<3 or g.observer_id.nunique()<3: continue
    lat,lon=np.deg2rad(g[["latitude","longitude"]].to_numpy()).T
    xyz=np.column_stack((np.cos(lat)*np.cos(lon),np.cos(lat)*np.sin(lon),np.sin(lat)))
    near=cdist(xyz,xyz)<=2*np.sin(2.5/(2*6371.0088))
    for centre in np.flatnonzero(near.sum(axis=1)>=3):
        window=g.iloc[near[centre]].drop_duplicates("observer_id")
        n=len(window)
        if n<3:continue
        key=(n,-int(g.photo_id.iloc[centre]))
        if tax not in records or key>records[tax][0]:
            records[tax]=(key,ym,int(window.white.sum()))
for tax,((n,minusid),ym,white) in records.items():
    s.loc[tax,["local_n","local_white","centre_photo","local_month"]]=[n,white,-minusid,ym]
s.to_csv(root/"Table_S6_white_frequencies.csv.gz",compression={"method":"gzip","mtime":0})
for name,ss,n,w in [("all",s,"n","white"),("repeated",s[s.repeated],"n","white"),("observer",s[s.repeated],"observer_n","observer_white")]+[
    ("local"+str(k),s[s.repeated & s.local_n.ge(k)],"local_n","local_white") for k in [3,5,8]]:
    ss=ss[ss[n]>0]
    frac=ss[w]/ss[n]
    print(name,len(ss),int(ss[w].sum()),int(ss[n].sum()),ss[w].sum()/ss[n].sum(),frac.quantile([.25,.5,.75]).to_list(),
          frac.between(.4,.6).sum(),((frac==0)|(frac==1)).sum(),ss[n].median())
