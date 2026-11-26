#!/usr/bin/env python3
"""Within-quadrat white-versus-purple fruit-set contrasts, original Takahashi et al. (2016).
The data are observational; this is not a manipulation of morph frequency.
Run: python scripts/analysis/fcp_sisyrinchium_paired_fitness_20261126.py ZIP OUTPUT_DIR
"""
import csv, hashlib, io, json, math, sys, zipfile
from collections import Counter, defaultdict
from pathlib import Path
import numpy as np
import statsmodels.api as sm

src=Path(sys.argv[1]); out=Path(sys.argv[2]); out.mkdir(parents=True,exist_ok=True)
raw=src.read_bytes()
assert len(raw)==127502
assert hashlib.md5(raw).hexdigest()=="008ec196aa5acb769dbfd36e835e32ba"
assert hashlib.sha256(raw).hexdigest()=="4859b803af66db1f8f18796fe08ce2581363258444de715f1d4851e9f2f988c5"
with zipfile.ZipFile(io.BytesIO(raw)) as z:
    assert z.testzip() is None
    def read(n):
        return list(csv.DictReader(io.StringIO(z.read("data/"+n).decode("utf-8-sig"),newline=""),delimiter="\t"))
    d=read("fitness.txt")
    census=read("frequency.txt")
assert len(d)==601 and len(census)==512
assert len(set(x["id"] for x in d))==601
def val(x,k):
    try:
        y=float(x[k])
        return y if math.isfinite(y) else None
    except (ValueError,TypeError):
        return None
def key(x):
    return tuple(x[k] for k in ("year2","population2","line","quad"))
def eligible(x):
    fields=("no_fruits","no_peduncles_without_fruits","bearing_rate","v","w","ok","ok_p","niwa_t")
    if x["sp"]!="n" or x["morph"] not in ("v","w") or any(val(x,k) is None for k in fields):
        return False
    f,p,v,w,o=[val(x,k) for k in ("no_fruits","no_peduncles_without_fruits","v","w","ok")]
    assert min(f,p,v,w,o)>=0
    assert v+w==val(x,"niwa_t")
    assert o==val(x,"o")+val(x,"k")
    if f+p<=0 or v+w<=0:
        return False
    assert abs(f/(f+p)-val(x,"bearing_rate"))<1e-6
    assert abs(o/(o+v+w)-val(x,"ok_p"))<1e-6
    return True
def bearing(x):
    f,p=[val(x,k) for k in ("no_fruits","no_peduncles_without_fruits")]
    return f/(f+p)
g=defaultdict(list)
for x in d:
    if eligible(x): g[key(x)].append(x)
paired=[]
for k,rows in g.items():
    if set(x["morph"] for x in rows)!={"w","v"}: continue
    # Census fields are stored with every author fitness record, not inferred via a
    # partial join with the separate frequency.txt.
    assert len({tuple(x[c] for c in ("v","w","ok","ok_p","niwa_t")) for x in rows})==1
    a=rows[0]; w=[bearing(x) for x in rows if x["morph"]=="w"]; v=[bearing(x) for x in rows if x["morph"]=="v"]
    paired.append(dict(key="|".join(k),site=k[1],n_white=len(w),n_purple=len(v),
        p_white=val(a,"w")/val(a,"niwa_t"),other_fraction=val(a,"ok_p"),
        total_count=val(a,"niwa_t")+val(a,"ok"),
        white_rate=float(np.mean(w)),purple_rate=float(np.mean(v)),
        delta_pp=100*float(np.mean(w)-np.mean(v))))
assert sum(map(len,g.values()))==570 and len(g)==86
assert len(paired)==32 and sum(x["n_white"]+x["n_purple"] for x in paired)==282
sites=sorted(set(x["site"] for x in paired))
assert len(sites)==5
names=sites+["white_frequency_per_10pp","other_fraction_per_10pp","ln_total_count"]
X=np.array([[int(s==x["site"]) for s in sites]+[10*x["p_white"],10*x["other_fraction"],math.log(x["total_count"])] for x in paired])
Y=np.array([x["delta_pp"] for x in paired])
assert np.linalg.matrix_rank(X)==len(names)
fit=sm.OLS(Y,X).fit(cov_type="HC3",use_t=True)
ci=fit.conf_int()
coefs=[dict(term=term,estimate_pp=float(fit.params[i]),hc3_se=float(fit.bse[i]),
    ci95=[float(z) for z in ci[i]],p_two_sided=float(fit.pvalues[i])) for i,term in enumerate(names)]
assert abs(coefs[-3]["estimate_pp"]-(-0.8871916709145705))<1e-7
assert abs(coefs[-2]["estimate_pp"]-(-1.0676761597212447))<1e-7
with (out/"paired_quadrats.csv").open("w") as f:
    writer=csv.DictWriter(f,fieldnames=list(paired[0]));writer.writeheader();writer.writerows(paired)
result=dict(schema="fcp_sisyrinchium_within_quadrat_fruit_set_frequency_v1",
    original_doi="10.5061/dryad.j8q34",article_doi="10.1371/journal.pone.0164381",
    source_sha256=hashlib.sha256(raw).hexdigest(),source_rows=len(d),frequency_rows=len(census),
    n_eligible=sum(map(len,g.values())),n_eligible_quadrats=len(g),
    n_paired_quadrats=len(paired),n_paired_records=sum(x["n_white"]+x["n_purple"] for x in paired),
    site_quadrats=dict(Counter(x["site"] for x in paired)),df_resid=int(fit.df_resid),
    response="100*(mean white individual-shoot fruit set - mean purple individual-shoot fruit set) within quadrat",
    formula="delta_pp ~ 0 + site + 10*p_white + 10*author_ok_p + log(total_count)",
    coefs=coefs,limitations=["Exploratory observational reanalysis, not frequency manipulation or selection coefficient.",
        "Equal quadrat weights; HC3 addresses unequal quadrat variance, not residual spatial dependence.",
        "Five sites, 32 paired quadrats; site-level generalization is not justified.",
        "Species filtered by original code sp=n and morph=w/v; census counts may not describe the identical plants sampled for fitness.",
        "Shore distance omitted because many paired records have missing values; unmeasured local environment remains a potential confounder.",
        "Female fruit set is not lifetime or male fitness. Original plant genotype and survival not observed."])
(out/"result.json").write_text(json.dumps(result,indent=2)+"\n")
print(json.dumps(result,indent=2))
