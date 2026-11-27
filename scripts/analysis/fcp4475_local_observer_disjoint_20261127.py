#!/usr/bin/env python3
"""Independent-observer local colour detection vs held-out geography and thermal white slope.

Post hoc source-reanalysis of PR #275. The two directions are separately
reported; no pooled treatment of repeated taxa as independent. All original
taxa remain in the ledger. Local PHOTO differences are not genetic morphs.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd

from fcp_full4475_hypothesis_census_20261021 import (
    SOURCE_SHA, BIO_ARCHIVE_SHA, TOTAL_ROWS, TOTAL_TAXA, CLASSIFIABLE_ROWS,
    PHOTO_LABELS, load_original, sample_climate, all_coords_ok, distance_km,
)
from fcp_all4475_seven_abiotic_axes_20261021 import add_srad, SRAD_ZIP_SHA
from fcp4475_crossfitted_observer_independent_thermal_white_20261024 import (
    observer_folds, withheld_slope,
)

SEED = 20261127
N_PERM = 199
N_BOOT = 1999
MIN_GROUP = 20
MIN_GENERA = 15

def local_measure(g, taxon, fold):
    """Exact prior <=5-km same-year/month distinct-observer rule on fold A only.
    Shuffle colours within calendar month (across years) on fixed photo vertices.
    Also record exact-year/month-conditioned expectation as a stricter diagnostic.
    """
    s = g.loc[all_coords_ok(g) & g.date.notna().to_numpy()].reset_index(drop=True)
    counts = g.morph.value_counts()
    repeated = int((counts >= 3).sum()) >= 2
    out = dict(local_classified_n=len(g), local_valid_n=len(s),
               local_observers=int(g.observer_id.nunique()),
               local_repeated3=bool(repeated), local_pairs=0, local_endpoints=0,
               local_strict_pairs=0, local_detected=None, local_mismatch=None,
               local_strict_detected=None, local_expected_detection_month=None,
               local_expected_mismatch_month=None, local_excess_month=None,
               local_expected_detection_yearmonth=None,
               local_expected_mismatch_yearmonth=None)
    if len(s) < 2:
        return out
    u, v = np.triu_indices(len(s), 1)
    dist = distance_km(s.latitude.to_numpy(float), s.longitude.to_numpy(float), u, v)
    year = s.date.dt.year.to_numpy()
    month = s.date.dt.month.to_numpy()
    ym = year * 12 + month
    obs = s.observer_id.to_numpy()
    obid = s.observation_id.to_numpy()
    pair = (dist <= 5) & (ym[u] == ym[v]) & (obs[u] != obs[v])
    uu, vv = u[pair], v[pair]
    out["local_pairs"] = len(uu)
    if not len(uu):
        return out
    labels = pd.Categorical(s.morph, categories=PHOTO_LABELS).codes
    assert (labels >= 0).all()
    diff = labels[uu] != labels[vv]
    out["local_endpoints"] = len(np.unique(np.r_[uu, vv]))
    out["local_detected"] = int(diff.any())
    out["local_mismatch"] = float(diff.mean())
    acc = s.positional_accuracy_m.to_numpy(float)
    strict = (acc[uu] <= 1000) & (acc[vv] <= 1000) & pd.notna(obid[uu]) & pd.notna(obid[vv]) & (obid[uu] != obid[vv])
    out["local_strict_pairs"] = int(strict.sum())
    if strict.any():
        out["local_strict_detected"] = int(diff[strict].any())
    for key, strata in (("month", month), ("yearmonth", ym)):
        rng = np.random.default_rng(SEED + int(taxon)*31 + fold*73 + (key == "yearmonth")*101)
        permutations = np.tile(labels, (N_PERM, 1))
        for m in np.unique(strata):
            ix = np.flatnonzero(strata == m)
            if len(ix) > 1:
                permutations[:, ix] = labels[ix][np.argsort(rng.random((N_PERM, len(ix))), axis=1)]
        z = permutations[:, uu] != permutations[:, vv]
        out["local_expected_detection_" + key] = float(z.any(axis=1).mean())
        out["local_expected_mismatch_" + key] = float(z.mean())
    out["local_excess_month"] = out["local_mismatch"] - out["local_expected_mismatch_month"]
    return out

def geographic_measure(g):
    """Original 50-km range-wide minus local four-label mismatch estimator on fold B."""
    s = g.loc[all_coords_ok(g)].reset_index(drop=True)
    out = dict(geo_valid_n=len(s), geo_pairs50=0, geo_global_mismatch=None,
               geo_delta=None, geo_independent_pairs50=0, geo_independent_delta=None)
    if len(s) < 2:
        return out
    u, v = np.triu_indices(len(s), 1)
    d = distance_km(s.latitude.to_numpy(float), s.longitude.to_numpy(float), u, v)
    labels = pd.Categorical(s.morph, categories=PHOTO_LABELS).codes
    assert (labels >= 0).all()
    mismatch = labels[u] != labels[v]
    near = d <= 50
    out["geo_pairs50"] = int(near.sum())
    out["geo_global_mismatch"] = float(mismatch.mean())
    if near.any():
        out["geo_delta"] = float(mismatch.mean() - mismatch[near].mean())
    obs = s.observer_id.to_numpy()
    distinct = obs[u] != obs[v]
    out["geo_independent_pairs50"] = int((near & distinct).sum())
    if (near & distinct).any():
        out["geo_independent_delta"] = float(mismatch[distinct].mean() - mismatch[near & distinct].mean())
    return out

def source_ledger(raw, classified):
    # Fold allocation uses all original photos BEFORE any colour/climate filtering.
    by_class = {int(t): g for t, g in classified.groupby("inat_taxon_id", sort=False)}
    empty = classified.iloc[:0]
    rows = []
    for k, (tid, allg) in enumerate(raw.groupby("inat_taxon_id", sort=True)):
        cl = by_class.get(int(tid), empty)
        g = cl.loc[cl.observer_id.notna()].copy()
        fold_map = observer_folds(allg)
        g["fold"] = g.observer_id.map(fold_map)
        assert g.fold.isin([0, 1]).all()
        parts = [g.loc[g.fold == f].reset_index(drop=True) for f in (0, 1)]
        assert not (set(parts[0].observer_id) & set(parts[1].observer_id))
        assert not (set(parts[0].photo_id) & set(parts[1].photo_id))
        # Observer separation is within a species. Source observation IDs
        # must not cross roles either; fail closed if malformed source metadata.
        assert not (set(parts[0].observation_id.dropna()) & set(parts[1].observation_id.dropna()))
        for a in (0, 1):
            b = 1 - a
            row = dict(inat_taxon_id=int(tid), species=str(allg.species.iloc[0]),
                       genus=str(allg.species.iloc[0]).split()[0],
                       source_cohort=str(allg.fcp_source_cohort.iloc[0]),
                       direction="forward" if a == 0 else "reverse",
                       local_fold=a, heldout_fold=b, source_classified_n=len(cl),
                       heldout_classified_n=len(parts[b]),
                       heldout_observers=int(parts[b].observer_id.nunique()))
            row.update(local_measure(parts[a], tid, a))
            row.update(geographic_measure(parts[b]))
            # The existing PR #271 climate estimator is reused without modification.
            # Skip climate fitting where local contrast cannot be evaluated;
            # the missing result is explicitly a scope gate, not a zero slope.
            if row["local_pairs"] > 0:
                thermal, slope, _ = withheld_slope(parts[b], tid, a, False)
                row["thermal_white"] = slope
                row["thermal_status"] = thermal["status"]
                row["thermal_complete"] = thermal["heldout_n_complete"]
                row["thermal_cells"] = thermal["heldout_n_cells"]
                row["thermal_unique_var_frac"] = thermal["thermal_unique_var_frac"]
            else:
                row.update(thermal_white=None, thermal_status="local_opportunity_absent_not_fitted",
                           thermal_complete=None, thermal_cells=None, thermal_unique_var_frac=None)
            rows.append(row)
        if (k + 1) % 1000 == 0:
            print("SOURCE_TAXA", k + 1, flush=True)
    d = pd.DataFrame(rows)
    assert len(d) == TOTAL_TAXA * 2
    assert d.groupby("direction").size().eq(TOTAL_TAXA).all()
    assert not d.duplicated(["inat_taxon_id", "direction"]).any()
    return d

def fit_contrast(d, outcome, adjusted, seed):
    """Taxon-equal regression, whole-genus bootstrap; marginal intervals."""
    x = d.loc[d[outcome].notna() & d.local_detected.notna()].copy()
    if adjusted:
        # Only estimated null probabilities with nontrivial detection opportunity.
        # The null expectation is itself estimated from the baseline colours.
        x = x.loc[x.local_expected_detection_month.between(.05, .95)].copy()
    n0 = int((x.local_detected == 0).sum())
    n1 = int((x.local_detected == 1).sum())
    G = int(x.genus.nunique())
    result = dict(n=len(x), n_detected=n1, n_not_detected=n0, n_genera=G,
                  mean_detected=float(x.loc[x.local_detected == 1, outcome].mean()) if n1 else None,
                  mean_not_detected=float(x.loc[x.local_detected == 0, outcome].mean()) if n0 else None)
    if min(n0,n1) < MIN_GROUP or G < MIN_GENERA:
        return dict(result, status="insufficient_group_or_genus_support")
    y = x[outcome].to_numpy(float)
    cols = [np.ones(len(x)), x.local_detected.to_numpy(float)]
    if adjusted:
        cov = np.column_stack([x.local_expected_detection_month,
                               np.log1p(x.local_pairs), np.log1p(x.local_valid_n),
                               np.log1p(x.heldout_classified_n)])
        for z in cov.T:
            if z.std() > 1e-10:
                cols.append((z-z.mean())/z.std())
        # Fixed source cohort controls, observed in this panel only.
        cohort = pd.get_dummies(x.source_cohort, drop_first=True, dtype=float)
        cols.extend(cohort.to_numpy().T)
    X = np.column_stack(cols)
    p = X.shape[1]
    if np.linalg.matrix_rank(X) < p:
        return dict(result, status="rank_deficient")
    beta = np.linalg.lstsq(X,y,rcond=None)[0]
    codes, genera = pd.factorize(x.genus, sort=True)
    gram = np.zeros((G,p,p))
    rhs = np.zeros((G,p))
    np.add.at(gram,codes,np.einsum("ni,nj->nij", X,X))
    np.add.at(rhs,codes,X*y[:,None])
    rng = np.random.default_rng(seed)
    estimates=[]
    for _ in range(N_BOOT):
        w = np.bincount(rng.integers(G,size=G), minlength=G)
        a = np.einsum("g,gij->ij",w,gram)
        b = w @ rhs
        if np.linalg.matrix_rank(a) == p:
            estimates.append(float(np.linalg.solve(a,b)[1]))
    if len(estimates) < .9*N_BOOT:
        return dict(result,status="bootstrap_rank_unstable",valid_boot=len(estimates))
    return dict(result,status="estimable",difference=float(beta[1]),
                ci95=[float(q) for q in np.quantile(estimates,[.025,.975])],
                valid_boot=len(estimates), adjusted=adjusted)

def summary(d):
    panels={}
    for i, direction in enumerate(("forward","reverse")):
        x=d.loc[d.direction == direction]
        local=x.loc[x.local_pairs > 0].copy()
        common=local.loc[local.geo_delta.notna() & local.thermal_white.notna()].copy()
        row=dict(source_taxa=len(x),local_opportunity=len(local),
                 local_detected=int(local.local_detected.sum()),
                 common_n=len(common),common_detected=int(common.local_detected.sum()),
                 common_both_positive_detected=int(((common.local_detected==1)&(common.geo_delta>0)&(common.thermal_white>0)).sum()),
                 exact_yearmonth_null_deterministic=int(local.local_expected_detection_yearmonth.isin([0,1]).sum()),
                 month_null_nontrivial=int(local.local_expected_detection_month.between(.05,.95).sum()))
        row["tests"]={}
        for j,(key,z) in enumerate((("all_local_opportunity",local),
                                  ("baseline_repeated3",local.loc[local.local_repeated3]),
                                  ("common_geo_thermal",common))):
            tests={}
            for k,outcome in enumerate(("geo_delta","thermal_white")):
                for adj in (False,True):
                    tests[outcome+("_adjusted" if adj else "_unadjusted")]=fit_contrast(z,outcome,adj,SEED+i*1000+j*100+k*10+int(adj))
            row["tests"][key]=tests
        panels[direction]=row
    return dict(source_photo_sha256=SOURCE_SHA, source_bio_sha256=BIO_ARCHIVE_SHA,
                source_srad_sha256=SRAD_ZIP_SHA, source_photos=TOTAL_ROWS,
                source_classifiable=CLASSIFIABLE_ROWS, source_taxa=TOTAL_TAXA,
                permutations=N_PERM, genus_bootstrap=N_BOOT,
                roles="Observer-disjoint within taxon, A local detection and B geography plus temperature; reverse roles separately. B geography and B temperature are not mutually independent.",
                scope="Exploratory photographic outcomes only. No genetic morph, balancing selection, pigment fitness or causal adaptation inference. No pooling of directions.",
                panels=panels)

def self_test():
    # A local pair is different-coloured, B-only geography is exactly zero
    # for monochromatic photos. The month null is degenerate if all dates
    # share a single year/month and exactly two eligible photo vertices.
    g=pd.DataFrame(dict(latitude=[0.,0.],longitude=[0.,.01],
                        positional_accuracy_m=[100.,100.],
                        date=pd.to_datetime(["2020-06-01","2020-06-05"]),
                        observer_id=[1,2],observation_id=[11,12],
                        morph=["white","red_pink"]))
    a=local_measure(g,1,0)
    assert a["local_pairs"]==1 and a["local_detected"]==1
    assert a["local_expected_detection_month"]==1
    assert a["local_expected_detection_yearmonth"]==1
    g["morph"]="white"
    b=geographic_measure(g)
    assert b["geo_pairs50"]==1 and b["geo_delta"]==0
    print("FCP_LOCAL_OBSERVER_DISJOINT_SELF_TEST_PASS")

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ("photos","worldclim-bio","worldclim-srad","out"):
        p.add_argument("--"+key,type=Path)
    p.add_argument("--self-test",action="store_true")
    a=p.parse_args()
    self_test()
    if a.self_test:
        return
    assert a.photos and a.worldclim_bio and a.worldclim_srad and a.out
    raw,d=load_original(a.photos)
    d=add_srad(sample_climate(d,a.worldclim_bio),a.worldclim_srad)
    assert len(raw)==TOTAL_ROWS and len(d)==CLASSIFIABLE_ROWS
    z=source_ledger(raw,d)
    report=summary(z)
    a.out.mkdir(parents=True,exist_ok=True)
    z.to_csv(a.out/"fcp4475_local_observer_disjoint_ledger.csv.gz",index=False,
             compression={"method":"gzip","mtime":0},float_format="%.12g")
    (a.out/"fcp4475_local_observer_disjoint_results.json").write_text(
        json.dumps(report,indent=2,allow_nan=False)+"\n",encoding="utf-8")
    print("FCP_LOCAL_OBSERVER_DISJOINT_REAL_SOURCE_COMPLETE",json.dumps(report,allow_nan=False),flush=True)

if __name__=="__main__":
    main()
