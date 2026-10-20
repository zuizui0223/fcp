#!/usr/bin/env python3
"""Descriptive WORLD photographic within-species flower-colour atlas, source frozen.

Not an FCP-population prevalence estimate, molecular result, spatial adaptation
claim, historical-range reconstruction, or a replacement for the frozen H2/W.
No exact coordinates or original photo IDs are exported.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd

SHA = "8e7f10f0427927208485b2da6fa9c46bb681d6f3611e53b8dfdc606cfc8dabe6"
TOTAL_ROWS = 447500
TOTAL_SPECIES = 4475
CLASSIFIED_ROWS = 227107
STATES = ("white", "yellow_orange", "red_pink", "blue_purple")
PALETTE = {
    "white": ("flower_fraction_white",),
    "yellow_orange": ("flower_fraction_yellow", "flower_fraction_orange"),
    "red_pink": ("flower_fraction_red", "flower_fraction_pink", "flower_fraction_magenta"),
    "blue_purple": ("flower_fraction_blue", "flower_fraction_purple"),
}
COLUMNS = [
    "inat_taxon_id", "photo_id", "observation_id", "observer_id",
    "latitude", "longitude", "positional_accuracy_m", "morph", "global_classifiable",
    *(f for group in PALETTE.values() for f in group),
]
MIN_PER_UNIT = 10
MIN_OBSERVERS = 3
MIN_GLOBAL = 40
BOOT = 399


def check(ok, msg):
    if not ok:
        raise RuntimeError("FCP_ATLAS_SOURCE_OR_METHOD_INVALID: " + msg)


def digest(p):
    h = hashlib.sha256()
    with Path(p).open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_original(path):
    check(digest(path) == SHA, "SHA256 source mismatch")
    raw = pd.read_csv(path, usecols=COLUMNS, low_memory=False)
    check(len(raw) == TOTAL_ROWS and raw.inat_taxon_id.nunique() == TOTAL_SPECIES,
          "complete real photograph/source species count mismatch")
    check(raw.photo_id.nunique() == TOTAL_ROWS,
          "photo IDs were reused")
    ok = raw.global_classifiable.fillna(False).astype(str).str.lower().isin(("true", "1", "yes"))
    d = raw.loc[ok & raw.morph.isin(STATES)].copy()
    check(len(d) == CLASSIFIED_ROWS, "classified real photo count mismatch")
    return d


def group_stats(d, keys, min_n=MIN_PER_UNIT, min_observers=MIN_OBSERVERS):
    """Unbiased photographed pair difference conditional on unit and inclusion."""
    if len(d) == 0:
        return pd.DataFrame(columns=[*keys,"n_photo","n_observers",*STATES,
                                     "D_pair", "p_white", "white_discordance_share", "n_states"])
    counter = d.groupby([*keys, "morph"], observed=True).size().unstack(fill_value=0)
    counter = counter.reindex(columns=STATES, fill_value=0)
    tally = d.groupby(keys, observed=True).agg(
        n_photo=("photo_id", "size"), n_observers=("observer_id", "nunique"))
    tab = tally.join(counter, how="left").reset_index()
    tab = tab.loc[(tab.n_photo >= min_n) & (tab.n_observers >= min_observers)].copy()
    if tab.empty:
        return tab
    n = tab.n_photo.to_numpy(float)
    counts = tab[list(STATES)].to_numpy(float)
    nonmatching = n * (n-1) - np.sum(counts * (counts-1), axis=1)
    tab["D_pair"] = nonmatching / (n * (n-1))
    tab["p_white"] = tab.white.to_numpy(float) / n
    tab["n_states"] = np.sum(counts > 0, axis=1)
    tab["white_discordance_share"] = np.divide(
        2 * tab.white.to_numpy(float) * (n - tab.white.to_numpy(float)),
        nonmatching, out=np.full(len(tab), np.nan), where=nonmatching > 0)
    return tab


def region_aggregate(tab, key):
    out=[]
    rng=np.random.default_rng(20261020)
    for region, g in tab.groupby(key, observed=True, sort=True):
        n=len(g)
        # Equal weighting of taxa in each geographic observational opportunity unit.
        vals=g.D_pair.to_numpy(float)
        samples=rng.integers(0,n,size=(BOOT,n)) if n>=20 else None
        ci=np.quantile(vals[samples].mean(axis=1), [.025,.975]).tolist() if n>=20 else [None,None]
        out.append({
            key: str(region), "n_species_eligible":int(n),
            "n_photographs_in_eligible_taxon_units":int(g.n_photo.sum()),
            "mean_D_pair_species_equal":float(vals.mean()),
            "mean_p_white_species_equal":float(g.p_white.mean()),
            "species_with_at_least_two_photo_labels":int((g.n_states>=2).sum()),
            "mean_white_share_of_different_colour_pairs":(
                float(g.white_discordance_share.mean()) if g.white_discordance_share.notna().any() else None),
            "bootstrap_D_pair_CI95_species_resampling":ci,
            "bootstrap_available":bool(n>=20),
            "descriptive_map_min_5_taxa":bool(n>=5),
        })
    return out


def pair_boot(differences):
    x=np.array(differences,float)
    if len(x)<15:return [None,None]
    r=np.random.default_rng(20261021)
    idx=r.integers(0,len(x),size=(BOOT,len(x)))
    return [float(t) for t in np.quantile(x[idx].mean(axis=1),[.025,.975])]


def paired_zones(zone):
    if zone.empty:return {"n_spanning_both_zones":0,"status":"not_estimable"}
    t=zone.pivot(index="inat_taxon_id",columns="zone",values="D_pair")
    if not {"tropical","extratropical"}.issubset(t.columns):
        return {"n_spanning_both_zones":0,"status":"not_estimable"}
    q=t[["tropical","extratropical"]].dropna()
    v=q.tropical.to_numpy(float)-q.extratropical.to_numpy(float)
    return {"n_spanning_both_zones":len(v),
            "mean_tropical_minus_extratropical_D_within_same_species":float(v.mean()) if len(v) else None,
            "CI95_species_bootstrap":pair_boot(v),
            "status":"descriptive_only" if len(v)>=15 else "underpowered_not_inferential"}


def geographic_units(d):
    lat=pd.to_numeric(d.latitude,errors="coerce")
    lon=pd.to_numeric(d.longitude,errors="coerce")
    accuracy=pd.to_numeric(d.positional_accuracy_m,errors="coerce")
    valid=(lat.between(-90,90) & lon.between(-180,180)
           & accuracy.ge(0) & accuracy.le(5000))
    z=d.loc[valid].copy()
    a=lat.loc[valid].to_numpy(float)
    b=lon.loc[valid].to_numpy(float)
    y=np.clip(np.floor(9*(np.sin(np.deg2rad(a))+1)/2).astype(int),0,8)
    x=np.clip(np.floor(18*(b+180)/360).astype(int),0,17)
    z["cell"]=[f"EA{u:02d}_{v:02d}" for u,v in zip(y,x)]
    z["abs_lat_band"]=pd.cut(np.abs(a),bins=[0,15,30,45,60,90.00001],
                              include_lowest=True,
                              labels=["0-15","15-30","30-45","45-60","60-90"]).astype(str)
    z["zone"]=np.where(np.abs(a)<23.5,"tropical","extratropical")
    return z


def cell_turnover(cell):
    species=[]
    for tid,g in cell.groupby("inat_taxon_id",observed=True,sort=True):
        if len(g)<2:continue
        fractions=g[list(STATES)].to_numpy(float) / g.n_photo.to_numpy(float)[:,None]
        values=[(float(np.abs(fractions[a]-fractions[b]).sum()/2),
                 float(abs(fractions[a,0]-fractions[b,0])))
                for a,b in itertools.combinations(range(len(g)),2)]
        species.append({"inat_taxon_id":int(tid),
                        "n_equal_area_cells_with_10plus_photos":len(g),
                        "n_cell_pairs_not_independent":len(values),
                        "mean_between_cell_total_variation":float(np.mean([v[0] for v in values])),
                        "mean_between_cell_absolute_white_frequency_difference":float(np.mean([v[1] for v in values]))})
    return pd.DataFrame(species)


def high_purity(d):
    a=np.column_stack([d[list(PALETTE[k])].sum(axis=1).to_numpy(float) for k in STATES])
    check(np.isfinite(a).all(),"ROI fractions nonfinite")
    label=d.morph.to_numpy()
    chosen=np.array([STATES.index(t) for t in label])
    primary=a[np.arange(len(a)),chosen]
    a[np.arange(len(a)),chosen]=-np.inf
    secondary=np.max(a,axis=1)
    return d.loc[(primary>=.65)&(secondary<=.20)].copy()


def selftest():
    rows=[]
    for tid in (1,2):
        for cell in ("EA00_00","EA08_17"):
            for k in range(10):
                rows.append({"inat_taxon_id":tid,"photo_id":len(rows)+1,
                             "observer_id":k//2,"morph":"white" if k<5 else "red_pink",
                             "cell":cell,"zone":"tropical" if cell=="EA00_00" else "extratropical"})
    d=pd.DataFrame(rows)
    v=group_stats(d,["inat_taxon_id","cell"])
    check(len(v)==4,"synthetic region units")
    check(np.allclose(v.D_pair,50/90),"unbiased D = 50/90")
    check(np.allclose(v.white_discordance_share,1),"white discordance bookkeeping")
    check(paired_zones(group_stats(d,["inat_taxon_id","zone"]))["n_spanning_both_zones"]==2,
          "paired-zone distinct species identities")
    print("FCP_WORLD_ATLAS_SYNTHETIC_TESTS_PASS")


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--photos",type=Path)
    p.add_argument("--out",type=Path)
    p.add_argument("--self-test",action="store_true")
    args=p.parse_args()
    if args.self_test:
        selftest();return
    check(args.photos is not None and args.out is not None,"required paths")
    selftest()
    d=read_original(args.photos)
    g=group_stats(d,["inat_taxon_id"],MIN_GLOBAL,MIN_OBSERVERS)
    geo=geographic_units(d)
    cell=group_stats(geo,["inat_taxon_id","cell"])
    band=group_stats(geo,["inat_taxon_id","abs_lat_band"])
    zone=group_stats(geo,["inat_taxon_id","zone"])
    tv=cell_turnover(cell)
    pure=high_purity(d)
    pureg=group_stats(pure,["inat_taxon_id"],MIN_GLOBAL,MIN_OBSERVERS)
    purecell=group_stats(geographic_units(pure),["inat_taxon_id","cell"])
    result={
        "schema":"fcp_original_measured_global_photo_pattern_exploratory_v1",
        "source_sha256":SHA,
        "original_measured_rows":TOTAL_ROWS,
        "original_measured_species":TOTAL_SPECIES,
        "original_classifiable_four_state_rows":CLASSIFIED_ROWS,
        "geo_accurate_le5000m_classifiable_rows":int(len(geo)),
        "global_D_evaluable_species_min40photos_3observers":int(len(g)),
        "global_D_mean_species_equal":float(g.D_pair.mean()),
        "global_species_with_2plus_detected_photo_colour_labels":int((g.n_states>=2).sum()),
        "global_species_with_both_white_and_nonwhite_labels":int(((g.white>0)&(g.white<g.n_photo)).sum()),
        "global_mean_p_white_species_equal":float(g.p_white.mean()),
        "global_mean_white_contribution_to_discordant_pairs":float(g.white_discordance_share.mean()),
        "colour_counts_classifiable_photos":{s:int((d.morph==s).sum()) for s in STATES},
        "n_taxon_cell_units_eligible_10photos_3observers":int(len(cell)),
        "n_species_eligible_any_cell":int(cell.inat_taxon_id.nunique()),
        "n_occupied_cells_by_classifiable_photos":int(geo.cell.nunique()),
        "n_cells_with_eligible_taxon_units":int(cell.cell.nunique()),
        "n_species_spanning_2plus_eligible_cells":int(len(tv)),
        "mean_species_equal_between_cell_total_variation":float(tv.mean_between_cell_total_variation.mean()) if len(tv) else None,
        "mean_species_equal_between_cell_white_delta":float(tv.mean_between_cell_absolute_white_frequency_difference.mean()) if len(tv) else None,
        "absolute_latitude_bands":region_aggregate(band,"abs_lat_band"),
        "tropical_extratropical_opportunities":region_aggregate(zone,"zone"),
        "paired_trop_extratrop_same_species":paired_zones(zone),
        "equal_area_grid_cells":region_aggregate(cell,"cell"),
        "quality_sensitivity":{
            "definition":"primary coarse colour ROI fraction >=0.65; every competing coarse group <=0.20; not a biological morph-validation score",
            "photos":int(len(pure)),"global_species_min40":int(len(pureg)),
            "global_D_mean_species_equal":float(pureg.D_pair.mean()) if len(pureg) else None,
            "equal_area_taxon_cell_units_min10":int(len(purecell)),
            "n_species_any_eligible_pure_cell":int(purecell.inat_taxon_id.nunique()) if len(purecell) else 0,
        },
        "strongest_allowed_claim":"observed photographic flower-colour contrast across sampled taxa/regions with equal taxon weight; NOT real local genetic FCP or unbiased global prevalence",
        "historical_vs_current_environment_attribution":"NOT IDENTIFIABLE without independent history/genetics and site-time environment data",
        "still_unverified":["native/wild status", "flower age/part identity", "true same-individual morphs", "independent camera/source replication", "genetic IBD", "genetic IBE", "phenotype fitness"],
        "scope":"post hoc world-pattern exploration from already-measured original-method photographs; distinct from pre-frozen W test",
    }
    args.out.mkdir(parents=True,exist_ok=True)
    with (args.out/"global_photo_pattern_result.json").open("w") as f:
        json.dump(result,f,indent=2,ensure_ascii=False,allow_nan=False);f.write("\n")
    tv.to_csv(args.out/"taxon_between_grid_cell_turnover_no_coordinates.csv",index=False)
    cell.to_csv(args.out/"taxon_equal_area_cell_photo_colour_no_coordinates.csv",index=False)
    print("FCP_GLOBAL_PATTERN_ACTUAL_SOURCE_COMPLETE",json.dumps({
        k:v for k,v in result.items() if k not in ("equal_area_grid_cells",)},
        ensure_ascii=False,sort_keys=True,allow_nan=False))
    print("FCP_GLOBAL_PATTERN_GRID_TOP_COVERAGE", json.dumps(sorted(result["equal_area_grid_cells"],
         key=lambda z:z["n_species_eligible"],reverse=True)[:8],allow_nan=False))


if __name__=="__main__":main()
