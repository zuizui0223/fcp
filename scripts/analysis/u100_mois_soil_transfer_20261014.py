"""Post hoc within-species FCP soil/abiotic opportunity transfer from mois.

Source authority: original 100,543-photo, 36 environmental-predictor archive,
with original photographed site coordinates (SoilGrids/WorldClim). Never merge
measurements by species name, average genus, or assume nearest-source soil is
an exact location measurement. Two prespecified proximity bands are reported;
not calibrated prospectively and not species-selection independent.
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import wilcoxon

ENV_SHA = "6e287294b3a05eaa383be3de644551825915bfecca9c4bde6950b4890f5e16ca"
FCP_MEASURED_SHA = "57630fc9f281bce94a0c40a70aaf7bce879dde93d6154175adcd021e8f5c1186"
SOIL = ("soil_pH", "soil_SOC", "soil_N", "soil_clay", "soil_available_water_proxy",
        "soil_cec_0_30cm_source_raw", "soil_sand_0_30cm_source_raw",
        "soil_silt_0_30cm_source_raw", "soil_bdod_0_30cm_source_raw",
        "soil_cfvo_0_30cm_source_raw")
BLOCKS = {
    "elevation": ("wc_elevation_m",),
    "temperature": tuple(f"wc_bio{i}" for i in range(1, 12)),
    "precipitation": tuple(f"wc_bio{i}" for i in range(12, 20)),
    "solar_radiation": ("wc_srad_annual_mean", "wc_srad_monthly_cv"),
    "wind": ("wc_wind_annual_mean", "wc_wind_monthly_cv"),
    "vapor_pressure": ("wc_vapr_annual_mean", "wc_vapr_monthly_cv"),
    "soil": SOIL,
}
FEATURES = tuple(x for block in BLOCKS.values() for x in block)
assert len(FEATURES) == 36
PHOTO_STATES = ("white", "yellow_orange", "red_pink", "blue_purple")
SOURCE_MEMBER = ("fcp42111_original_photo_all_climate_soil/"
                 "all_original_100543_photo_sites_expanded_environment.csv.gz")
EARTH_RADIUS_KM = 6371.0088
MIN_COLOUR = 5
MIN_SPECIES = 30
RADII_KM = (1.0, 5.0)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def bool_column(s: pd.Series) -> pd.Series:
    return s.astype(str).str.strip().str.lower().isin({"true", "1", "yes"})


def source_table(path: Path, verify: bool = True) -> pd.DataFrame:
    if verify and sha256(path) != ENV_SHA:
        raise ValueError("mois original 36-abiotic source archive SHA256 mismatch")
    with zipfile.ZipFile(path) as z:
        if SOURCE_MEMBER not in z.namelist():
            raise ValueError("original photo-coordinate soil source member is missing")
        cols = ["inat_taxon_id", "photo_id", "latitude", "longitude",
                "positional_accuracy_m", "site_geo_status", *FEATURES]
        d = pd.read_csv(io.BytesIO(z.read(SOURCE_MEMBER)), compression="gzip",
                        usecols=cols, low_memory=False, dtype={"photo_id": str})
    if d.photo_id.duplicated().any() or (verify and len(d) != 100543):
        raise ValueError("mois source ID / 100543-photo denominator changed")
    if verify and d.inat_taxon_id.nunique() != 42111:
        raise ValueError("mois opportunity taxa changed")
    return d


def photo_table(measured: Path, technical: Path, join_key: Path, high_clip: Path,
                verify: bool = True) -> tuple[pd.DataFrame, dict]:
    if verify and sha256(measured) != FCP_MEASURED_SHA:
        raise ValueError("FCP frozen prospective photo SHA256 mismatch")
    wanted = ["photo_id", "inat_taxon_id", "species", "observer_id", "morph",
              "global_classifiable", "latitude", "longitude", "positional_accuracy_m"]
    d = pd.read_csv(measured, usecols=wanted, dtype={"photo_id": str}, low_memory=False)
    if d.photo_id.duplicated().any():
        raise ValueError("FCP prospective photograph duplicate")
    d = d.loc[bool_column(d.global_classifiable) & d.morph.isin(PHOTO_STATES)].copy()
    tech = pd.read_csv(technical, compression="gzip", dtype={"measurement_id": str})
    keys = pd.read_csv(join_key, dtype={"measurement_id": str, "photo_id": str})
    high = pd.read_csv(high_clip, dtype={"measurement_id": str})
    if tech.measurement_id.duplicated().any() or keys.measurement_id.duplicated().any() or keys.photo_id.duplicated().any():
        raise ValueError("technical measurement join is not one-to-one")
    tech = tech.merge(keys[["measurement_id", "photo_id"]], on="measurement_id", how="left", validate="one_to_one")
    tech["high_clip"] = tech.measurement_id.isin(set(high.measurement_id.astype(str)))
    if d.merge(tech[["photo_id", "near_clip_fraction"]], on="photo_id", how="left")["near_clip_fraction"].notna().mean() < .85:
        raise ValueError("near-clipping technical coverage below 85%")
    d = d.merge(tech[["photo_id", "near_clip_fraction", "high_clip"]], on="photo_id", how="left", validate="one_to_one")
    d = d.loc[d.near_clip_fraction.notna() & ~d.high_clip.fillna(False)].copy()
    d["white"] = d.morph.eq("white")
    count = d.groupby("inat_taxon_id").white.agg(["sum", "count"])
    good = count.index[(count["sum"] >= MIN_COLOUR) & ((count["count"] - count["sum"]) >= MIN_COLOUR)]
    d = d.loc[d.inat_taxon_id.isin(good)].copy()
    receipt = {"photo_rows_after_original_gates": int(len(d)),
               "eligible_taxa_after_original_gates": int(d.inat_taxon_id.nunique()),
               "source_photo_sha256": sha256(measured)}
    if verify and (receipt["photo_rows_after_original_gates"] != 12583 or
                   receipt["eligible_taxa_after_original_gates"] != 281):
        raise ValueError("FCP 2026-09 frozen 281-species eligibility mismatch")
    return d, receipt


def unit_xyz(lat: np.ndarray, lon: np.ndarray) -> np.ndarray:
    latr, lonr = np.deg2rad(lat), np.deg2rad(lon)
    return np.column_stack((np.cos(latr) * np.cos(lonr),
                            np.cos(latr) * np.sin(lonr), np.sin(latr)))


def nearest_site_map(photos: pd.DataFrame, env: pd.DataFrame, max_accuracy_m: float = 1000.0) -> pd.DataFrame:
    """Nearest geographically screened original photo site, independent of its soil values.

    A site is chosen BEFORE checking environment availability: avoids a
    soil-completeness-driven nearest-neighbour search. Exact photo matches
    additionally require the same taxon ID and compatible coordinates.
    """
    source = env.reset_index(drop=True)
    source_ids = source.photo_id.astype(str)
    valid = (source.latitude.between(-90, 90) & source.longitude.between(-180, 180) &
             pd.to_numeric(source.positional_accuracy_m, errors="coerce").between(0, max_accuracy_m))
    eligible_indices = np.flatnonzero(valid.to_numpy(bool))
    if len(eligible_indices) < 2:
        raise ValueError("insufficient original georeferenced source sites")
    coords = unit_xyz(source.loc[valid, "latitude"].to_numpy(float),
                      source.loc[valid, "longitude"].to_numpy(float))
    geo_good = (photos.latitude.between(-90, 90) & photos.longitude.between(-180, 180) &
                pd.to_numeric(photos.positional_accuracy_m,errors="coerce").between(0,max_accuracy_m))
    best_index = np.full(len(photos), -1, dtype=int)
    best_dist = np.full(len(photos), np.nan, dtype=float)
    if geo_good.any():
        targets = unit_xyz(photos.loc[geo_good, "latitude"].to_numpy(float),
                           photos.loc[geo_good, "longitude"].to_numpy(float))
        chord, near = cKDTree(coords).query(targets, workers=2)
        best_index[np.flatnonzero(geo_good.to_numpy(bool))] = eligible_indices[near]
        best_dist[np.flatnonzero(geo_good.to_numpy(bool))] = (2 * EARTH_RADIUS_KM *
            np.arcsin(np.clip(chord / 2.0, 0, 1)))
    exact_index = {str(pid): i for i, pid in enumerate(source_ids)}
    exact_photo = np.zeros(len(photos), dtype=bool)
    for i, (pid, taxon) in enumerate(zip(photos.photo_id.astype(str), photos.inat_taxon_id)):
        j = exact_index.get(pid)
        if j is None:
            continue
        if int(source.iloc[j].inat_taxon_id) != int(taxon):
            raise ValueError("exact photo_id has conflicting taxon across sources")
        if valid.iloc[j] and geo_good.iloc[i]:
            x1 = unit_xyz(np.array([float(photos.iloc[i].latitude)]),
                          np.array([float(photos.iloc[i].longitude)]))
            x2 = coords[np.where(eligible_indices == j)[0][0]]
            sep = float(2 * EARTH_RADIUS_KM * np.arcsin(np.clip(np.linalg.norm(x1[0] - x2) / 2, 0, 1)))
            if sep > 0.1:
                raise ValueError("same photo ID has materially different source coordinates")
            best_index[i], best_dist[i], exact_photo[i] = j, 0.0, True
    mapped = photos.reset_index(drop=True).copy()
    mapped["source_index"] = best_index
    mapped["distance_to_source_km"] = best_dist
    mapped["source_match_exact_photo"] = exact_photo
    mapped["transfer_1km"] = exact_photo | (best_dist <= 1)
    mapped["transfer_5km"] = exact_photo | (best_dist <= 5)
    for col in FEATURES:
        buf = np.full(len(mapped), np.nan, dtype=float)
        good = best_index >= 0
        buf[good] = pd.to_numeric(source.iloc[best_index[good]][col], errors="coerce").to_numpy(float)
        mapped[col] = buf
    return mapped


def holm_adjust(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, dtype=float)
    p = np.where(np.isfinite(p), p, 1.0)
    ix = np.argsort(p, kind="stable")
    adj = np.minimum(1, np.maximum.accumulate((len(p) - np.arange(len(p))) * p[ix]))
    out = np.empty_like(adj)
    out[ix] = adj
    return out


def per_species_effect(data: pd.DataFrame, var: str) -> pd.DataFrame:
    lines = []
    for (taxon, sp), g in data.loc[np.isfinite(data[var]), ["inat_taxon_id", "species", "observer_id", "latitude", "longitude", "white", var]].groupby(["inat_taxon_id", "species"], sort=True):
        wn = int(g.white.sum()); cn = int((~g.white).sum())
        if min(wn, cn) < MIN_COLOUR:
            continue
        v = g[var].to_numpy(float)
        scale = float(np.std(v, ddof=0))
        if scale <= 1e-12:
            continue
        g = g.copy()
        g["z"] = (g[var] - np.mean(v)) / scale
        delta = float(g.loc[g.white, "z"].mean() - g.loc[~g.white, "z"].mean())
        paired_obs=[]
        for _, h in g.dropna(subset=["observer_id"]).groupby("observer_id"):
            if h.white.any() and (~h.white).any():
                paired_obs.append(float(h.loc[h.white, "z"].mean()-h.loc[~h.white,"z"].mean()))
        cell_lat = np.floor(g.latitude.to_numpy(float))
        cell_lon = np.floor(g.longitude.to_numpy(float))
        g["cell"] = [f"{int(x)}:{int(y)}" for x,y in zip(cell_lat,cell_lon)]
        paired_cells=[]
        for _, h in g.groupby("cell"):
            if h.white.any() and (~h.white).any():
                paired_cells.append(float(h.loc[h.white,"z"].mean()-h.loc[~h.white,"z"].mean()))
        lines.append({"inat_taxon_id": int(taxon), "species": str(sp), "variable": var,
                      "n_white": wn, "n_nonwhite": cn, "delta_white_minus_nonwhite_sd": delta,
                      "paired_observers": len(paired_obs),
                      "observer_paired_delta": float(np.mean(paired_obs)) if len(paired_obs)>=3 else np.nan,
                      "paired_1degree_cells": len(paired_cells),
                      "same_cell_delta": float(np.mean(paired_cells)) if paired_cells else np.nan})
    return pd.DataFrame(lines)


def all_species_effects(data: pd.DataFrame, features: tuple[str, ...] = FEATURES) -> pd.DataFrame:
    """Vectorize 36 environmental effects within each species, with observer/cell checks."""
    out = []
    for (taxon, sp), g in data.groupby(["inat_taxon_id", "species"], sort=True):
        x = g[list(features)].to_numpy(dtype=float)
        white = g["white"].to_numpy(bool)
        finite = np.isfinite(x)
        nw = finite[white].sum(axis=0)
        nc = finite[~white].sum(axis=0)
        good = (nw >= MIN_COLOUR) & (nc >= MIN_COLOUR)
        n_all = finite.sum(axis=0)
        sum_all = np.where(finite, x, 0.0).sum(axis=0)
        means = np.divide(sum_all, n_all, out=np.zeros(len(features)), where=n_all>0)
        variance = np.divide(np.where(finite, (x - means)**2, 0.0).sum(axis=0),
                             n_all, out=np.zeros(len(features)), where=n_all>0)
        sd = np.sqrt(variance)
        good &= np.isfinite(sd) & (sd > 1e-12)
        if not good.any():
            continue
        wmean = np.divide(np.where(finite & white[:,None],x,0.0).sum(axis=0),nw,
                          out=np.zeros(len(features)),where=nw>0)
        cmean = np.divide(np.where(finite & (~white)[:,None],x,0.0).sum(axis=0),nc,
                          out=np.zeros(len(features)),where=nc>0)
        delta = np.divide(wmean-cmean,sd,out=np.zeros(len(features)),where=sd>0)

        def matched_groups(codes: np.ndarray, min_groups: int) -> tuple[np.ndarray,np.ndarray]:
            valid_codes = codes>=0
            n_groups = int(np.max(codes[valid_codes])+1) if np.any(valid_codes) else 0
            if n_groups == 0:
                return np.zeros(len(features), dtype=int), np.full(len(features), np.nan)
            accum = np.zeros((2,n_groups,len(features)),float)
            num = np.zeros((2,n_groups,len(features)),int)
            for j,which in enumerate((white & valid_codes,(~white) & valid_codes)):
                np.add.at(accum[j],codes[which],np.where(finite[which],x[which],0.0))
                np.add.at(num[j],codes[which],finite[which].astype(int))
            matched = (num[0] > 0) & (num[1] > 0)
            n_match = matched.sum(axis=0)
            obs_w = np.divide(accum[0],num[0],out=np.zeros_like(accum[0]),where=num[0]>0)
            obs_c = np.divide(accum[1],num[1],out=np.zeros_like(accum[1]),where=num[1]>0)
            total = np.where(matched, (obs_w-obs_c),0.0).sum(axis=0)
            mean = np.divide(total,n_match,out=np.zeros(len(features)),where=n_match>0)
            z = np.divide(mean,sd,out=np.full(len(features),np.nan),where=(sd>0))
            z[n_match<min_groups] = np.nan
            return n_match,z

        observers = g["observer_id"].astype("string")
        obs_codes,uniques = pd.factorize(observers, sort=True)
        n_obs, obs_delta = matched_groups(obs_codes,3)
        cells = np.floor(g[["latitude","longitude"]].to_numpy(float))
        cell_labels = [f"{int(a)}:{int(b)}" for a,b in cells]
        cell_codes,_ = pd.factorize(np.asarray(cell_labels,dtype=str),sort=True)
        n_cells,cell_delta = matched_groups(np.asarray(cell_codes),1)
        for j,var in enumerate(features):
            if good[j]:
                out.append({"inat_taxon_id":int(taxon),"species":str(sp),"variable":var,
                    "n_white":int(nw[j]),"n_nonwhite":int(nc[j]),
                    "delta_white_minus_nonwhite_sd":float(delta[j]),
                    "paired_observers":int(n_obs[j]),
                    "observer_paired_delta":float(obs_delta[j]) if np.isfinite(obs_delta[j]) else np.nan,
                    "paired_1degree_cells":int(n_cells[j]),
                    "same_cell_delta":float(cell_delta[j]) if np.isfinite(cell_delta[j]) else np.nan})
    return pd.DataFrame(out,columns=["inat_taxon_id","species","variable","n_white","n_nonwhite",
                                     "delta_white_minus_nonwhite_sd","paired_observers",
                                     "observer_paired_delta","paired_1degree_cells","same_cell_delta"])


def stats_effects(rows: pd.DataFrame, var: str, group: str, band: str) -> dict:
    if rows.empty:
        return {"variable":var,"block":group,"transfer_band":band,"species":0,
                "status":"not_estimable","nominal_p":float("nan"),
                "observer_paired_species":0,"observer_paired_nominal_p":float("nan"),
                "same_cell_species":0,"same_cell_nominal_p":float("nan")}
    a=rows["delta_white_minus_nonwhite_sd"].to_numpy(float)
    a=a[np.isfinite(a)]
    def paired_summary(name):
        vals=rows[name].dropna().to_numpy(float)
        nominal=float(wilcoxon(vals,alternative="two-sided").pvalue) if len(vals)>=MIN_SPECIES and np.any(vals!=0) else float("nan")
        return {"n":int(len(vals)),
                "median":float(np.median(vals)) if len(vals) else None,"p":nominal}
    obs=paired_summary("observer_paired_delta")
    geo=paired_summary("same_cell_delta")
    p=float(wilcoxon(a,alternative="two-sided").pvalue) if len(a)>=MIN_SPECIES and np.any(a!=0) else float("nan")
    return {"variable":var,"block":group,"transfer_band":band,
            "species":int(len(a)), "median_sd":float(np.median(a)) if len(a) else None,
            "mean_sd":float(np.mean(a)) if len(a) else None,
            "positive_fraction":float(np.mean(a>0)) if len(a) else None,
            "nominal_p":p,"status":"exploratory" if np.isfinite(p) else "not_estimable",
            "observer_paired_species":obs["n"],
            "observer_paired_median_sd":obs["median"],
            "observer_paired_nominal_p":obs["p"],
            "same_cell_species":geo["n"],
            "same_cell_median_sd":geo["median"],
            "same_cell_nominal_p":geo["p"]}


def run_analysis(photo: pd.DataFrame, source: pd.DataFrame, outdir: Path,
                 max_accuracy_m: float=1000.0, radii=RADII_KM) -> dict:
    mapped=nearest_site_map(photo,source,max_accuracy_m=max_accuracy_m)
    bands=[]; effects=[]
    for radius in (0.0,*radii):
        label = "exact" if radius == 0 else f"{radius:g}km"
        sub = mapped.loc[mapped["source_match_exact_photo"] if radius == 0 else
                         mapped[f"transfer_{int(radius)}km"]].copy()
        byvar = {v:int(np.isfinite(sub[v]).sum()) for v in FEATURES}
        count_soil=int(sub[list(SOIL)].notna().all(axis=1).sum())
        bands.append({"band":label,"matched_photo_rows":int(len(sub)),
            "matched_species":int(sub.inat_taxon_id.nunique()),
            "full_10soil_rows":count_soil,
            "unique_source_site_indices":int(sub.source_index.nunique()),
            "median_source_distance_km":float(sub.distance_to_source_km.median()) if len(sub) else None,
            "per_variable_photo_support":byvar})
        species_rows=all_species_effects(sub) if len(sub) else pd.DataFrame()
        for block,vs in BLOCKS.items():
            for var in vs:
                eff=(species_rows.loc[species_rows.variable.eq(var)].copy() if len(species_rows)
                     else pd.DataFrame())
                effects.append((label,var,block,eff,stats_effects(eff,var,block,label)))
    outcome=[]; species=[]
    for band,var,block,eff,summ in effects:
        if not eff.empty:
            eff["transfer_band"]=band
            species.append(eff)
        outcome.append(summ)
    table=pd.DataFrame(outcome)
    for band in ["exact",*[f"{r:g}km" for r in radii]]:
        idx=table.transfer_band.eq(band)
        table.loc[idx,"holm_36_p"]=holm_adjust(pd.to_numeric(table.loc[idx,"nominal_p"],errors="coerce").to_numpy(float))
        table.loc[idx,"holm_36_supported"]=(table.loc[idx,"status"].eq("exploratory") &
                                              table.loc[idx,"holm_36_p"].lt(.05))
        for kind in ("observer_paired", "same_cell"):
            col=f"{kind}_nominal_p"
            table.loc[idx,f"{kind}_holm36_p"] = holm_adjust(pd.to_numeric(table.loc[idx,col],errors="coerce").to_numpy(float))
            table.loc[idx,f"{kind}_holm36_supported"] = table.loc[idx,f"{kind}_holm36_p"].lt(.05) & table.loc[idx,col].notna()
    outdir.mkdir(parents=True,exist_ok=True)
    table.to_csv(outdir/"mois_36_environment_transfer_screen.csv",index=False)
    if species:
        pd.concat(species,ignore_index=True).to_csv(outdir/"mois_36_environment_transfer_species.csv",index=False)
    receipt={"status":"EXPLORATORY_NEAREST_SOURCE_SOIL_NOT_EXACT_PHOTO_SOIL",
       "source_sites":int(len(source)),"photo_rows_after_frozen_quality_gate":int(len(photo)),
       "original_eligible_species":int(photo.inat_taxon_id.nunique()),
       "exact_photo_matches":int(mapped.source_match_exact_photo.sum()),
       "source_accuracy_max_m":max_accuracy_m,
       "target_accuracy_max_m":max_accuracy_m,
       "target_geo_accurate_rows":int((photo.latitude.between(-90,90) & photo.longitude.between(-180,180) & photo.positional_accuracy_m.between(0,max_accuracy_m)).sum()),
       "bands":bands,
       "num_soil_variables":len(SOIL),"num_total_variables":len(FEATURES),
       "tests_per_band":len(FEATURES),
       "holm_exact_supported":int(table.loc[table.transfer_band.eq("exact"),"holm_36_supported"].sum()),
       "holm_1km_supported":int(table.loc[table.transfer_band.eq("1km"),"holm_36_supported"].sum()) if "1km" in table.transfer_band.values else None,
       "holm_5km_supported":int(table.loc[table.transfer_band.eq("5km"),"holm_36_supported"].sum()) if "5km" in table.transfer_band.values else None,
       "observer_paired_all36_supported_by_band":{b:int(z["observer_paired_holm36_supported"].sum()) for b,z in table.groupby("transfer_band")},
       "same_cell_all36_supported_by_band":{b:int(z["same_cell_holm36_supported"].sum()) for b,z in table.groupby("transfer_band")},
       "source_environment_not_photo_site_exact_for_neighbours":True,
       "exact_photo_soil_5plus_species":int(pd.concat(species).loc[lambda d:d.transfer_band.eq("exact") & d.variable.eq("soil_pH")].shape[0]) if species else 0,
       "selected_nearest_site_before_checking_soil_completeness":True,
       "original_frozen_environment_mechanism_modified":False,
       "within_species_effects_have_new_posthoc_proximity_selection":True,
       "limitations":["5-km modelled SoilGrids inputs are nearest-other-source-site proxies, not sampled at each FCP high-depth photograph",
                      "observer and fine-spatial within-species effects are not identified by geographic proxies",
                      "36-variable exploratory familywise correction does not adjust for choosing the 1-km/5-km cutoff after seeing earlier FCP results"]}
    (outdir/"mois_36_environment_transfer_receipt.json").write_text(json.dumps(receipt,indent=2,ensure_ascii=False)+"\n")
    return receipt


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source-archive",type=Path,required=True)
    p.add_argument("--prospective",type=Path,required=True)
    p.add_argument("--technical",type=Path,required=True)
    p.add_argument("--join-key",type=Path,required=True)
    p.add_argument("--high-clip",type=Path,required=True)
    p.add_argument("--outdir",type=Path,required=True)
    a=p.parse_args(argv)
    source=source_table(a.source_archive)
    photos,gate=photo_table(a.prospective,a.technical,a.join_key,a.high_clip)
    receipt=run_analysis(photos,source,a.outdir)
    receipt["input_sha256"]={"source_zip":sha256(a.source_archive),"photo_csv":sha256(a.prospective),
        "technical_table":sha256(a.technical),"join_key":sha256(a.join_key),"high_clip":sha256(a.high_clip)}
    receipt["frozen_panel_reproduced"] = gate
    (a.outdir/"mois_36_environment_transfer_receipt.json").write_text(json.dumps(receipt,indent=2,ensure_ascii=False)+"\n")
    print(json.dumps({k:receipt[k] for k in ("status","exact_photo_matches","bands","holm_exact_supported","holm_1km_supported","holm_5km_supported")},ensure_ascii=False))
    return 0

if __name__=="__main__":
    raise SystemExit(main())
