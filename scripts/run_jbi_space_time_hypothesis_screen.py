#!/usr/bin/env python3
"""Exploratory space-time hypothesis screen for strict v6 FCP organization states.

Builds literature-motivated spatial/temporal niche diagnostics from occupied climate cells,
then compares C-only, S-only, and C+S states with the same documentation-IPW / family-
bootstrap multinomial design used by the rebuilt FCP analysis.

This is an exploratory screen, not confirmatory inference.
"""
from __future__ import annotations

import argparse, json
from pathlib import Path
import numpy as np
import pandas as pd
import statsmodels.api as sm
from numpy.linalg import pinv, slogdet
from scipy.spatial import ConvexHull, QhullError
from sklearn.cluster import DBSCAN, KMeans
from sklearn.linear_model import LogisticRegression

CLASSES = ["local_coexistence_only", "spatial_segregation_only", "coexistence_and_segregation"]
R_EARTH_KM = 6371.0088


def z(v: pd.Series) -> pd.Series:
    x = pd.to_numeric(v, errors="coerce").astype(float)
    sd = float(x.std(ddof=0))
    return (x - float(x.mean())) / sd if np.isfinite(sd) and sd > 0 else pd.Series(np.nan, index=x.index)


def haversine(lat1, lon1, lat2, lon2):
    p1 = np.radians(lat1); p2 = np.radians(lat2)
    dp = np.radians(lat2 - lat1); dl = np.radians(lon2 - lon1)
    a = np.sin(dp / 2) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(dl / 2) ** 2
    return 2 * R_EARTH_KM * np.arcsin(np.sqrt(np.clip(a, 0, 1)))


def spherical_xyz(lat, lon):
    lat = np.radians(lat); lon = np.radians(lon)
    return np.c_[np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)]


def bhattacharyya_overlap(x1: np.ndarray, x2: np.ndarray, eps: float = 1e-4) -> float:
    d = x1.shape[1]
    mu1, mu2 = x1.mean(0), x2.mean(0)
    s1 = np.cov(x1, rowvar=False) + eps * np.eye(d)
    s2 = np.cov(x2, rowvar=False) + eps * np.eye(d)
    s = (s1 + s2) / 2
    delta = mu1 - mu2
    sg, ld = slogdet(s); sg1, ld1 = slogdet(s1); sg2, ld2 = slogdet(s2)
    if min(sg, sg1, sg2) <= 0:
        return np.nan
    db = 0.125 * delta @ pinv(s) @ delta + 0.5 * (ld - 0.5 * (ld1 + ld2))
    return float(np.exp(-db))


def build_species_metrics(cells: pd.DataFrame, core: pd.DataFrame, seed: int, hv_draws: int) -> pd.DataFrame:
    cells = cells[cells.canonical_name.isin(core.canonical_name)].copy()
    state = core.set_index("canonical_name").organization_state.to_dict()
    rng_pairs = np.random.default_rng(seed)
    rng_hv = np.random.default_rng(seed)
    rows = []
    for sp, g in cells.groupby("canonical_name"):
        g = g.dropna(subset=["decimalLatitude", "decimalLongitude", "pc1", "pc2", "pc3"]).copy()
        n = len(g)
        if n < 20:
            continue
        lat = g.decimalLatitude.to_numpy(); lon = g.decimalLongitude.to_numpy()
        x = g[["pc1", "pc2", "pc3"]].to_numpy(float)
        mean_bio1 = float(g.bio1.mean())
        mean_bio4 = float(g.bio4.mean())
        mean_bio15 = float(g.bio15.mean())
        median_abs_lat = float(np.nanmedian(np.abs(lat)))

        coords = np.radians(g[["decimalLatitude", "decimalLongitude"]].to_numpy())
        frag = []
        for km in (100, 250, 500):
            labels = DBSCAN(eps=km / R_EARTH_KM, min_samples=1, metric="haversine", algorithm="ball_tree").fit_predict(coords)
            counts = np.bincount(labels)
            frag.append(1 - counts.max() / len(labels))
        frag_multiscale = float(np.mean(frag))

        total = n * (n - 1) // 2
        if total <= 20000:
            ii, jj = np.triu_indices(n, 1)
        else:
            ii = rng_pairs.integers(0, n, size=25000); jj = rng_pairs.integers(0, n, size=25000)
            ok = ii != jj; ii = ii[ok]; jj = jj[ok]
        gd = haversine(lat[ii], lon[ii], lat[jj], lon[jj])
        ed = np.linalg.norm(x[ii] - x[jj], axis=1)
        q25, q75 = np.quantile(gd, [0.25, 0.75])
        short = np.median(ed[gd <= q25]); long = np.median(ed[gd >= q75])
        spatial_turnover = float((long - short) / (np.median(ed) + 1e-9))

        vols = []
        for _ in range(hv_draws):
            idx = rng_hv.choice(n, 20, replace=False)
            try:
                vols.append(ConvexHull(x[idx]).volume)
            except QhullError:
                pass
        hv3d = float(np.median(vols)) if vols else np.nan

        regional_centroid = np.nan; regional_overlap = np.nan
        if n >= 40:
            labels = KMeans(n_clusters=2, random_state=seed, n_init=20).fit_predict(spherical_xyz(lat, lon))
            counts = np.bincount(labels)
            if len(counts) == 2 and counts.min() >= 15:
                x1, x2 = x[labels == 0], x[labels == 1]
                centroid = np.linalg.norm(x1.mean(0) - x2.mean(0))
                rms = np.sqrt(np.mean(np.sum((x - x.mean(0)) ** 2, axis=1)))
                regional_centroid = float(centroid / (rms + 1e-9))
                regional_overlap = bhattacharyya_overlap(x1, x2)

        rows.append({
            "canonical_name": sp,
            "organization_state": state.get(sp),
            "n_climate_cells_rebuilt": n,
            "mean_bio1": mean_bio1,
            "mean_bio4": mean_bio4,
            "mean_bio15": mean_bio15,
            "median_abs_lat": median_abs_lat,
            "frag_multiscale": frag_multiscale,
            "spatial_niche_turnover": spatial_turnover,
            "hv3d_rarefied20": hv3d,
            "regional_centroid_sep": regional_centroid,
            "regional_gaussian_overlap": regional_overlap,
        })
    return pd.DataFrame(rows)


def add_documentation_weights(d: pd.DataFrame) -> pd.DataFrame:
    d = d.copy()
    d["D_documented"] = ((d.C_local_coexistence_documented > 0) | (d.S_spatial_segregation_documented > 0)).astype(int)
    d["log_sources_z"] = z(np.log1p(d.n_FCP_eligible_sources))
    d["attention_z"] = z(np.log1p(d.n_v22_exact_name_records.fillna(0)))
    d["geo_z"] = z(np.log1p(d.geographic_radius_95_km))
    q = d.dropna(subset=["D_documented", "log_sources_z", "attention_z"]).copy()
    X = sm.add_constant(q[["log_sources_z", "attention_z"]], has_constant="add")
    fit = sm.GLM(q.D_documented.astype(int), X, family=sm.families.Binomial()).fit()
    p = np.clip(np.asarray(fit.predict(X), float), 0.05, 0.95)
    prev = float(q.D_documented.mean())
    q["stabilized_ipw"] = np.where(q.D_documented.eq(1), prev / p, (1 - prev) / (1 - p))
    return d.merge(q[["canonical_name", "stabilized_ipw"]], on="canonical_name", how="left", validate="one_to_one")


def contrasts(model: LogisticRegression) -> np.ndarray:
    c = {k: float(model.coef_[list(model.classes_).index(k), 0]) for k in CLASSES}
    return np.array([
        c["spatial_segregation_only"] - c["local_coexistence_only"],
        c["coexistence_and_segregation"] - c["local_coexistence_only"],
        c["coexistence_and_segregation"] - c["spatial_segregation_only"],
    ])


def fit_metric(d: pd.DataFrame, metric: str, controls: list[str], nboot: int, seed: int) -> dict:
    x = d[d.organization_state.isin(CLASSES)].copy()
    x["metric_z"] = z(x[metric])
    ccols = []
    for c in controls:
        name = "z_" + c
        x[name] = z(x[c]); ccols.append(name)
    cols = ["metric_z", "geo_z"] + ccols
    x = x.dropna(subset=cols + ["stabilized_ipw", "family"])
    model = LogisticRegression(C=1.0, solver="lbfgs", max_iter=2000)
    model.fit(x[cols].to_numpy(float), x.organization_state.astype(str).to_numpy(), sample_weight=x.stabilized_ipw.to_numpy(float))
    point = contrasts(model)
    families = sorted(x.family.astype(str).unique())
    rng = np.random.default_rng(seed)
    boots = []
    for _ in range(nboot):
        draws = rng.choice(families, size=len(families), replace=True)
        idx = []
        for f in draws:
            idx.extend(x.index[x.family.astype(str).eq(str(f))].tolist())
        q = x.loc[idx]
        if set(q.organization_state) != set(CLASSES):
            continue
        try:
            m = LogisticRegression(C=1.0, solver="lbfgs", max_iter=2000)
            m.fit(q[cols].to_numpy(float), q.organization_state.astype(str).to_numpy(), sample_weight=q.stabilized_ipw.to_numpy(float))
            boots.append(contrasts(m))
        except Exception:
            pass
    b = np.asarray(boots, float)
    out = {"metric": metric, "controls": "+".join(controls) if controls else "none", "n_species": len(x), "n_families": len(families), "valid_family_bootstraps": len(b)}
    for j, name in enumerate(["S_vs_C", "mixed_vs_C", "mixed_vs_S"]):
        out[name + "_OR"] = float(np.exp(point[j]))
        out[name + "_ci_low"] = float(np.exp(np.quantile(b[:, j], 0.025)))
        out[name + "_ci_high"] = float(np.exp(np.quantile(b[:, j], 0.975)))
    return out


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--analysis", required=True)
    p.add_argument("--attention", required=True)
    p.add_argument("--occupied-cells", required=True)
    p.add_argument("--core-species", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--family-bootstraps", type=int, default=499)
    p.add_argument("--hypervolume-draws", type=int, default=99)
    p.add_argument("--seed", type=int, default=20260827)
    a = p.parse_args()

    analysis = pd.read_csv(a.analysis)
    attention = pd.read_csv(a.attention)
    cells = pd.read_csv(a.occupied_cells)
    core = pd.read_csv(a.core_species)
    if len(core) != 74:
        raise SystemExit(f"Expected display-core-v6 with 74 species; got {len(core)}")

    metrics = build_species_metrics(cells, core, a.seed, a.hypervolume_draws)
    d = analysis.merge(core[["canonical_name", "organization_state", "C_local_coexistence_documented", "S_spatial_segregation_documented", "n_FCP_eligible_sources"]], on="canonical_name", how="inner", suffixes=("_analysis", ""), validate="one_to_one")
    for c in ["organization_state", "C_local_coexistence_documented", "S_spatial_segregation_documented", "n_FCP_eligible_sources"]:
        if c + "_analysis" in d.columns:
            d.drop(columns=[c + "_analysis"], inplace=True)
    d = d.merge(attention[["canonical_name", "n_v22_exact_name_records"]], on="canonical_name", how="left", validate="one_to_one")
    d = d.merge(metrics, on=["canonical_name", "organization_state"], how="inner", validate="one_to_one")
    d = add_documentation_weights(d)

    hypotheses = {
        "mean_bio4": ("H1a", "temporal_climatic_heterogeneity", "higher temperature seasonality -> C relative to S"),
        "mean_bio15": ("H1b", "temporal_climatic_heterogeneity", "higher precipitation seasonality -> C relative to S"),
        "frag_multiscale": ("H2", "spatial_fragmentation", "greater fragmentation -> S or mixed relative to C"),
        "spatial_niche_turnover": ("H3a", "spatial_environmental_turnover", "greater geographic structuring of climate space -> S"),
        "regional_centroid_sep": ("H3b", "regional_niche_separation", "greater climate centroid separation between major geographic sectors -> S"),
        "regional_gaussian_overlap": ("H3c", "regional_niche_overlap", "lower regional climate-space overlap -> S"),
        "hv3d_rarefied20": ("H0", "total_niche_hypervolume", "total climatic hypervolume size alone does not distinguish C/S organization"),
    }
    rows = []
    for metric, meta in hypotheses.items():
        r = fit_metric(d, metric, [], a.family_bootstraps, a.seed)
        r.update({"hypothesis_id": meta[0], "mechanism": meta[1], "prediction": meta[2]})
        rows.append(r)
    results = pd.DataFrame(rows)
    front = ["hypothesis_id", "mechanism", "prediction", "metric", "controls"]
    results = results[front + [c for c in results.columns if c not in front]]

    robustness = pd.DataFrame([
        fit_metric(d, "mean_bio4", controls, a.family_bootstraps, a.seed)
        for controls in [[], ["mean_bio1"], ["median_abs_lat"], ["mean_bio1", "median_abs_lat"]]
    ])

    out = Path(a.outdir); out.mkdir(parents=True, exist_ok=True)
    d.to_csv(out / "jbi_space_time_species_metrics_v1.csv", index=False)
    results.to_csv(out / "jbi_space_time_hypothesis_screen_v1.csv", index=False)
    robustness.to_csv(out / "jbi_bio4_robustness_v1.csv", index=False)
    manifest = {
        "status": "complete",
        "protocol": "space-time-hypothesis-screen-v1",
        "core_n": int(len(core)),
        "climate_eligible_n": int(len(d)),
        "informative_cs_n": int(d.organization_state.isin(CLASSES).sum()),
        "seed": a.seed,
        "family_bootstraps": a.family_bootstraps,
        "hypervolume_draws": a.hypervolume_draws,
        "conditional_model": "documentation-IPW L2 multinomial; metric + geographic_radius_95; family bootstrap",
        "hypervolume": "PC1-PC3 convex-hull volume rarefied to 20 cells; regional diagnostics use two spherical-KMeans geographic sectors and Gaussian Bhattacharyya overlap",
        "boundary": "exploratory literature-driven screen; BIO4 climatological seasonality is not year-specific temporal climate",
    }
    (out / "jbi_space_time_hypothesis_manifest_v1.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({"core_n": len(core), "analysis_n": len(d), "informative_n": int(d.organization_state.isin(CLASSES).sum()), "result_rows": len(results)}, indent=2))


if __name__ == "__main__":
    main()
