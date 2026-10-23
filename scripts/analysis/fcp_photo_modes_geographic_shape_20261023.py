#!/usr/bin/env python3
"""FCP photographic phenotype modes vs broad-cline/local-nonlinear geography.

No hard colour classes enter mode fitting. K=1..4 Gaussian components are
compared with an explicit unimodal skew-normal competitor. Components are not
genetic morphs. Spatial signatures are out-of-sample predictions, not proof
of individual populations, selection, clines, or mosaics.
"""
from __future__ import annotations
import argparse
import json
import warnings
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy.signal import find_peaks
from scipy.stats import skewnorm, norm
from sklearn.cluster import KMeans
from sklearn.linear_model import Ridge
from sklearn.mixture import GaussianMixture

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "analysis"))
from fcp_full4475_hypothesis_census_20261021 import (
    load_original, SOURCE_SHA, PHOTO_LABELS
)

RADIUS_KM = 50.0
EARTH = 6371.0088
SEED = 20261023
THRESHOLD = 10.0
OBS_BOOT = 12


def check(test, text):
    if not test:
        raise RuntimeError("FCP_PHOTOGRAPHIC_MODE_GEOGRAPHY: " + text)


def sphere(lat, lon):
    a = np.deg2rad(np.asarray(lat, float))
    b = np.deg2rad(np.asarray(lon, float))
    return np.stack((np.cos(a)*np.cos(b), np.cos(a)*np.sin(b), np.sin(a)), axis=1)


def distance_km(x, y):
    dot = np.clip(np.asarray(x) @ np.asarray(y).T, -1, 1)
    return EARTH * np.arccos(dot)


def geometry(g):
    valid = (g.latitude.between(-90, 90) & g.longitude.between(-180, 180)
             & g.positional_accuracy_m.between(0, 5000))
    a = g.loc[valid].sort_values("photo_id", kind="stable").reset_index(drop=True)
    if len(a) < 2:
        return None
    xyz = sphere(a.latitude, a.longitude)
    pairs = np.triu_indices(len(a), k=1)
    geo = distance_km(xyz, xyz)[pairs]
    is_local = geo <= RADIUS_KM
    if np.count_nonzero(is_local) < 30:
        return None
    labels = pd.Categorical(a.morph, categories=PHOTO_LABELS).codes
    check(np.all(labels >= 0), "invalid frozen photographic coarse label")
    d = float(np.mean(labels[pairs[0]] != labels[pairs[1]])
              - np.mean((labels[pairs[0]] != labels[pairs[1]])[is_local]))
    hues = np.stack(a.palette_norm9.to_numpy()).astype(float)
    check(hues.shape == (len(a), 9) and np.isfinite(hues).all(),
          "invalid original measured nine-palette")
    Y = np.sqrt(hues / hues.sum(axis=1)[:, None])
    month = a.date.dt.month.fillna(0).astype(int).to_numpy()
    observers = a.observer_id.fillna("").astype(str).to_numpy()
    return a, xyz, Y, month, observers, d


def mode_scores(Y):
    z = Y - np.mean(Y, axis=0)
    u, sing, vh = np.linalg.svd(z, full_matrices=False)
    denom = float(np.sum(sing**2))
    if denom < 1e-9:
        return None, 0.
    component = z @ vh[0]
    sd = float(component.std())
    if sd < 1e-8:
        return None, 0.
    return component / sd, float(sing[0]**2 / denom)


def gaussian_peaks(model, z):
    a, b = min(z.min()-1, -4.), max(z.max()+1, 4.)
    grid = np.linspace(a, b, 512)
    loc = model.means_.reshape(-1)
    sd = np.sqrt(model.covariances_.reshape(-1))
    density = np.sum(
        model.weights_[:, None] * norm.pdf(grid[None, :], loc[:, None], sd[:, None]),
        axis=0
    )
    peaks, _ = find_peaks(density, prominence=float(density.max())*0.1, distance=15)
    return int(len(peaks))


def fit_mode(z, observers, n_boot=OBS_BOOT):
    """Conservative photographic density-mode candidate; not inherited morph K."""
    out = {"mode_status": "unresolved", "bic_selected_components": None,
           "density_peaks": None, "delta_bic_vs_unimodal": None,
           "stable_observer_subsamples": None, "mode_cluster_min_photos": None}
    if z is None or len(z) < 40 or np.std(z) < 0.1:
        return out
    def fit_once(x, observer=None, seed=0, bootstrap=False):
        x = np.asarray(x, float)
        vals = x[:, None]
        gm = [GaussianMixture(n_components=k, covariance_type="full",
                              random_state=seed, reg_covar=.02,
                              n_init=2, max_iter=120).fit(vals)
              for k in (1, 2, 3, 4)]
        bic = np.array([m.bic(vals) for m in gm])
        chosen = int(np.argmin(bic))
        fit = gm[chosen]
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                shape, loc, scale = skewnorm.fit(x)
                skew_bic = float(-2 * np.sum(skewnorm.logpdf(x, shape, loc, scale))
                                 + 3 * np.log(len(x)))
            if not np.isfinite(skew_bic):
                skew_bic = np.inf
        except (ValueError, FloatingPointError):
            skew_bic = np.inf
        margin = min(float(bic[0]), skew_bic) - float(bic[chosen])
        peaks = gaussian_peaks(fit, x)
        hard = fit.predict(vals)
        sizes = np.bincount(hard, minlength=chosen+1)
        supports = all(
            sizes[k] >= 5 and
            (observer is None or len(np.unique(observer[hard == k])) >= 3)
            for k in range(chosen+1)
        )
        candidate = (chosen >= 1 and peaks >= 2 and margin >= THRESHOLD
                     and supports and np.min(fit.weights_) >= .10)
        return candidate, chosen+1, peaks, margin, int(sizes.min())

    candidate, k, peaks, margin, min_size = fit_once(z, observers, SEED)
    out.update(bic_selected_components=k, density_peaks=peaks,
               delta_bic_vs_unimodal=margin,
               mode_cluster_min_photos=min_size)
    if not candidate:
        out["mode_status"] = "not_resolved_as_multimodal"
        return out
    if n_boot <= 0:
        out["mode_status"] = "candidate_unverified"
        return out
    rng = np.random.default_rng(SEED + int(len(z)*13))
    idx_by_observer = {}
    for j, o in enumerate(observers):
        idx_by_observer.setdefault(o, []).append(j)
    units = list(idx_by_observer.values())
    stable = 0
    for b in range(n_boot):
        keep = rng.choice(len(units), size=max(1, int(.8*len(units))),
                          replace=False)
        idx = np.array([i for k2 in keep for i in units[k2]], int)
        if len(idx) < 35:
            continue
        flag, *_ = fit_once(z[idx], observers[idx], SEED+b+1)
        stable += int(flag)
    out["stable_observer_subsamples"] = stable
    out["mode_status"] = ("recurrent_photographic_mode_candidate"
                          if stable >= int(np.ceil(.7*n_boot))
                          else "unstable_photo_mode_candidate")
    return out


def observer_folds(obs, n):
    if any(x == "" for x in obs):
        return []
    unique = np.unique(obs)
    if len(unique) < 30:
        return []
    rng = np.random.default_rng(SEED+n)
    unique = unique[rng.permutation(len(unique))]
    groups = {v: j % 3 for j, v in enumerate(unique)}
    code = np.array([groups[v] for v in obs])
    return [np.flatnonzero(code == j) for j in range(3)]


def geo_folds(xyz, n):
    # No geography-specific training/validation after sparse or tiny coverage.
    g = distance_km(xyz, xyz)
    if float(g.max()) < 100.:
        return []
    ids = KMeans(n_clusters=3, random_state=SEED+n,
                 n_init=4, max_iter=80).fit_predict(xyz)
    folds = [np.flatnonzero(ids == k) for k in range(3)]
    return folds if all(len(v) >= 8 and n-len(v) >= 30 for v in folds) else []


def model_predictions(xyz, Y, month, train, test):
    tr, te = xyz[train], xyz[test]
    mt = np.eye(13)[month[train], 1:]
    mv = np.eye(13)[month[test], 1:]
    base = Ridge(alpha=3., fit_intercept=True).fit(mt, Y[train])
    p0 = base.predict(mv)
    center = tr.mean(axis=0)
    scale = np.sqrt(np.sum((tr-center)**2, axis=1).mean())
    scale = max(scale, 1e-5)
    xt = np.concatenate([mt, (tr-center)/scale], axis=1)
    xv = np.concatenate([mv, (te-center)/scale], axis=1)
    linear = Ridge(alpha=3., fit_intercept=True).fit(xt, Y[train])
    p1 = linear.predict(xv)
    residual = Y[train] - linear.predict(xt)
    span = float(np.max(distance_km(tr, tr)))
    bandwidth = max(30., .18*span)
    weights = np.exp(-.5 * (distance_km(te, tr)/bandwidth)**2)
    p2 = p1 + (weights @ residual)/(3. + weights.sum(axis=1))[:, None]
    return p0, p1, p2


def crossval(xyz, Y, month, folds):
    n = len(Y)
    if not folds:
        return None
    errors = np.zeros(3)
    for test in folds:
        train = np.setdiff1d(np.arange(n), test, assume_unique=True)
        if len(test) < 8 or len(train) < 30:
            return None
        preds = model_predictions(xyz, Y, month, train, test)
        for i, val in enumerate(preds):
            errors[i] += np.square(val - Y[test]).sum()
    if errors[0] < 1e-9:
        return None
    return {
        "baseline_sse": float(errors[0]),
        "cline_gain_vs_month": float(1-errors[1]/errors[0]),
        "local_nonlin_gain_vs_month": float(1-errors[2]/errors[0]),
        "nonlin_increment_vs_cline": float((errors[1]-errors[2])/errors[0]),
    }


def status_shape(observer, blocked):
    if observer is None or blocked is None:
        return "geography_underidentified"
    cline = (observer["cline_gain_vs_month"] > .08
             and blocked["cline_gain_vs_month"] > .08)
    local = observer["nonlin_increment_vs_cline"] > .08
    if cline and local:
        return "broad_trend_plus_local_nonlinearity_candidate"
    if cline:
        return "transferable_broad_cline_candidate"
    if local:
        return "local_nonlinearity_only_candidate"
    return "no_clear_predictive_geometry"


def run(frame, do_boot=True):
    rows = []
    coarse_depletion = []
    for j, (tid, group) in enumerate(frame.groupby("inat_taxon_id", sort=True)):
        if len(group) < 40:
            continue
        prepared = geometry(group)
        if prepared is None:
            continue
        g, xyz, Y, month, observers, d = prepared
        coarse_depletion.append(d)
        score, variance = mode_scores(Y)
        mode = fit_mode(score, observers, n_boot=OBS_BOOT if do_boot else 0)
        ocv = crossval(xyz, Y, month, observer_folds(observers, len(g)))
        bcv = crossval(xyz, Y, month, geo_folds(xyz, len(g)))
        row = {
            "inat_taxon_id": int(tid), "species": str(g.species.iloc[0]),
            "genus": str(g.species.iloc[0]).split()[0],
            "n_original_classifiable": int(len(group)),
            "n_spatial_classifiable": int(len(g)),
            "photographic_coarse_label_count": int(g.morph.nunique()),
            "photographic_pc1_explained_variance": variance,
            "coarse_50km_depletion": d,
            **mode,
            "spatial_signature": status_shape(ocv, bcv),
        }
        for pre, cv in (("observer", ocv), ("spatial_block", bcv)):
            row[pre+"_estimable"] = cv is not None
            for metric in ("cline_gain_vs_month", "local_nonlin_gain_vs_month",
                           "nonlin_increment_vs_cline"):
                row[pre+"_"+metric] = cv[metric] if cv is not None else None
        rows.append(row)
        if (j+1) % 500 == 0:
            print("FCP_PHOTO_MODE_GEO_TAXA_SCANNED", j+1, flush=True)
    if do_boot:
        check(len(rows) == 1761, "source original strict 50km eligible species changed")
        check(abs(float(np.mean(coarse_depletion))-.0176149021685411) < 1e-10,
              "historical 50km coarse depletion changed")
    d = pd.DataFrame(rows)
    result = {
        "n_species": int(len(d)),
        "n_source_registry_species": 4475,
        "photographic_mode_status": {str(k): int(v) for k, v in
                                    d.mode_status.value_counts().items()},
        "spatial_signature": {str(k): int(v) for k, v in
                              d.spatial_signature.value_counts().items()},
        "n_spatial_block_estimable": int(d.spatial_block_estimable.sum()),
        "n_observer_estimable": int(d.observer_estimable.sum()),
        "cross_tab_modes_by_geography": pd.crosstab(
            d.mode_status, d.spatial_signature).to_dict(),
        "confirmed_old_p2_strict_n": len(rows),
        "confirmed_old_p2_strict_delta": float(np.mean(coarse_depletion)),
        "photographic_not_genetic_morphs": True,
        "cline_and_local_categories_are_cv_signatures_only": True,
        "not_an_independent_photo_provider": True,
        "posthoc_exploratory": True,
    }
    return result, d


def selftest():
    rng = np.random.default_rng(1066)
    x1 = rng.normal(0, 1, 120)
    x2 = np.concatenate([rng.normal(-3, .28, 60),
                         rng.normal(3, .28, 60)])
    ids = np.arange(120).astype(str)
    a = fit_mode(x1, ids, 0)
    b = fit_mode(x2, ids, 0)
    check(a["mode_status"] != "candidate_unverified",
          "ordinary unimodal normal was called multimodal")
    check(b["mode_status"] == "candidate_unverified",
          "two well separated photographic peaks not detected")

    lng = np.concatenate([np.linspace(-1.5, -.5, 30),
                          np.linspace(-.5, .5, 30), np.linspace(.5, 1.5, 30)])
    xyz = sphere(np.full(90, 35.), 135.+lng)
    mont = np.full(90, 6)
    grad = np.column_stack([.5+.25*lng, .5-.25*lng])
    grad = np.sqrt(grad/grad.sum(axis=1)[:, None])
    cv = crossval(xyz, grad, mont, geo_folds(xyz, 90))
    check(cv is not None and cv["cline_gain_vs_month"] > .05,
          "smooth large-scale colour gradient unrecognized")
    y_mosaic = np.where(lng < -.5, 1., np.where(lng < .5, -1., 1.))
    shade = np.sqrt(np.column_stack([.5+.35*y_mosaic,
                                     .5-.35*y_mosaic]))
    idx = np.array([str(j) for j in range(90)])
    local = crossval(xyz, shade, mont, observer_folds(idx, 90))
    check(local is not None and local["nonlin_increment_vs_cline"] > .05,
          "nonlinear geographic patch pattern not recognized")
    print("FCP_PHOTOGRAPHIC_MODE_AND_GEO_SHAPE_SYNTHETIC_PASS")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--photos", type=Path)
    parser.add_argument("--out", type=Path)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    selftest()
    if args.self_test:
        return
    check(args.photos and args.out, "authentic photos/out needed")
    raw, classified = load_original(args.photos)
    result, species = run(classified)
    args.out.mkdir(parents=True, exist_ok=True)
    species.to_csv(args.out/"per_species_mode_geo_metrics.csv", index=False)
    report = {
        "schema": "fcp_photographic_modes_geographic_shape_v1",
        "source_original_photos": len(raw),
        "source_classifiable_photos": len(classified),
        "original_source_sha256": SOURCE_SHA,
        "mode_method": "Hellinger-palette PC1, GMM K1..4 BIC versus skewnorm, actual density peaks >=2, 12 observer-subsamples",
        "spatial_method": "month baseline versus linear sphere-XYZ versus locally smoothed residual, observer-disjoint and three-region held-out CV",
        "candidate_not_proof": "a GMM component, density peak, broad-cline or local nonlinear signature is not an inherited population morph, cline or mosaic",
        "results": result,
    }
    (args.out/"summary.json").write_text(json.dumps(report, indent=2, allow_nan=False)+"\n")
    print("FCP_MODE_SHAPE_AUTHENTIC_RESULTS",json.dumps(result, allow_nan=False))


if __name__ == "__main__":
    main()
