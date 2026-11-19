#!/usr/bin/env python3
"""FCP: diagnose estimand-specific attrition, structural overlap and selection.

Source: frozen 4,475-species ledger, Table S10 and Table S2.
Output: a result JSON and 4,475-row audit table; NOT an original-photo refit.
Selection correction is explicitly a sensitivity restricted to *observable*
geographic, repeated-colour, sufficient-photo opportunities (NOT all species).

python scripts/analysis/fcp4475_selection_transport_diagnostic_20261119.py \
  --out-dir /tmp/fcp4475-selection-transport
"""
from __future__ import annotations
import argparse
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import brier_score_loss, roc_auc_score
from sklearn.model_selection import GroupKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

ROOT = Path(__file__).resolve().parents[2]
MASTER = ROOT / "results/fcp_np_all4475_unified_species_ledger_20261021.csv"
COMMON = ROOT / "manuscripts/climate_colour/Table_S10_common_species.csv.gz"
S2 = ROOT / "manuscripts/climate_colour/Table_S2_species_coefficients.csv.gz"
SOURCE_SHA = "23d3daadc9799c2e7b5867d54ccc1ed98369ba79a81b72dc7f831e5197e63c49"
SEED = 20261119
REPS = 2500


def check(condition, label):
    if not bool(condition):
        raise AssertionError("FCP_TRANSPORT_DIAGNOSTIC: " + label)


def load(path):
    return pd.read_csv(path, dtype={"inat_taxon_id": "int64"}, low_memory=False)


def quant(a, p):
    return float(np.quantile(np.asarray(a, dtype=float), p))


def score_data():
    check(hashlib.sha256(MASTER.read_bytes()).hexdigest() == SOURCE_SHA,
          "Master species ledger does not match verified source SHA256")
    m, c, t = [load(p) for p in (MASTER, COMMON, S2)]
    for name, x in [("master", m), ("common", c), ("thermal", t)]:
        check(len(x) == 4475 and x.inat_taxon_id.notna().all() and
              x.inat_taxon_id.is_unique, "4475 unique IDs: " + name)
    d = m.merge(c, on="inat_taxon_id", validate="one_to_one", suffixes=("", "_common"))
    d = d.merge(t[["inat_taxon_id", "species", "source_cohort",
                   "thermal_level_white", "n_original_classifiable", "status"]],
                on="inat_taxon_id", validate="one_to_one", suffixes=("", "_S2"))
    check(len(d) == 4475 and
          (d.species == d.species_S2).all() and
          (d.source_cohort == d.source_cohort_S2).all(), "source names/cohorts drift")
    check((d.n_classified_photo == d.n_original_classifiable).all(),
          "source classified photo count changed")
    check(np.allclose(d.loc[d.geo.notna(), "geo"],
                      d.loc[d.geo.notna(), "global_vs_local_50km_colour_mismatch_delta"],
                      atol=1e-8) and
          d.geo.notna().equals(d.global_vs_local_50km_colour_mismatch_delta.notna()),
          "geographic source coefficient mismatch")
    check(np.allclose(d.loc[d.thermal.notna(), "thermal"],
                      d.loc[d.thermal.notna(), "thermal_level_white"],
                      atol=1e-8) and d.thermal.notna().equals(d.thermal_level_white.notna()),
          "shared thermal source coefficient mismatch")
    d["R"] = d.repeated.astype(bool)
    d["G"] = d.geo.notna()
    d["T"] = d.thermal.notna()
    d["C"] = d.common.astype(bool)
    check(((d.R & d.G & d["T"]) == d.C).all(), "common equals explicit intersection")
    check((int(d.R.sum()), int(d.G.sum()), int(d["T"].sum()), int(d.C.sum())) ==
          (2400, 3600, 2837, 1070), "denominator drift")
    check((d.n_local50km_colour_photo_pairs.gt(0) == d.G).all(),
          "geographic coefficient unexpectedly absent when 50km pair exists")
    d["G_R"] = d.G & d.R
    d["T_R"] = d["T"] & d.R
    # Hard eligibility for geo-to-thermal selection, outcome-blind opportunity only.
    d["G_H"] = d.G_R & d.n_classified_photo.ge(40) & d.n_distinct_observers_source.ge(20)
    check((d.C <= d.G_H).all(), "common taxa outside hard geographic support")
    d["T_near"] = d.T_R & d.n_local50km_colour_photo_pairs.gt(0)
    check((d.T_near == d.C).all(), "near-pair structural gate changed")
    return d


def summary(d):
    out = {
        "registry": len(d), "G": int(d.G.sum()), "T": int(d["T"].sum()),
        "R": int(d.R.sum()), "G_R": int(d.G_R.sum()),
        "T_R": int(d.T_R.sum()), "C": int(d.C.sum()),
        "G_H": int(d.G_H.sum()), "T_near": int(d.T_near.sum()),
        "G_R_below_photo40_or_source_observers20":
            int((d.G_R & ~d.G_H).sum()),
        "T_R_no_50km_pair":
            int((d.T_R & ~d.T_near).sum()),
    }
    dsets = {
        "geography_natural": d.loc[d.G, "geo"],
        "geography_repeated": d.loc[d.G_R, "geo"],
        "geography_hard_opportunity": d.loc[d.G_H, "geo"],
        "geography_paired": d.loc[d.C, "geo"],
        "thermal_natural": d.loc[d["T"], "thermal"],
        "thermal_repeated": d.loc[d.T_R, "thermal"],
        "thermal_paired": d.loc[d.C, "thermal"],
    }
    out["means"] = {k: float(s.mean()) for k, s in dsets.items()}
    out["source_selection_coverage"] = {
        "geography_paired_fraction_of_natural": float(d.C.sum() / d.G.sum()),
        "thermal_paired_fraction_of_natural": float(d.C.sum() / d["T"].sum()),
        "geography_common_fraction_after_hard_gate": float(d.C.sum() / d.G_H.sum()),
        "thermal_no_pair_fraction_of_repeated": float((d.T_R & ~d.T_near).sum() / d.T_R.sum()),
    }
    return out


def covariates(frame):
    # Sampling opportunity only; NEVER white fractions, R condition or
    # geographic/thermal fitted outcomes as covariates.
    n = len(frame)
    arr = pd.DataFrame(index=frame.index)
    for col in ("n_classified_photo", "n_distinct_observers_source",
                "n_spatial_valid_classified", "n_local50km_colour_photo_pairs",
                "n_local50km_diff_observer_pairs"):
        arr["log1p_" + col] = np.log1p(frame[col].fillna(0).clip(lower=0))
    arr["classifiability"] = (
        frame.n_classified_photo / frame.n_source_photo)
    for col in ("ob_bio5_observed_SD_C", "ob_bio5_observed_span_C"):
        x = pd.to_numeric(frame[col], errors="coerce")
        arr[col + "_missing"] = x.isna().astype(float)
        median = float(x.median()) if x.notna().any() else 0.
        arr[col] = np.log1p(x.fillna(median).clip(lower=0))
    cat = pd.get_dummies(frame.source_cohort, prefix="cohort", dtype=float)
    arr = pd.concat([arr, cat], axis=1).astype(float)
    check(len(arr) == n and np.isfinite(arr.to_numpy()).all(),
          "invalid selection predictors")
    return arr


def crossfit_propensity(d):
    # Our only possible target is G_H, not all 3,600 G-eligible species.
    # The full 4,475 species cannot be 'corrected' where G / T is undefined.
    target = d[d.G_H].copy()
    X = covariates(target).to_numpy()
    y = target.C.to_numpy(dtype=int)
    g = target.genus.to_numpy()
    check(y.sum() == 1070 and (y == 0).sum() > 0, "no mixed selection cases")
    score = np.full(len(target), np.nan)
    # A source genus may only appear in one test fold, preventing genus leakage.
    for train, test in GroupKFold(n_splits=5).split(X, y, groups=g):
        model = make_pipeline(StandardScaler(),
                              LogisticRegression(C=0.5, max_iter=2500,
                                                 solver="lbfgs"))
        model.fit(X[train], y[train])
        score[test] = model.predict_proba(X[test])[:, 1]
    check(np.isfinite(score).all(), "crossfit propensity contains NA")
    score = np.clip(score, 1e-6, 1 - 1e-6)
    target["propensity_oof"] = score
    selected = target.C.to_numpy()
    estimand = target.geo.to_numpy(dtype=float)
    result = {
        "target": "geography repeated-colour taxa with >=40 classified source photos and >=20 original observers, and <=50km colour pair",
        "n_target": len(target),
        "n_selected": int(selected.sum()),
        "folds": 5, "grouping": "genus",
        "features": ["source cohort", "classified photo count/fraction",
                     "source observer count", "spatial-valid classified count",
                     "50km local pair opportunity", "different-observer 50km pairs",
                     "observed BIO5 span and SD where available"],
        "auc_oof": float(roc_auc_score(y, score)),
        "brier_oof": float(brier_score_loss(y, score)),
        "base_rate": float(y.mean()),
        "predicted_p_quantiles": {
            str(x): quant(score, x) for x in (.01, .05, .5, .95, .99)
        },
        "selected_predicted_p_below_0_05": int((score[selected] < .05).sum()),
        "selected_predicted_p_below_0_10": int((score[selected] < .10).sum()),
    }
    calib = []
    for indexes in np.array_split(np.argsort(score, kind="stable"), 10):
        calib.append({"n": len(indexes),
                      "mean_predicted": float(score[indexes].mean()),
                      "observed_selected_fraction": float(y[indexes].mean())})
    result["decile_calibration"] = calib
    weighted = {}
    for floor in (.02, .05, .10):
        p = np.clip(score, floor, .98)
        w = np.where(selected, 1 / p, 0.)
        ess = float(w.sum() ** 2 / np.sum(w * w))
        weighted[f"floor_{floor:.2f}"] = {
            "mean": float(np.sum(w * estimand) / np.sum(w)),
            "ess_selected": ess,
            "ess_over_n_selected": ess / selected.sum(),
            "max_weight": float(np.max(w)),
            "target_observed_mean": float(np.mean(estimand)),
        }
        if floor == .05:
            target["weight"] = w
    result["weight_sensitivity"] = weighted
    # 95% CI for weighted vs actual target using paired whole-genus bootstrap,
    # holding crossfit propensity scores fixed. Does not include fit uncertainty.
    cluster = []
    for _, frame in target.groupby("genus", sort=True):
        cluster.append(np.array([frame.geo.sum(), len(frame),
                                 (frame.geo * frame.weight).sum(),
                                 frame.weight.sum()], dtype=float))
    mat = np.stack(cluster)
    rng = np.random.default_rng(SEED + 34)
    boot = []
    for _ in range(REPS):
        take = mat[rng.integers(len(mat), size=len(mat))].sum(axis=0)
        boot.append(take[2] / take[3] - take[0] / take[1])
    result["weighted_minus_observed_target"] = float(
        weighted["floor_0.05"]["mean"] - float(np.mean(estimand)))
    result["weighted_minus_observed_target_genus95_fixed_propensities"] = [
        quant(boot, .025), quant(boot, .975)]
    result["n_genera_target"] = len(mat)
    d["propensity_oof"] = np.nan
    d["weight_transport_geo"] = np.nan
    d.loc[target.index, "propensity_oof"] = target["propensity_oof"]
    d.loc[target.index, "weight_transport_geo"] = target["weight"]
    return result


def gate_bootstrap(d):
    # Samples whole genera from original 4,475-species registry and estimates
    # each nested gate difference using exact, paired bootstrap replicates.
    groups = []
    for _, f in d.groupby("genus", sort=True):
        arr = []
        for flag, key in (("G", "geo"), ("G_R", "geo"),
                          ("G_H", "geo"), ("C", "geo"),
                          ("T", "thermal"), ("T_R", "thermal"),
                          ("C", "thermal")):
            z = f.loc[f[flag], key].to_numpy(float)
            arr += [float(z.sum()), float(len(z))]
        groups.append(arr)
    a = np.array(groups)
    rng = np.random.default_rng(SEED + 15)
    outcome = {
        "geo_repeated_minus_natural": [],
        "geo_hard_minus_repeated": [],
        "geo_paired_minus_hard": [],
        "geo_paired_minus_repeated": [],
        "thermal_repeated_minus_natural": [],
        "thermal_paired_minus_repeated": [],
    }
    def evaluate(v):
        av = [v[2 * i] / v[2 * i + 1] for i in range(7)]
        return {
            "geo_repeated_minus_natural": av[1] - av[0],
            "geo_hard_minus_repeated": av[2] - av[1],
            "geo_paired_minus_hard": av[3] - av[2],
            "geo_paired_minus_repeated": av[3] - av[1],
            "thermal_repeated_minus_natural": av[5] - av[4],
            "thermal_paired_minus_repeated": av[6] - av[5],
        }
    point = evaluate(a.sum(axis=0))
    for _ in range(REPS):
        g = a[rng.integers(len(a), size=len(a))].sum(axis=0)
        for k, v in evaluate(g).items():
            outcome[k].append(v)
    return {
        "point": point, "genus95": {
            k: [quant(v, .025), quant(v, .975)] for k, v in outcome.items()},
        "n_genera": len(groups), "reps": REPS,
        "confidence_scope": "species coefficients fixed; entire genera resampled"
    }


def by_cohort(d):
    values = []
    for key, f in d.groupby("source_cohort", sort=True):
        c = f[f.C]
        values.append({
            "cohort": key, "n_registry": len(f), "n_geo": int(f.G.sum()),
            "n_thermal": int(f.T.sum()), "n_repeated": int(f.R.sum()),
            "n_common": len(c), "geo_common_mean": float(c.geo.mean()),
            "thermal_common_mean": float(c.thermal.mean()),
        })
    return values


def main(output):
    d = score_data()
    info = summary(d)
    boot = gate_bootstrap(d)
    weights = crossfit_propensity(d)
    out = {
        "schema": "fcp_4475_eligibility_selection_transport_sensitivity_20261119",
        "source": "frozen species-level original-photo estimands; not a photo-level refit",
        "source_sha256_master": SOURCE_SHA,
        "seed": SEED,
        "primary": info,
        "gate_decomposition_bootstrap": boot,
        "geography_common_to_hard_opportunity_transport": weights,
        "cohort_replication": by_cohort(d),
        "identifiability": [
            "Repeated-colour is defined by observed colour labels and therefore is an outcome-dependent inclusion gate: differences across R/no-R cannot be removed by opportunity weighting alone.",
            "Species without >=1 classified <=50km conspecific pair have no geographic effect; species without identifiable thermal contrast have no climate coefficient. Do not impute zero or extrapolate these effects.",
            "Among thermal + repeated-colour eligible species, no 50km classified pair deterministically precludes entry into the paired common set. No transport is identified to this stratum.",
            "Opportunity weights are cross-fit by genus in the geographic repeated-colour and >=40 classified / >=20 source-observer support only; no claims about entire 4,475 taxa.",
            "Propensity features include classification success, which itself can depend on photographic colour or floral organs; weighted estimates are sensitivity checks, not MAR proofs.",
            "Whole-genus intervals use fixed archived species coefficients and fixed learned weights; omit original photo-level fitting uncertainty and propensity-model refit uncertainty.",
            "Cross-validated propensities are model-dependent; in-sample observed target means are supplied as a direct calibration benchmark.",
        ]
    }
    output.mkdir(parents=True, exist_ok=True)
    (output / "fcp4475_selection_transport_diagnostic_20261119.json").write_text(
        json.dumps(out, ensure_ascii=False, indent=2) + "\n", encoding="utf8")
    fields = ["inat_taxon_id", "species", "genus", "source_cohort",
              "n_source_photo", "n_classified_photo", "n_distinct_observers_source",
              "n_local50km_colour_photo_pairs", "R", "G", "T", "C", "G_H",
              "G_R", "T_R", "T_near", "geo", "thermal",
              "propensity_oof", "weight_transport_geo"]
    d[fields].to_csv(output / "fcp4475_selection_transport_species_20261119.csv", index=False)
    print(json.dumps({
        "counts": {k: v for k, v in info.items() if isinstance(v, int)},
        "means": info["means"],
        "gate_changes": boot["point"],
        "gate_CI": boot["genus95"],
        "propensity": {k: v for k, v in weights.items()
                       if k not in ("decile_calibration", "features")},
        "cohorts": out["cohort_replication"],
    }, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--out-dir", type=Path,
                        default=ROOT / "output/fcp_eligibility_selection")
    args = parser.parse_args()
    main(args.out_dir)
