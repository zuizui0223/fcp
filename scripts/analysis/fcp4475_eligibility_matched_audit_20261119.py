#!/usr/bin/env python3
"""Source-joined 4,475-species FCP estimability ledger and paired-genus audit.

Reads frozen species-level outputs; it does not refit original photograph data.
No exact photo, observer or coordinate identifiers are exported.
Standard library only. Run from repository root:
  python scripts/analysis/fcp4475_eligibility_matched_audit_20261119.py \
      --out-dir /tmp/fcp-eligibility --check-archive
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
MASTER = ROOT / "results/fcp_np_all4475_unified_species_ledger_20261021.csv"
COMMON = ROOT / "manuscripts/climate_colour/Table_S10_common_species.csv.gz"
THERMAL = ROOT / "manuscripts/climate_colour/Table_S2_species_coefficients.csv.gz"
ARCHIVE = ROOT / "results/fcp4475_matched_eligibility_effect_audit_20261119.json"
EXPECTED_MASTER_SHA256 = "23d3daadc9799c2e7b5867d54ccc1ed98369ba79a81b72dc7f831e5197e63c49"
KEY = "inat_taxon_id"
TOTAL = 4475
SEED = 20261119
REPS = 2500


def require(test: bool, explanation: str) -> None:
    if not test:
        raise AssertionError("FCP_ELIGIBILITY_AUDIT: " + explanation)


def load(path: Path) -> list[dict[str, str]]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as stream:
        result = list(csv.DictReader(stream))
    require(len(result) == TOTAL, f"{path}: expected {TOTAL} species, got {len(result)}")
    ids = [r[KEY] for r in result]
    require(len(set(ids)) == TOTAL and all(ids), f"{path}: duplicated/missing ID")
    return result


def numeric(value: str | None) -> float | None:
    if value is None or value.strip() in ("", "NA", "nan", "NaN"):
        return None
    out = float(value)
    return out if math.isfinite(out) else None


def flag(s: str) -> bool:
    require(s in ("True", "False"), f"unexpected boolean {s!r}")
    return s == "True"


def q(values: list[float], fraction: float) -> float:
    sorted_values = sorted(values)
    f = (len(sorted_values) - 1) * fraction
    i = math.floor(f)
    j = min(i + 1, len(sorted_values) - 1)
    return sorted_values[i] + (sorted_values[j] - sorted_values[i]) * (f - i)


def means(s: list[float]) -> dict[str, float]:
    geo_all = s[0] / s[1]
    geo_repeated = s[2] / s[3]
    thermal_all = s[4] / s[5]
    thermal_repeated = s[6] / s[7]
    geo_common = s[8] / s[10]
    thermal_common = s[9] / s[10]
    rho = ((s[13] - s[8] * s[9] / s[10]) /
           math.sqrt((s[11] - s[8] ** 2 / s[10]) *
                     (s[12] - s[9] ** 2 / s[10])))
    return {
        "geo_all": geo_all, "geo_repeated": geo_repeated, "geo_common": geo_common,
        "thermal_all": thermal_all, "thermal_repeated": thermal_repeated,
        "thermal_common": thermal_common,
        "geo_common_minus_all": geo_common - geo_all,
        "thermal_common_minus_all": thermal_common - thermal_all,
        "geo_common_minus_repeated": geo_common - geo_repeated,
        "thermal_common_minus_repeated": thermal_common - thermal_repeated,
        "geo_thermal_common_pearson_unadjusted": rho,
    }


class Xorshift:
    """Match the unsigned 32-bit JS xorshift audit bootstrap exactly."""

    def __init__(self, seed: int) -> None:
        self.seed = seed & 0xFFFFFFFF

    def sample_index(self, n: int) -> int:
        s = self.seed
        s = (s ^ (s << 13)) & 0xFFFFFFFF
        s ^= s >> 17
        s = (s ^ (s << 5)) & 0xFFFFFFFF
        self.seed = s
        return (s * n) // 4294967296


def bootstrap(groups: list[list[float]]) -> dict:
    rand = Xorshift(SEED)
    sims: dict[str, list[float]] = {}
    for _ in range(REPS):
        sums = [0.0] * 14
        for _ in range(len(groups)):
            group = groups[rand.sample_index(len(groups))]
            for k in range(14):
                sums[k] += group[k]
        point = means(sums)
        for name, value in point.items():
            sims.setdefault(name, []).append(value)
    return {
        "method": "whole-genus replacement over fixed registry; paired changes in the same resamples",
        "seed": SEED, "reps": REPS, "n_genera": len(groups),
        "interval_type": "percentile_95",
        "intervals": {name: [q(samples, 0.025), q(samples, 0.975)]
                      for name, samples in sims.items()},
    }


def main(out_dir: Path, check_archive: bool) -> None:
    require(hashlib.sha256(MASTER.read_bytes()).hexdigest() == EXPECTED_MASTER_SHA256,
            "master CSV SHA-256 does not match frozen receipt")
    master, common, s2 = load(MASTER), load(COMMON), load(THERMAL)
    mi = {r[KEY]: r for r in master}
    ti = {r[KEY]: r for r in s2}
    require({r[KEY] for r in common} == set(mi) == set(ti),
            "not the same source 4475 taxon IDs")
    require(sum(int(r["n_classified_photo"]) for r in master) == 227107,
            "source classified-photo count drifted")
    counts: Counter[str] = Counter()
    failures_geo: Counter[str] = Counter()
    failures_thermal: Counter[str] = Counter()
    groups: dict[str, list[float]] = {}
    output = []
    per_cohort: dict[str, dict[str, float]] = {}
    for c in common:  # preserve Table S10 order for exact seeded genus bootstrap
        key = c[KEY]
        m, t = mi[key], ti[key]
        require(m["species"] == t["species"] and
                m["source_cohort"] == t["source_cohort"], "taxon/cohort join drift")
        require(int(m["n_classified_photo"]) == int(t["n_original_classifiable"]),
                f"classified-photo gate drift: {key}")
        g, warm = numeric(c["geo"]), numeric(c["thermal"])
        source_g = numeric(m["global_vs_local_50km_colour_mismatch_delta"])
        source_t = numeric(t["thermal_level_white"])
        require((g is None) == (source_g is None) and
                (warm is None) == (source_t is None) and
                (g is None or abs(g - source_g) <= 1e-8) and
                (warm is None or abs(warm - source_t) <= 1e-8),
                f"archive coefficient changed: {key}")
        repeated, paired = flag(c["repeated"]), flag(c["common"])
        require(paired == (repeated and g is not None and warm is not None),
                f"non-intersection common flag: {key}")
        local5 = flag(c["local5_opportunity"])
        classified = int(m["n_classified_photo"])
        original = int(m["n_source_photo"])
        near = int(m["n_local50km_colour_photo_pairs"])
        bio = numeric(m["ob_BIO5_Firth_beta_logodds_per_C"])
        reason_g = ("eligible" if g is not None else
                    "no_50km_classifiable_pair" if near == 0 else
                    "geo_coefficient_unavailable_despite_pairs")
        reason_t = ("eligible" if warm is not None else
                    "less_than_40_classifiable" if classified < 40 else
                    t["status"] if t["status"] != "source_photo_macroclimate_estimable"
                    else "white_thermal_coefficient_unavailable")
        failures_geo[reason_g] += 1
        failures_thermal[reason_t] += 1
        counts["all"] += 1
        counts["repeated"] += int(repeated)
        counts["geo"] += int(g is not None)
        counts["thermal"] += int(warm is not None)
        counts["bio5_firth"] += int(bio is not None)
        counts["geo_repeated"] += int(repeated and g is not None)
        counts["thermal_repeated"] += int(repeated and warm is not None)
        counts["common"] += int(paired)
        counts["local5_opportunity"] += int(local5)
        counts["n_classified_photos"] += classified
        cohort = m["source_cohort"]
        cb = per_cohort.setdefault(cohort, {"n": 0, "common": 0, "geo_sum": 0., "thermal_sum": 0.})
        cb["n"] += 1
        cb["common"] += int(paired)
        if paired:
            cb["geo_sum"] += g
            cb["thermal_sum"] += warm
        a = groups.setdefault(c["genus"], [0.0] * 14)
        if g is not None:
            a[0] += g
            a[1] += 1
            if repeated:
                a[2] += g
                a[3] += 1
        if warm is not None:
            a[4] += warm
            a[5] += 1
            if repeated:
                a[6] += warm
                a[7] += 1
        if paired:
            a[8] += g
            a[9] += warm
            a[10] += 1
            a[11] += g * g
            a[12] += warm * warm
            a[13] += g * warm
        output.append({
            "inat_taxon_id": key, "species": m["species"], "genus": c["genus"],
            "source_cohort": cohort,
            "n_source_photo": original, "n_classified_photo": classified,
            "n_unclassified_photo": int(m["n_missing_colour_label"]),
            "n_distinct_observer_source": int(m["n_distinct_observers_source"]),
            "classifiability_fraction": classified / original,
            "n_local50km_photo_pairs": near,
            "n_local5km_same_year_month_diff_observer_pairs":
                int(m["n_local5km_same_year_month_diff_observer_pairs"]),
            "repeated_colour": repeated, "geo_estimable": g is not None,
            "thermal_estimable": warm is not None, "bio5_firth_estimable": bio is not None,
            "geo_repeated_estimable": repeated and g is not None,
            "thermal_repeated_estimable": repeated and warm is not None,
            "paired_common": paired, "local5_opportunity": local5,
            "geo_effect": g, "thermal_white_slope": warm,
            "bio5_firth_logodds_per_C": bio,
            "geo_gate_reason": reason_g, "thermal_gate_reason": reason_t,
            "s2_fit_status": t["status"],
        })

    require(dict((k, counts[k]) for k in
                 ("all", "repeated", "geo", "thermal", "bio5_firth",
                  "geo_repeated", "thermal_repeated", "common", "local5_opportunity",
                  "n_classified_photos")) ==
            dict(zip(("all", "repeated", "geo", "thermal", "bio5_firth",
                      "geo_repeated", "thermal_repeated", "common",
                      "local5_opportunity", "n_classified_photos"),
                     (4475, 2400, 3600, 2837, 4112, 1902, 1376, 1070, 840, 227107))),
            "registered denominator/checksum drift")
    require(len(groups) == 1634, "genera count changed")
    sums = [sum(a[k] for a in groups.values()) for k in range(14)]
    point = means(sums)
    conf = bootstrap(list(groups.values()))
    frac_all = sum(r["classifiability_fraction"] for r in output) / TOTAL
    common_output = [r for r in output if r["paired_common"]]
    frac_common = sum(r["classifiability_fraction"] for r in common_output) / len(common_output)
    for cb in per_cohort.values():
        cb["geo_common_mean"] = cb.pop("geo_sum") / cb["common"]
        cb["thermal_common_mean"] = cb.pop("thermal_sum") / cb["common"]
    receipt = {
        "schema": "fcp_4475_matched_species_eligibility_refit_free_reproduction_20261119",
        "status": "source_coefficients_joined_no_original_photo_refit",
        "inputs": {"registry": str(MASTER.relative_to(ROOT)),
                   "common": str(COMMON.relative_to(ROOT)),
                   "thermal": str(THERMAL.relative_to(ROOT)),
                   "registry_sha256": EXPECTED_MASTER_SHA256},
        "n": dict(counts), "n_genera": len(groups),
        "means": point, "bootstrap": conf,
        "selection_descriptives": {
            "classified_photo_fraction_mean_all": frac_all,
            "classified_photo_fraction_mean_common": frac_common,
            "median_local50_pairs_all": q([r["n_local50km_photo_pairs"] for r in output], .5),
            "median_local50_pairs_common":
                q([r["n_local50km_photo_pairs"] for r in common_output], .5),
        },
        "cohorts": per_cohort, "eligibility_reason_geo": dict(failures_geo),
        "eligibility_reason_thermal": dict(failures_thermal),
        "limitations": [
            "Selection and reporting audit, not full propensity adjustment or transport",
            "Genus intervals hold archived per-species estimates fixed",
            "Unadjusted paired-species correlation is not a causal or mediated association",
            "Photographic label is not validated heritable floral morph",
        ],
    }
    if check_archive:
        frozen = json.loads(ARCHIVE.read_text(encoding="utf8"))
        for key in ("geo_all", "geo_repeated", "geo_common",
                    "thermal_all", "thermal_repeated", "thermal_common",
                    "geo_common_minus_all", "thermal_common_minus_all",
                    "geo_thermal_common_pearson_unadjusted"):
            require(abs(point[key] - frozen["means"][key]) < 1e-11,
                    "frozen species estimate mismatch: " + key)
            reference = frozen["bootstrap"]["intervals"][key]
            result = conf["intervals"][key]
            require(all(abs(a - b) < 1e-10 for a, b in zip(reference, result)),
                    "frozen genus bootstrap mismatch: " + key)
    out_dir.mkdir(parents=True, exist_ok=True)
    with (out_dir / "fcp4475_species_eligibility_20261119.csv").open("w", encoding="utf8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(output[0]))
        writer.writeheader()
        writer.writerows(output)
    (out_dir / "fcp4475_eligibility_bootstrap_20261119.json").write_text(
        json.dumps(receipt, indent=2, ensure_ascii=False) + "\n", encoding="utf8"
    )
    print(json.dumps({"n": dict(counts), "point": point,
                      "geo_common_minus_all_CI":
                      conf["intervals"]["geo_common_minus_all"],
                      "thermal_common_minus_all_CI":
                      conf["intervals"]["thermal_common_minus_all"],
                      "eligible_reasons_geo": dict(failures_geo),
                      "eligible_reasons_thermal": dict(failures_thermal)},
                     ensure_ascii=False, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--out-dir", type=Path, default=ROOT / "output/fcp_eligibility")
    parser.add_argument("--check-archive", action="store_true")
    args = parser.parse_args()
    main(args.out_dir, args.check_archive)
