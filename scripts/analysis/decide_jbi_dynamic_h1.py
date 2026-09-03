#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

DIRECT = [
    "seasonal_centroid_cycle_mean",
    "interannual_centroid_drift_mean",
    "interannual_overlap_loss_mean",
    "annual_hypervolume_log_sd",
    "temporal_variance_component",
]
DIAGNOSTIC = ["spatial_variance_component", "space_time_variance_ratio"]


def finite(x) -> bool:
    try:
        return bool(np.isfinite(float(x)))
    except Exception:
        return False


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--models", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--reassembly-qc", required=True)
    a = p.parse_args()

    qc = json.loads(Path(a.reassembly_qc).read_text(encoding="utf-8"))
    d = pd.read_csv(a.models)

    result = {
        "protocol": "jbi-dynamic-h1-decision-v2",
        "decision_rule_source": "PR #17 comment 5506045557, frozen before dynamic model results",
        "scientific_definition_changed_from_v1": False,
        "direct_temporal_metrics": DIRECT,
        "diagnostic_S_metrics": DIAGNOSTIC,
        "reassembly_status": qc.get("status"),
        "decision": "not_evaluable",
        "supported_metrics": [],
        "metric_checks": [],
    }

    required_cols = {
        "metric",
        "status",
        "pure_S_vs_C_OR",
        "pure_S_vs_C_cluster_p_holm",
        "S_vs_C_OR_ci_high",
    }
    if qc.get("status") != "complete" or not required_cols <= set(d.columns):
        Path(a.out).write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
        print(json.dumps(result, indent=2))
        return

    if set(d.metric) != set(DIRECT + DIAGNOSTIC) or len(d) != 7 or set(d.status) != {"complete"}:
        Path(a.out).write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
        print(json.dumps(result, indent=2))
        return

    direct_evaluable = True
    for metric in DIRECT:
        r = d.loc[d.metric.eq(metric)].iloc[0]
        vals = {
            "metric": metric,
            "pure_S_vs_C_OR": None if not finite(r.pure_S_vs_C_OR) else float(r.pure_S_vs_C_OR),
            "pure_S_vs_C_cluster_p_holm": None
            if not finite(r.pure_S_vs_C_cluster_p_holm)
            else float(r.pure_S_vs_C_cluster_p_holm),
            "S_vs_C_OR_ci_high": None if not finite(r.S_vs_C_OR_ci_high) else float(r.S_vs_C_OR_ci_high),
        }
        evaluable = all(vals[k] is not None for k in ["pure_S_vs_C_OR", "pure_S_vs_C_cluster_p_holm", "S_vs_C_OR_ci_high"])
        direct_evaluable &= evaluable
        supports = bool(
            evaluable
            and vals["pure_S_vs_C_OR"] < 1.0
            and vals["pure_S_vs_C_cluster_p_holm"] <= 0.05
            and vals["S_vs_C_OR_ci_high"] < 1.0
        )
        vals["evaluable"] = evaluable
        vals["satisfies_frozen_support_conjunction"] = supports
        result["metric_checks"].append(vals)
        if supports:
            result["supported_metrics"].append(metric)

    if not direct_evaluable:
        result["decision"] = "not_evaluable"
    elif result["supported_metrics"]:
        result["decision"] = "supported"
    else:
        result["decision"] = "unsupported"

    result["claim_ceiling"] = {
        "supported": "comparative concordance between documented C/S organization and dynamic temporal environmental heterogeneity; not causal temporal selection or a temperature-specific mechanism",
        "unsupported": "the seven frozen TerraClimate space-time metrics did not discriminate C/S organization under this design; not evidence that temporal or fine-scale environmental selection is absent",
        "not_evaluable": "no biological decision; incomplete exact input, model or prespecified estimable output",
    }[result["decision"]]

    Path(a.out).write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
