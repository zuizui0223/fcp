#!/usr/bin/env python3
"""Toy dynamical conditions, NOT fitted FCP photo/temperature/genetic parameters.

True inherited flower-morph frequency p (NOT photo pixel/share):
dp/dt = p(1-p) * [s(T)-c_route + k*(1-2p)]
For k>0, local polymorphism protected by mutual rare-morph invasion iff
|s(T)-c_route| < k. Both fixation boundaries and within-population mixtures
exist for the SAME biochemical route under different selection/migration.
Actual population history/gene flow and reproductive isolation are omitted.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path


def phase(*, s:float, c_route:float, k:float)->dict:
    if k<=0:
        raise ValueError("k>0 required for stabilizing rare-morph advantage")
    net=s-c_route
    rare_A=net+k
    rare_B=k-net
    if rare_A>0 and rare_B>0:
        x=.5*(1+net/k)
        verdict="protected_intrapopulation_polymorphism"
    elif net>=k:
        x=1.
        verdict="locally_A_fixed_or_dominant"
    elif net<=-k:
        x=0.
        verdict="locally_B_fixed_or_dominant"
    else:
        raise AssertionError("unreachable")
    return {
        "s":s,"route_pleiotropy_cost":c_route,"negative_frequency_dependence":k,
        "rare_A_invasion_growth":rare_A,"rare_B_invasion_growth":rare_B,
        "stable_p_true_genetic_morph_A":x,
        "local_outcome":verdict,
        "photo_frequency_not_equivalent_to_genetic_frequency":True,
        "reproductive_isolation_not_in_model":True,
    }


def same_static_frequency_multiple_mechanisms(p:float,k:float,m:float,bar_p:float)->float:
    """Solve for an s that yields arbitrary observed p under migration+NFD.

    dp/dt = p(1-p)[s-c+k(1-2p)] + m*(bar_p-p).
    Here returns s-c, NOT causal selection inferred from observed photo p.
    """
    if not (0<p<1) or k<0 or m<0 or not 0<=bar_p<=1:
        raise ValueError("invalid true morph frequency or model parameters")
    return -k*(1-2*p)-m*(bar_p-p)/(p*(1-p))


def workflow()->dict:
    # Hypothetical illustration only: all numbers are arbitrary, not derived
    # from BIO5/fitness RNAseq comparisons and not a pollinator selection model.
    temps=(15,20,25,30,35)
    delta_s_per_C=.02
    center_C=25
    scenarios=(
        {"route":"petal_only_regulatory_loss","illustrative_c":.02},
        {"route":"whole_plant_systemic_loss","illustrative_c":.24},
    )
    s=[]
    for a in scenarios:
        states=[]
        for t in temps:
            selection=delta_s_per_C*(t-center_C)
            states.append({"illustrative_T_C":t,**phase(s=selection,
                 c_route=a["illustrative_c"],k=.20)})
        s.append({"route":a["route"],"results":states})
    return {
        "schema":"fcp_tissue_modularity_BIO5_selection_maintenance_phase_illustration_v1",
        "toy_parameter_values_not_fitted":True,
        "not_a_BIO5_association_significance_test":True,
        "not_morph_evolution_empirically_observed":True,
        "not_a_reproductive_isolation_model":True,
        "genetic_morph_A_not_necessarily_white":True,
        "frequency_dependent_local_protection_criterion":"abs(s(T)-pleiotropic_cost)<k",
        "standing_migration_can_produce_spatial_frequency_structure_without_speciation":True,
        "selection_from_environment_and_pathway_tissue_impact_are_separable_parameters":True,
        "example_parameterization":{"temp_center_C":center_C,
            "slope_per_C":delta_s_per_C,"frequency_dependence_k":.20,
            "tissue_pleiotropy_costs_illustrative":[.02,.24]},
        "route_scenarios":s,
        "unidentifiable_from_single_cross_section":
            "Any given genetic morph-frequency p in (0,1) can be stationary under multiple combinations of frequency dependence k, migration m and locally signed selection s; photographic morph frequencies contain even less information.",
    }


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--out",required=True,type=Path)
    opts=ap.parse_args()
    output=workflow()
    opts.out.parent.mkdir(parents=True,exist_ok=True)
    opts.out.write_text(json.dumps(output,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print("FCP_THEORETICAL_NONIDENTIFIABLE_MAINTENANCE_FILTER",output["schema"])


if __name__=="__main__":
    main()
