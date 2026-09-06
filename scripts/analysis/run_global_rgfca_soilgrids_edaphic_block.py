#!/usr/bin/env python3
"""Run the fixed SoilGrids edaphic member of the expanded environmental panel.

Four 0-30 cm properties are built from the exact 5-km SoilGrids source layers
using the frozen 5/10/15-cm thickness weights.  The block reports a raw p-value
only; final support remains pending Holm correction across the full fixed panel.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
from rasterio.warp import transform as warp_transform

from fcp_pipeline.global_repeated_atlas import build_repeated_atlas_schedule
from fcp_pipeline.global_rgfca_engine import (
    _rank_edges, _rank_edges_matrix,
    _raw_jsd_from_source_rows, _raw_jsd_matrix_from_source_rows,
    build_pairwise_jsd_cache, canonical_colour_pool, null_source_row_matrix,
)

ROOT=Path(__file__).resolve().parents[2]
PANEL_CONTRACT=ROOT/'docs/supporting/global_rgfca_expanded_environmental_process_panel_contract_v1.json'
PANEL_EXECUTION=ROOT/'docs/supporting/global_rgfca_expanded_environmental_panel_execution_v1.json'
SOIL_SOURCE=ROOT/'docs/supporting/global_rgfca_soilgrids_expanded_source_verification_v1.json'
WORLDCLIM_SCRIPT=ROOT/'scripts/analysis/run_global_rgfca_worldclim_edge_secondary.py'
PROPERTIES=('phh2o','clay','soc','cec')
DEPTHS=(('0-5cm',5.0),('5-15cm',10.0),('15-30cm',15.0))


def sha256_file(path:Path)->str:
    h=hashlib.sha256()
    with path.open('rb') as fh:
        for chunk in iter(lambda:fh.read(1<<20),b''): h.update(chunk)
    return h.hexdigest()


def load_worldclim_module():
    name='rgfca_worldclim_secondary_for_soil'
    spec=importlib.util.spec_from_file_location(name,WORLDCLIM_SCRIPT)
    if spec is None or spec.loader is None: raise RuntimeError('cannot import frozen RGFCA edge implementation')
    module=importlib.util.module_from_spec(spec); sys.modules[name]=module; spec.loader.exec_module(module)
    return module


def physical(raw:np.ndarray, ds)->np.ndarray:
    scale=float(ds.scales[0]) if ds.scales else 1.0
    offset=float(ds.offsets[0]) if ds.offsets else 0.0
    return raw.astype(np.float64,copy=False)*scale+offset


def valid_raw(raw:np.ndarray, ds)->np.ndarray:
    keep=np.isfinite(raw)
    if ds.nodata is not None: keep &= raw != float(ds.nodata)
    return keep


def derived_property_z_at_photos(paths:list[Path], latitude:np.ndarray, longitude:np.ndarray)->tuple[np.ndarray,dict[str,object]]:
    datasets=[rasterio.open(p) for p in paths]
    try:
        ref=datasets[0]
        for ds in datasets[1:]:
            if (ds.width,ds.height)!=(ref.width,ref.height): raise RuntimeError('SoilGrids depth raster shape mismatch')
            if str(ds.crs)!=str(ref.crs): raise RuntimeError('SoilGrids depth raster CRS mismatch')
            if not np.allclose(np.asarray(tuple(ds.transform)),np.asarray(tuple(ref.transform)),atol=0,rtol=0):
                raise RuntimeError('SoilGrids depth raster transform mismatch')
        weights=np.asarray([w for _,w in DEPTHS],dtype=float); weights/=weights.sum()
        count=0; total=0.0; total_sq=0.0; dmin=float('inf'); dmax=float('-inf')
        for _,window in ref.block_windows(1):
            raws=[ds.read(1,window=window) for ds in datasets]
            masks=[valid_raw(raw,ds) for raw,ds in zip(raws,datasets)]
            keep=masks[0]&masks[1]&masks[2]
            if not np.any(keep): continue
            vals=np.vstack([physical(raw[keep],ds) for raw,ds in zip(raws,datasets)])
            derived=np.sum(vals*weights[:,None],axis=0)
            count += int(len(derived)); total += float(np.sum(derived,dtype=np.float64)); total_sq += float(np.sum(derived*derived,dtype=np.float64))
            dmin=min(dmin,float(np.min(derived))); dmax=max(dmax,float(np.max(derived)))
        if count<2: raise RuntimeError('insufficient finite SoilGrids derived support')
        mean=total/count; variance=max(0.0,total_sq/count-mean*mean); sd=float(np.sqrt(variance))
        if not np.isfinite(sd) or sd<=0: raise RuntimeError('cannot standardize SoilGrids derived property')

        xs=np.asarray(longitude,dtype=float); ys=np.asarray(latitude,dtype=float)
        if ref.crs is None: raise RuntimeError('SoilGrids raster CRS missing')
        if str(ref.crs)!='EPSG:4326':
            x2,y2=warp_transform('EPSG:4326',ref.crs,xs.tolist(),ys.tolist()); coords=list(zip(x2,y2))
        else: coords=list(zip(xs.tolist(),ys.tolist()))
        samples=[]; masks=[]
        for ds in datasets:
            raw=np.asarray([v[0] for v in ds.sample(coords)],dtype=np.float64)
            masks.append(valid_raw(raw,ds)); samples.append(physical(raw,ds))
        complete=masks[0]&masks[1]&masks[2]
        out=np.full(len(latitude),np.nan,dtype=float)
        if np.any(complete):
            vals=np.vstack([s[complete] for s in samples])
            derived=np.sum(vals*weights[:,None],axis=0)
            out[complete]=(derived-mean)/sd
        audit={
            'global_finite_derived_cells':int(count),'global_mean_derived':float(mean),'global_sd_ddof0_derived':float(sd),
            'global_min_derived':dmin,'global_max_derived':dmax,
            'photo_complete_count':int(np.count_nonzero(np.isfinite(out))),'photo_complete_fraction':float(np.mean(np.isfinite(out))),
            'width':int(ref.width),'height':int(ref.height),'crs':str(ref.crs),
            'depth_weights':{'0-5cm':5/30,'5-15cm':10/30,'15-30cm':15/30},
        }
        return out,audit
    finally:
        for ds in datasets: ds.close()


def main()->int:
    parser=argparse.ArgumentParser(); parser.add_argument('--source-dir',type=Path,required=True); parser.add_argument('--output-result',type=Path,required=True); parser.add_argument('--output-species',type=Path,required=True); parser.add_argument('--output-null',type=Path,required=True); parser.add_argument('--null-batch-size',type=int,default=64); args=parser.parse_args()
    panel=json.loads(PANEL_CONTRACT.read_text()); order=json.loads(PANEL_EXECUTION.read_text()); source=json.loads(SOIL_SOURCE.read_text())
    if panel.get('status')!='postoutcome_secondary_panel_frozen_before_any_expanded_environmental_panel_result_or_new_source_payload_is_opened': raise RuntimeError('panel contract drift')
    if order.get('status')!='technical_execution_order_frozen_before_any_expanded_panel_block_colour_alignment_result': raise RuntimeError('panel execution drift')
    if source.get('status')!='pass_exact_expanded_soilgrids_source_acquisition_before_any_edaphic_colour_alignment': raise RuntimeError('SoilGrids source verification missing')
    if source.get('official_checksums_verified') is not True: raise RuntimeError('SoilGrids official checksums not verified')

    by_key={(x['property'],x['depth']):x for x in source['files']}
    property_paths={}
    for prop in PROPERTIES:
        paths=[]
        for depth,_ in DEPTHS:
            entry=by_key[(prop,depth)]; p=args.source_dir/prop/entry['filename']
            if not p.exists() or sha256_file(p)!=entry['sha256']: raise RuntimeError(f'SoilGrids source SHA drift: {prop} {depth}')
            paths.append(p)
        property_paths[prop]=paths

    wc=load_worldclim_module(); frame,inference,_secondary,_primary=wc.load_frozen_pool(); pool=canonical_colour_pool(frame)
    outer=inference['outer_schedule']; g1=inference['g1_primary']; schedule=build_repeated_atlas_schedule(pool.photo_ids,pool.species,n_outer=int(outer['observed_resamples']),species_per_outer=int(outer['species_per_resample']),photos_per_species=int(outer['photos_per_species']),minimum_pool_photos_per_species=int(inference['input_gate']['minimum_classifiable_photos_per_species']),species_seed=int(outer['species_seed']),photo_master_seed=int(outer['photo_master_seed']))
    blocks=wc.build_edge_blocks(pool,schedule,k=int(g1['k'])); total_occurrences=int(sum(len(b.edge_nodes) for b in blocks)); print(json.dumps({'stage':'edges_built','block':'edaphic_regime','edge_occurrences':total_occurrences}),flush=True)

    columns=[]; source_audit={}
    for prop in PROPERTIES:
        z,audit=derived_property_z_at_photos(property_paths[prop],pool.latitude,pool.longitude); columns.append(z); source_audit[prop]=audit; print(json.dumps({'stage':'soil_axis_ready','property':prop,'photo_complete':audit['photo_complete_count']}),flush=True)
    photo_z=np.column_stack(columns); joint=np.isfinite(photo_z).all(axis=1); print(json.dumps({'stage':'edaphic_join_ready','joint_photo_complete':int(np.count_nonzero(joint)),'photo_total':int(len(photo_z)),'joint_photo_complete_fraction':float(np.mean(joint))}),flush=True)

    species_parts=[]; distance_parts=[]; external_parts=[]; cursor=0
    for b in blocks:
        left=photo_z[b.edge_nodes[:,0]]; right=photo_z[b.edge_nodes[:,1]]; keep=np.isfinite(left).all(axis=1)&np.isfinite(right).all(axis=1); external=np.full(len(b.edge_nodes),np.nan,dtype=float)
        if np.any(keep):
            delta=left[keep]-right[keep]; external[keep]=np.sqrt(np.mean(delta*delta,axis=1))
        b.keep_external=keep; n_keep=int(np.count_nonzero(keep)); species_parts.append(b.global_species_index[keep]); distance_parts.append(b.edge_distance_km[keep]); external_parts.append(external[keep]); b.sorted_target_positions=np.arange(cursor,cursor+n_keep,dtype=np.int64); cursor+=n_keep
    if cursor==0: raise RuntimeError('no edaphic-evaluable edges')
    species_unsorted=np.concatenate(species_parts); distance_unsorted=np.concatenate(distance_parts); external_unsorted=np.concatenate(external_parts); sort_order=np.argsort(species_unsorted,kind='stable'); inverse=np.empty(len(sort_order),dtype=np.int64); inverse[sort_order]=np.arange(len(sort_order),dtype=np.int64); species_sorted=species_unsorted[sort_order]; distance_sorted=distance_unsorted[sort_order]; external_sorted=external_unsorted[sort_order]; species_counts=np.bincount(species_sorted,minlength=len(pool.species_labels)).astype(np.int64)
    cursor=0
    for b in blocks:
        keep=b.keep_external; assert keep is not None; n_keep=int(np.count_nonzero(keep)); b.sorted_target_positions=inverse[np.arange(cursor,cursor+n_keep,dtype=np.int64)]; cursor+=n_keep

    cache=build_pairwise_jsd_cache(pool); identity=np.arange(len(pool.photo_ids),dtype=np.int64); observed_scores=np.empty(len(species_sorted),dtype=float)
    for b in blocks:
        keep=b.keep_external; target=b.sorted_target_positions; assert keep is not None and target is not None; raw=_raw_jsd_from_source_rows(cache,b.edge_nodes,identity); rank=_rank_edges(raw,b.edge_species_slices); observed_scores[target]=rank[keep]
    minimum_edges=int(panel['primary_edge_test']['minimum_edges_per_species']); minimum_species=int(panel['primary_edge_test']['minimum_evaluable_species']); observed_values,observed_species_ids=wc._species_partial_values(observed_scores,external_sorted,distance_sorted,species_counts,minimum_edges=minimum_edges)
    if len(observed_values)<minimum_species: raise RuntimeError(f'only {len(observed_values)} edaphic-evaluable species')
    observed_mean=float(np.mean(observed_values)); observed_median=float(np.median(observed_values)); observed_positive=float(np.mean(observed_values>0)); print(json.dumps({'stage':'observed_block','block':'edaphic_regime','n_species':int(len(observed_values)),'mean_partial_rho':observed_mean,'median_partial_rho':observed_median,'positive_fraction':observed_positive}),flush=True)

    batch_size=int(args.null_batch_size); null_indices=np.arange(999,dtype=np.int64); null_stat=np.empty(999,dtype=float)
    for batch_start in range(0,999,batch_size):
        batch_indices=null_indices[batch_start:batch_start+batch_size]; source_rows=null_source_row_matrix(pool,batch_indices,master_seed=int(g1['null_master_seed'])); batch_scores=np.empty((len(batch_indices),len(species_sorted)),dtype=float)
        for b in blocks:
            keep=b.keep_external; target=b.sorted_target_positions; assert keep is not None and target is not None; raw=_raw_jsd_matrix_from_source_rows(cache,b.edge_nodes,source_rows); rank=_rank_edges_matrix(raw,b.edge_species_slices); batch_scores[:,target]=rank[:,keep]
        for local_index,null_index in enumerate(batch_indices):
            values,ids=wc._species_partial_values(batch_scores[local_index],external_sorted,distance_sorted,species_counts,minimum_edges=minimum_edges)
            if not np.array_equal(ids,observed_species_ids): raise RuntimeError(f'null {int(null_index)} changed edaphic species identity')
            null_stat[int(null_index)]=float(np.mean(values))
        print(json.dumps({'stage':'null_batch_complete','block':'edaphic_regime','first_null':int(batch_indices[0]),'last_null':int(batch_indices[-1]),'completed':int(batch_indices[-1])+1,'total':999}),flush=True)
    p_upper=float((1+np.count_nonzero(null_stat>=observed_mean))/1000); alpha=float(panel['primary_edge_test']['alpha']); raw_positive=bool(observed_mean>0 and p_upper<alpha)
    args.output_result.parent.mkdir(parents=True,exist_ok=True); args.output_species.parent.mkdir(parents=True,exist_ok=True); args.output_null.parent.mkdir(parents=True,exist_ok=True)
    pd.DataFrame({'species':[pool.species_labels[int(i)] for i in observed_species_ids],'species_index':observed_species_ids,'retained_edge_occurrences':species_counts[observed_species_ids],'partial_rho_colour_vs_edaphic_regime_given_distance':observed_values}).to_csv(args.output_species,index=False,lineterminator='\n'); pd.DataFrame({'null_index':null_indices,'global_mean_species_partial_rho':null_stat}).to_csv(args.output_null,index=False,lineterminator='\n')
    payload={'protocol':panel['protocol'],'status':'complete_edaphic_regime_block_raw_pending_full_panel_holm','block':'edaphic_regime','inferential_role':panel['inferential_role'],'final_panel_support_decision_available':False,'reason_final_support_pending':'Holm adjustment requires all evaluable fixed blocks in the nonselective five-block panel','predictor_source':'exact official-checksum-verified SoilGrids 2.0 5-km aggregated mean, thickness-weighted 0-30 cm','axes':['phh2o_0_30cm','clay_0_30cm','soc_0_30cm','cec_0_30cm'],'observed_mean_species_partial_rho':observed_mean,'observed_median_species_partial_rho':observed_median,'observed_positive_species_fraction':observed_positive,'n_evaluable_species':int(len(observed_values)),'minimum_evaluable_species':minimum_species,'total_scheduled_edge_occurrences':total_occurrences,'block_evaluable_edge_occurrences':int(len(species_sorted)),'block_edge_occurrence_coverage_fraction':float(len(species_sorted)/total_occurrences),'joint_photo_complete_count':int(np.count_nonzero(joint)),'joint_photo_complete_fraction':float(np.mean(joint)),'source_audit':source_audit,'null_permutations':999,'null_indices_exact_0_998':True,'null_mean':float(np.mean(null_stat)),'null_q025':float(np.quantile(null_stat,.025)),'null_q975':float(np.quantile(null_stat,.975)),'raw_p_upper':p_upper,'raw_nominal_positive_before_panel_Holm':raw_positive,'alpha':alpha,'p_holm':None,'panel_supported':None,'individual_variable_decomposition_open':False,'causal_language_allowed':False,'execution_audit':{'all_200_outer_realizations_used':len(blocks)==200,'same_RGFCA_edge_geometry':True,'repeated_edge_occurrences_deduplicated':False,'missing_environment_recoded_to_zero':False,'four_fixed_derived_soil_axes_only':True,'same_999_RGFCA_colour_nulls':True,'block_result_cannot_stop_later_fixed_blocks':True},'lineage':{'panel_contract_sha256':sha256_file(PANEL_CONTRACT),'panel_execution_sha256':sha256_file(PANEL_EXECUTION),'soil_source_verification_sha256':sha256_file(SOIL_SOURCE),'measurement_manifest_sha256':sha256_file(wc.MEASUREMENT),'measured_table_sha256':sha256_file(wc.MEASURED),'primary_g1_result_sha256':sha256_file(wc.PRIMARY_G1),'species_result_sha256':sha256_file(args.output_species),'null_result_sha256':sha256_file(args.output_null)},'files':{'species':str(args.output_species),'null':str(args.output_null)}}; args.output_result.write_text(json.dumps(payload,indent=2)+'\n'); print(json.dumps(payload,indent=2),flush=True); return 0

if __name__=='__main__': raise SystemExit(main())
