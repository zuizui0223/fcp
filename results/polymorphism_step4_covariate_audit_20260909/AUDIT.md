# Step 4 covariate availability audit

Header-only audit; no D/covariate association was computed.

CSV headers scanned: **118**

## taxonomy
- `data/frozen/frozen_34species_coexistence_segregation_v22.csv`: family
- `data/frozen/frozen_34species_five_metric_dataset.csv`: family
- `data/frozen/global_monte_carlo_candidate_photos_v1.csv`: inat_taxon_id
- `data/frozen/global_monte_carlo_candidate_species_audit_v1.csv`: inat_taxon_id, wrong_taxon
- `data/frozen/global_monte_carlo_capacity_scan_selected_species_v2.csv`: inat_taxon_id
- `data/frozen/global_monte_carlo_capacity_scan_selected_species_v3.csv`: inat_taxon_id
- `data/frozen/global_monte_carlo_capacity_scan_species_audit_v2.csv`: inat_taxon_id, wrong_taxon
- `data/frozen/global_monte_carlo_capacity_scan_species_audit_v3.csv`: inat_taxon_id, wrong_taxon
- `data/frozen/global_monte_carlo_species_discovery_combined_v2_species_v1.csv`: inat_taxon_id
- `data/frozen/global_monte_carlo_species_discovery_species_v1.csv`: inat_taxon_id
- `data/frozen/global_monte_carlo_species_discovery_v2_species_v1.csv`: inat_taxon_id
- `data/frozen/jbi_ch1_photo_source_manifest.csv`: inat_taxon_id, inat_taxon_name
- `data/frozen/jbi_ch1_photo_split_v1.csv`: inat_taxon_id, inat_taxon_name
- `data/frozen/jbi_cs_climate_analysis_v22.csv`: family
- `data/frozen/jbi_cs_evidence_freeze_v22.csv`: family
- `data/frozen/random_photo_first_candidate_pool_v1.csv`: inat_taxon_id
- `data/frozen/random_photo_first_h7_fresh_metadata_species_support_v1.csv`: inat_taxon_id
- `data/frozen/random_photo_first_h7_fresh_metadata_target_audit_v1.csv`: inat_taxon_id, target_cell_order, wrong_taxon
- `data/frozen/random_photo_first_h7_fresh_metadata_v1.csv`: inat_taxon_id, target_cell_order
- `data/frozen/random_photo_first_h7_species_cells_v1.csv`: inat_taxon_id, target_cell_order
- `data/frozen/random_photo_first_h8_design_species_cells_v1.csv`: inat_taxon_id, h8_target_order
- `data/frozen/random_photo_first_h8_design_species_v1.csv`: inat_taxon_id
- `data/frozen/random_photo_first_h9_fresh_metadata_species_audit_v1.csv`: inat_taxon_id, wrong_taxon
- `data/frozen/random_photo_first_h9_fresh_metadata_v1.csv`: inat_taxon_id, h9_selection_order
- `data/frozen/random_photo_first_h9_metadata_species_v1.csv`: inat_taxon_id
- `data/frozen/rgfca_independent_global_axis_prequalified_photo_manifest_20260909.csv`: inat_taxon_id, selected_group_order
- `data/frozen/rgfca_independent_global_axis_species_audit_20260909.csv`: inat_taxon_id
- `data/frozen/rgfca_reserve_geometry_audit_v1.csv`: inat_taxon_id
- `data/frozen/rgfca_sharedness_v2_geometry_photos_v1.csv`: species_order, inat_taxon_id, photo_order
- `data/frozen/rgfca_sharedness_v2_geometry_species_audit_v1.csv`: species_order, inat_taxon_id
- `data/frozen/rgfca_sharedness_v2_high_depth_capacity_pilot_v1.csv`: inat_taxon_id
- `data/derived/global_monte_carlo_measured_photos_v1.csv`: inat_taxon_id
- `data/derived/global_monte_carlo_measurement_species_support_v1.csv`: inat_taxon_id
- `data/derived/global_rgfca_environmental_holm_power_s3_v1.csv`: family, family_any_discovery_probability
- `data/derived/global_rgfca_environmental_interaction_power_supplement_v1.csv`: family_any_discovery_probability
- `data/derived/random_photo_first_h7_metadata_feasibility_species_v1.csv`: inat_taxon_id
- `data/derived/random_photo_first_measured_photos_v1.csv`: inat_taxon_id
- `data/derived/rgfca_reserve_replication_measured_photos_v1.csv`: inat_taxon_id
- `data/derived/rgfca_reserve_replication_species_support_v1.csv`: inat_taxon_id
- `data/derived/rgfca_reserve_replication_species_v1.csv`: inat_taxon_id

## life_form
- no tracked CSV header match

## pollination
- no tracked CSV header match

## range
- `data/frozen/frozen_34species_five_metric_dataset.csv`: pca_hull_area
- `data/frozen/global_monte_carlo_candidate_species_audit_v1.csv`: maximum_span_km_after_observer_cap
- `data/frozen/global_monte_carlo_capacity_scan_selected_species_v2.csv`: maximum_span_km
- `data/frozen/global_monte_carlo_capacity_scan_selected_species_v3.csv`: maximum_span_km
- `data/frozen/global_monte_carlo_capacity_scan_species_audit_v2.csv`: maximum_span_km
- `data/frozen/global_monte_carlo_capacity_scan_species_audit_v3.csv`: maximum_span_km
- `data/frozen/jbi_cs_climate_analysis_v22.csv`: pca_hull_area
- `data/frozen/random_photo_first_h9_metadata_species_v1.csv`: h7_maximum_span_km
- `data/frozen/rgfca_reserve_geometry_audit_v1.csv`: maximum_span_km
- `data/derived/global_monte_carlo_measured_photos_v1.csv`: palette_count_orange, flower_fraction_orange, colour_yellow_orange
- `data/derived/random_photo_first_h4_within_species_capacity_v1.csv`: yellow_orange_n
- `data/derived/random_photo_first_measured_photos_v1.csv`: palette_count_orange, flower_fraction_orange
- `data/derived/rgfca_reserve_replication_measured_photos_v1.csv`: palette_count_orange, background_palette_count_orange, flower_fraction_orange, colour_yellow_orange

## latitude
- `data/frozen/global_monte_carlo_candidate_photos_v1.csv`: latitude
- `data/frozen/global_monte_carlo_platform_effort_cells_v1.csv`: swlat, nelat
- `data/frozen/global_monte_carlo_sampling_availability_cells_v1.csv`: swlat, nelat
- `data/frozen/global_monte_carlo_species_discovery_round_audit_v1.csv`: cumulative_species
- `data/frozen/global_monte_carlo_species_discovery_v2_round_audit_v1.csv`: cumulative_v1_plus_v2_species
- `data/frozen/jbi_ch1_photo_source_manifest.csv`: latitude
- `data/frozen/jbi_ch1_photo_split_v1.csv`: latitude
- `data/frozen/random_photo_first_candidate_pool_cell_audit_v1.csv`: swlat, nelat
- `data/frozen/random_photo_first_candidate_pool_v1.csv`: latitude
- `data/frozen/random_photo_first_h7_fresh_metadata_v1.csv`: latitude, target_cell_center_latitude
- `data/frozen/random_photo_first_h7_species_cells_v1.csv`: cell_center_latitude
- `data/frozen/random_photo_first_h9_fresh_metadata_v1.csv`: latitude
- `data/frozen/rgfca_independent_global_axis_prequalified_photo_manifest_20260909.csv`: latitude
- `data/frozen/rgfca_independent_global_axis_species_audit_20260909.csv`: max_qualifying_group_centroid_separation_km
- `data/frozen/rgfca_sharedness_v2_geometry_photos_v1.csv`: latitude
- `data/derived/global_monte_carlo_measured_photos_v1.csv`: latitude
- `data/derived/global_rgfca_design_power_s1_v1.csv`: simulation_replicates
- `data/derived/global_rgfca_environmental_holm_power_s3_v1.csv`: correlation
- `data/derived/global_rgfca_worldclim_scale_moderation_replication_species_v1.csv`: mean_pairwise_H1_cell_centroid_great_circle_km, log1p_mean_pairwise_H1_cell_centroid_great_circle_km
- `data/derived/random_photo_first_h1_sensitivities_v1.csv`: n_sinlat
- `data/derived/random_photo_first_h5_species_primary_v1.csv`: mean_pairwise_h1_cell_centroid_great_circle_km, max_pairwise_h1_cell_centroid_great_circle_km
- `data/derived/random_photo_first_h6_species_primary_v1.csv`: mean_pairwise_h1_cell_centroid_great_circle_km
- `data/derived/random_photo_first_measured_photos_v1.csv`: latitude
- `data/derived/rgfca_global_signal_spatial_recurrence_20260909/exploratory_flower_side_candidate_cells.csv`: lat
- `data/derived/rgfca_reserve_replication_measured_photos_v1.csv`: latitude

## Freeze implications

- genus can be derived from the frozen binomial species name.
- range-size and latitude-centroid proxies can be built from colour-blind frozen metadata; their exact estimators must be frozen before opening D associations.
- family, life form, and pollination mode require a tracked source if none is listed above; they must not be back-filled after Step 4 outcomes are opened.
