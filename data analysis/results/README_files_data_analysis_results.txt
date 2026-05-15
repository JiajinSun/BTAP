# Data Analysis Results CSV Files

This file explains the CSV outputs in `BTAP/data analysis/results`.

General notes:
- `filtered_node` = node index after node filtering in `1_preprocess.R` (the node indexing used by analysis/plots).
  Specifically, nodes are kept when their symmetric-network degree is > 1 in BOTH analyzed hour slices (hour1 = 8 and hour2 = 18); the symmetric-network here is constructed from nybike_data.rds.
- `original_node` = node index before filtering in the raw station network.
- `hour1` and `hour2` correspond to the two analyzed time windows (morning 8-9 vs evening 18-19).

------------------------------------------------------------
1) nybike_08_vs_18_alpha.csv
------------------------------------------------------------
Produced by:
- `data analysis/3_two_sample_tests.R`

What each row is:
- one row for each filtered station; the row gives node-wise two-sample test for baseline parameter alpha.

Columns:
- `filtered_node`: station index in the filtered network.
- `original_node`: corresponding station index in the original network.
- `borough`: borough label of the station.
- `hour1`: first hour label (e.g., "08:00-09:00").
- `hour2`: second hour label (e.g., "18:00-19:00").
- `alpha_hat_1`: estimated baseline alpha at `hour1`.
- `alpha_hat_2`: estimated baseline alpha at `hour2`.
- `diff`: estimated difference (`alpha_hat_1 - alpha_hat_2`).
- `se`: estimated standard error of `diff`.
- `z_stat`: z statistic for testing equality of the two alphas, computed as `diff / se`
- `p_value`: raw two-sided p-value.
- `q_value_bh`: Benjamini-Hochberg adjusted p-value.

------------------------------------------------------------
2) nybike_08_vs_18_innerprod.csv
------------------------------------------------------------
Produced by:
- `data analysis/3_two_sample_tests.R`

What each row is:
- one row for each unordered station pair `(i, j)` with `i < j`; each row gives two-sample test for latent inner-product difference across hours.

Columns:
- `filtered_node_i`: first station index in filtered network.
- `original_node_i`: corresponding first-station index in original network.
- `borough_i`: borough label for first station.
- `filtered_node_j`: second station index in filtered network.
- `original_node_j`: corresponding second-station index in original network.
- `borough_j`: borough label for second station.
- `hour1`: first hour label (e.g., "08:00-09:00").
- `hour2`: second hour label (e.g., "18:00-19:00").
- `innerprod_hat_1`: estimated latent inner product for pair `(i, j)` at `hour1`.
- `innerprod_hat_2`: estimated latent inner product for pair `(i, j)` at `hour2`.
- `diff`: estimated difference (`innerprod_hat_1 - innerprod_hat_2`).
- `se`: estimated standard error of `diff`.
- `z_stat`: z statistic for testing equality of the two inner products, computed as `diff / se`
- `p_value`: raw two-sided p-value.
- `q_value_bh`: Benjamini-Hochberg adjusted p-value.

------------------------------------------------------------
3) nybike_innerprod_rejections_by_node.csv
------------------------------------------------------------
Produced by:
- `data analysis/3_two_sample_tests.R`

What each row is:
- One filtered station, summarizing how often pairwise inner-product tests involving this station are rejected.

Columns:
- `filtered_node`: station index in the filtered network.
- `original_node`: corresponding station index in the original network.
- `borough`: borough label.
- `rejected_count`: number of rejected pairwise tests involving this station.
- `total_tests`: total number of pairwise tests involving this station.
- `rejected_frac`: rejection fraction (`rejected_count / total_tests`).

------------------------------------------------------------
4) nybike_station_locations_latlon.csv
------------------------------------------------------------
Produced by:
- `data analysis/4_z_visualize.R`

What each row is:
- One filtered station with its map coordinates.

Columns:
- `filtered_node`: station index in the filtered network.
- `station_id`: Citi Bike station ID from raw trip data.
- `lat`: station latitude.
- `lon`: station longitude.

------------------------------------------------------------
5) nybike_innerprod_heatmap_order.csv
------------------------------------------------------------
Produced by:
- `data analysis/5_innerprod_heatmap.R`

What each row is:
- One filtered station, with its ordering used to draw the latent-inner-product heatmaps.

Columns:
- `heatmap_rank`: plotting order index used in the heatmap axes.
- `filtered_node`: station index in the filtered network.
- `original_node`: corresponding station index in the original network.
- `borough`: borough label.
- `z1_hour1`: first latent coordinate (hour 1 fit), used for ordering/tie-breaking.
- `z2_hour1`: second latent coordinate (hour 1 fit), used for ordering/tie-breaking.
