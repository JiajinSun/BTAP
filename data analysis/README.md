# Data Analysis of "Bridging Theory and Practice: Statistical Inference of Latent Space Models for Networks"

## Project Description
This folder contains code and output files for data analysis in the paper.

## Working Directory Requirement
Set working directory to `BTAP` before running scripts in this folder; `getwd()` should return the `BTAP` folder path.

## Data Files

### `NYbike.rda`
- Network object used for model fitting, contains hourly bike riding records on August 1, 2019. This rda is derived from raw data by https://github.com/yatian20/DLSM.
- Contains:
  - `NYbike1`: hourly network array (node x node x hour, 782 x 782 x 24 dimension),
  - `borough`: borough label for each original node.

### `201908-citibike-tripdata_aug01.csv`

- Raw Citi Bike trip records for August 1, 2019.
- Used to rebuild station geolocation information (latitude, longitude).

## R Scripts

### `1_preprocess.R`
- Purpose: preprocess the network data in `NYbike.rda` to the two-hour-data of interest (`hour1 = 8`, `hour2 = 18`).
- Output:
  - `data analysis/results/nybike_data.rds`

### `2_nybike_fit.R`
- Purpose: fit the latent-space model for hour 08 and hour 18 using RA-SVT initialization + PGDwBLS.
- Output:
  - `data analysis/results/nybike_fit.rds`

### `3_two_sample_tests.R`
- Purpose: perform two-sample inference between hour 08 and hour 18.
- Outputs:
  - results of two-sample test for alpha: `data analysis/results/nybike_08_vs_18_alpha.csv`
  - results of two-sample test for latent inner-product (each row is a pair of nodes): `data analysis/results/nybike_08_vs_18_innerprod.csv`
  - node-level  summarized results of two-sample test for latent inner-product (count of total tests, reject tests, rejection rate for each node): `data analysis/results/nybike_innerprod_rejections_by_node.csv`
  
  - these CSVs are used by scripts 5, 6, and 8.

### `4_z_visualize.R`
- Purpose: make latent position plots and true geolocation plots.
- Outputs:
  - `data analysis/results/nybike_station_locations_latlon.csv`
  - `data analysis/plots/nybike_true_locations.pdf`
  - `data analysis/plots/nybike_hour_08_latent.pdf`
  - `data analysis/plots/nybike_hour_18_latent.pdf`
- Paper figures produced:
  - Figure 5(a)
  - Figure S9

### `5_innerprod_heatmap.R`
- Purpose: produce heatmaps for latent inner-product structures in the two hour windows, and differences of latent inner-product between the two hours (all pairs and rejected only); also produce heat map node order to be used for later codes (for alpha two-sample test heat maps).
- Outputs:
  - `data analysis/plots/nybike_innerprod_hour1.pdf`
  - `data analysis/plots/nybike_innerprod_hour2.pdf`
  - `data analysis/plots/nybike_innerprod_difference.pdf`
  - `data analysis/plots/nybike_innerprod_rejected_only.pdf`
  - `data analysis/results/nybike_innerprod_heatmap_order.csv`
- Paper figures produced:
  - Figure 5(b), Figure 5(c)
  - Figure 6(a), Figure 6(b)

### `6_innerprod_reject_map.R`
- Purpose: show station-level fraction of rejected pairwise inner-product tests in geolocation map.
- Output:
  - `data analysis/plots/nybike_station_locations_innerprod_fraction.pdf`
- Paper figure produced:
  - Figure 6(c)

### `7_alpha_visualize.R`
- Purpose: show estimated baseline intensity (`exp(alpha)`) for each hour in geolocation map.
- Outputs:
  - `data analysis/plots/nybike_station_locations_alpha_08.pdf`
  - `data analysis/plots/nybike_station_locations_alpha_18.pdf`
- Paper figures produced:
  - Figure S10(a), Figure S10(b)

### `8_alpha_reject_map.R`
- Purpose: show BH-significant node-wise alpha differences with p-value-based highlighting in geolocation map.
- Output:
  - `data analysis/plots/nybike_station_locations_alpha_pvalue.pdf`
- Paper figures produced:
  - Figure S10(c)
