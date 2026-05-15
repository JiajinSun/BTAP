# Simulation Studies and Data Analysis of "Bridging Theory and Practice: Statistical Inference of Latent Space Models for Networks"

## Project Description
This folder contains code and outputs for the simulation studies and data analysis in the paper
"Bridging Theory and Practice: Statistical Inference of Latent Space Models for Networks".

## Working Directory Requirement
Set the working directory to `BTAP` before running any script in this project.
Example: `getwd()` should return the full path of the `BTAP` folder.

## Repository Structure

### Core Algorithm Functions (`functions_RASVT_PGDwBLS.R`)
- `functions_RASVT_PGDwBLS.R`: Main function library for this project (stored in the BTAP root). It includes:
  - `RASVT` (Range-Adaptive Singular Value Thresholding, Algorithm 2 in the paper) for initialization.
  - `RASVT_PGDwBLS` (Projected Gradient Descent with Backtracking Line Search, Algorithm 1 in the paper). 
    - "method" option can be specified to {'fixed','linesearch'} for fixed step size and line search. 
  - Utility functions used in the above two algorithm functions, including
    - `get_link_functions`: return link and inverse-link functions for each distribution family;
    - `negloglik_term`, `nll_term_mpfr_vec`: compute edgewise negative log-likelihood in ordinary precision and in high precision (MPFR);
    - `check_linesearch_conditions`: evaluate line-search conditions in ordinary precision;
    - `recompute_cond_global_mpfr_stablediff`, `recompute_cond_individual_mpfr_stablediff`: recompute global line-search and node-wise line-search condition in high precision.
  - Utility functions used by the simulation and plotting scripts, including     
    - `eta_cornodiv6_from_init`: calculate an initial step size eta_0 from initial estimate;
    - `DLii`: compute the i-th diagonal expected negative Hessian block;
    - `build_asymptotic_record`: extract (Z,alpha) true, (Z,alpha) estimate, and calculate t(z_{q,11}) from the GD estimate.

### Simulation Folder (`simulation/`)
- `simulation/linesearch_vs_fixed_4etas.R`:
  Runs fixed-step vs line-search GD behavior studies (8 cases per replication: 2 methods x 4 `eta_init/eta_0` values) across families and sample sizes.
- `simulation/linesearch_eta0.R`:
  Runs asymptotic-distribution studies using line-search GD with `eta_init = eta_0`.
- `simulation/results for GD behavior/`:
  Stores run-level and summary CSV outputs for GD convergence/backtracking diagnostics.
- `simulation/results for asymptotic distribution/`:
  Stores per-replication asymptotic records (`.rds`) used by QQ/histogram plotting scripts.
- `simulation/Fig 1/`, `simulation/Fig 2 boxplot/`, `simulation/Fig 3/`, `simulation/Fig 4/`:
  Figure-specific scripts and output files for trajectory plots, boxplots, and asymptotic QQ/hist plots.
- More detailed file description for the full simulation pipeline is given in `simulation/README.md`.

### Data Analysis Folder (`data analysis/`)
- `data analysis/NYbike.rda` and `data analysis/201908-citibike-tripdata_aug01.csv`: Citi Bike trip records on Aug 1, 2019, with network data, borough labels, and station geolocation information.
- `data analysis/1_preprocess.R` and`data analysis/2_nybike_fit.R`:
  Preprocesses the network data and fits the latent-space model for the two hours.
- `data analysis/3_two_sample_tests.R`:
  Runs z inner product and alpha two-sample tests between the two hours.
- `data analysis/4_z_visualize.R`, `data analysis/5_innerprod_heatmap.R`, `data analysis/6_innerprod_reject_map.R`, `data analysis/7_alpha_visualize.R`, `data analysis/8_alpha_reject_map.R`:
  Visualization scripts that generate latent-position plots, inner-product heatmaps, and rejection geolocation maps in the paper and supplement.
- `data analysis/results/`:
  Stores preprocessing outputs, model-fit outputs, test-result tables, and helper CSVs for plotting.
- `data analysis/plots/`:
  Stores generated PDF figures for the data-analysis section.
- More detailed file description for the full data-analysis pipeline is given in `data analysis/README.md`.
  
