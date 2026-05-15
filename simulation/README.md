# Simulation Studies of "Bridging Theory and Practice: Statistical Inference of Latent Space Models for Networks"

## Project Description
This folder contains code and output files for simulation studies in the paper.

## Working Directory Requirement
Set working directory to `BTAP` before running scripts in this folder; `getwd()` should return the `BTAP` folder path.

## Repository Structure

### Simulation Studies for Fixed-Step vs Line-Search GD Behavior
- `linesearch_vs_fixed_4etas.R`: Runs simulations corresponding to Section 6.1 and Supplementary Material Section J.1 of the paper. For each setup with distribution in `{poisson, gaussian, bernoulli}` and `n` in `{500, 1000, 2000, 4000}`, each replication runs 8 GD cases:
  - method in `{fixed, linesearch}`
  - $\eta_{init} / \eta_0$ in `{10, 5, 1, 1/5}`
  The script summarizes GD behavior, especially convergence and backtracking.
- `results for GD behavior/`: Output files produced by `linesearch_vs_fixed_4etas.R`.
  - `fs_vs_ls_4etas_{family}_n{n}_run.csv`: One row per `(rep, method, eta_mult)` case with GD convergence and backtracking diagnostics.
  - `fs_vs_ls_4etas_{family}_n{n}_summary.csv`: Aggregated summaries by `(method, eta_mult)`.
  - `README_files_GDbehavior.txt`: File-format notes for csv files in this folder.

### Simulation Studies for Asymptotic Distribution
- `linesearch_eta0.R`: Runs simulations corresponding to Section 6.2 and Supplementary Material Section J.2 of the paper. For each setup with distribution in `{poisson, gaussian, bernoulli}` and `n` in `{500, 1000, 2000, 4000}`, each replication runs line-search GD with $\eta_{init} = \eta_0$.
- `results for asymptotic distribution/`: Output files produced by `linesearch_eta0.R`.
  - `eta0_only_{family}_n{n}_asymptotic_records.rds`: Per-replication records, which includes setup information as well as estimator information such as $(Z^\star,\alpha^\star, \hat Z, \hat\alpha)$ used for asymptotic plotting scripts.
  - `README_files_asymptoticdistribution.txt`: File-format notes for rds data in this folder.

### Figure Scripts and Outputs

#### Figure 1 (`Fig 1/`)
- `fs_vs_ls_trajectory.R`: Generates score and likelihood trajectory CSV files for one-seed run with $\eta_{init} / \eta_0$ in $\{10, 5\}$ for fixed-step and line-search GD.
- `plot_score_trajectory-P.R`, `plot_score_trajectory-B.R`, `plot_score_trajectory-G.R`: Plot score trajectories for Poisson, Bernoulli, and Gaussian settings. Produce Figures 1, S1, S2 in the paper.
- `scoreY_*.pdf`: Generated trajectory figures.

#### Figure 2 Boxplots (`Fig 2 boxplot/`)
- `boxplot_gd_steps_n1000.R`: Produces boxplots of GD iteration counts versus $\eta_{init} / \eta_0$ (line-search runs only, `n = 1000`), i.e, Figures 2, S3, S4 panel (a) in the paper.
- `boxplot_avg_backtracks_n1000.R`: Produces boxplots of average backtracking steps per GD iteration versus $\eta_{init} / \eta_0$ (line-search runs only, `n = 1000`), i.e, Figures 2, S3, S4 panel (b) in the paper.
- `boxplot_*.pdf`: Generated boxplot figures by distribution.

#### Figure 3 (`Fig 3/`)
- `plot_AN_t11.R`: Reads asymptotic records from the folder `results for asymptotic distribution` and generates QQ plots for standardized $t(\hat z_{q,11})$, i.e., Figures 3, S5, S6 in the paper.
- `qq_t_z11_*.pdf`, `hist_t_z11_*.pdf`: Generated Figures.

#### Figure 4 (`Fig 4/`)
- `plot_AN_tmu12.R`: Reads asymptotic records from the folder `results for asymptotic distribution` and generates QQ plots for standardized $t(\hat \mu_{12})$, i.e., Figures 4, S7, S8 in the paper.
- `qq_t_muTheta12_*.pdf`, `hist_t_muTheta12_*.pdf`: Generated Figures.
