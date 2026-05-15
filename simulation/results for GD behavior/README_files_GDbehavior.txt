This folder stores outputs from:
  simulation/linesearch_vs_fixed_4etas.R

Working directory requirement
- Set working directory to BTAP before running the producing script.

Simulation scope
- Methods: fixed-step GD and line-search GD.
- Initial step sizes: eta_init = eta_mult * eta0, with eta_mult in {10, 5, 1, 0.2}.
- Families: bernoulli, poisson, gaussian.
- Network sizes: n in {500, 1000, 2000, 4000}.
- For each replication, the script runs 8 cases:
  2 methods x 4 eta_mult values.

File naming pattern (one pair per family x n)
1) fs_vs_ls_4etas_{family}_n{n}_run.csv
2) fs_vs_ls_4etas_{family}_n{n}_summary.csv

================================================================================
run.csv: one row per (rep, method, eta_mult) case
================================================================================
Columns and meanings
- rep:
  Replication index.
- seed:
  Random seed used for this replication.
- n:
  Number of nodes in the simulated network.
- k:
  Latent dimension.
- family:
  Edge distribution family: bernoulli / poisson / gaussian.
- method:
  GD update type: fixed or linesearch.
- eta0:
  Baseline step size computed from RA-SVT initialization.
- eta_mult:
  Step-size multiplier in {10, 5, 1, 0.2}.
- eta_init:
  Initial step size used in this case, eta_init = eta_mult * eta0.
- eta_times_n:
  n-scaled step size, eta_init * n.
- eps_grad:
  GD stopping threshold on max(abs(grad_Y)).
- max_gd_steps:
  Maximum allowed GD iterations.
- beta_ls:
  Line-search backtracking shrink factor.
- c_ls:
  Line-search condition constant.
- Rprime:
  Backtracking cap per GD iteration.
  - For fixed-step GD: 0 (no line search).
  - For line-search GD: positive integer cap.
- error:
  Run-level error flag.
  - 1: optimizer returned an error for this case.
  - 0: normal return.
- error_message:
  Error message text when error = 1; empty string otherwise.
- elapsed_sec:
  runtime (seconds) for this case.
- num_of_steps:
  Number of GD iterations executed.
- converged:
  Convergence indicator (1/0), based on eps_grad stopping criterion.
- exploded:
  Explosion indicator (1/0) from optimizer diagnostics.
- score_Y_max_abs:
  Terminal max absolute entry of grad_Y.
- bt_mean:
  Mean number of backtracks per GD iteration.
  NA for fixed-step GD.
- bt_max:
  Maximum number of backtracks observed in one GD iteration.
  NA for fixed-step GD.
- bt_frac0:
  Fraction of GD iterations with zero backtracking.
  NA for fixed-step GD.
- bt_frac_hit_cap:
  Fraction of GD iterations that hit exactly Rprime backtracks.
  NA for fixed-step GD.
- ls_fail_frac:
  Fraction of GD iterations where line-search conditions were still unsatisfied
  after hitting the Rprime cap.
  NA for fixed-step GD.
- record_id:
  Sequential row id in this file (1-based). Used only as a record index.

Notes on NA values
- If error = 1, many algorithm-output columns can be NA (for example
  num_of_steps, score_Y_max_abs, backtracking metrics).
- For fixed-step GD, line-search-specific metrics (bt_*, ls_fail_frac) are NA by design.

================================================================================
summary.csv: aggregated statistics by (method, eta_mult)
================================================================================
Each row summarizes all run.csv rows in one (method, eta_mult) cell, i.e., averaged over all replications.

Columns and meanings
- method:
  fixed or linesearch.
- eta_mult:
  Step-size multiplier in {10, 5, 1, 0.2}.
- reps:
  Number of rows aggregated in this cell.
- error_rate:
  Mean of error (fraction with error = 1).
- conv_rate:
  Mean of converged, convergence fraction.
- explode_rate:
  Mean of exploded, explosion fraction.
- mean_elapsed:
  Mean elapsed_sec (seconds).
- median_elapsed:
  Median elapsed_sec (seconds).
- mean_steps:
  Mean num_of_steps.
- mean_score_Y:
  Mean score_Y_max_abs.
- mean_bt_frac0:
  Mean bt_frac0.
- mean_bt_max:
  Mean bt_max.
- mean_bt_hit_cap:
  Mean bt_frac_hit_cap.
- mean_ls_fail_frac:
  Mean ls_fail_frac.

Notes on interpretation
- Backtracking summaries are meaningful for method = linesearch. For method = fixed, bt/ls-related summaries can be NA because source columns are NA.
