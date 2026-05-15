This folder stores scripts and outputs for Figure 1 style trajectory plots.

Working directory requirement
- Set working directory to BTAP before running the Figure 1 scripts.
- Example: getwd() should return the BTAP folder path.

Main scripts
- fs_vs_ls_trajectory.R
  Generates per-iteration trajectory CSV and per-case summary CSV.
- plot_score_trajectory-P.R / -B.R / -G.R
  Read trajectory CSVs and generate score trajectory PDFs.

File naming pattern
1) fs_vs_ls_{family}_n{n}_trajectory.csv
2) fs_vs_ls_{family}_n{n}_summary.csv
3) scoreY_fsGD_{family}_n{n}.pdf
4) scoreY_lsGD_{family}_n{n}.pdf

================================================================================
trajectory.csv: one row per GD iteration
================================================================================
Columns and meanings
- iter:
  GD iteration index (starts at 1).
- loglik:
  log-likelihood at this iterate.
- score_max_abs:
  max(abs(grad_Y)) at this iterate.
- finite_score:
  Indicator for finite score_max_abs (1 = finite, 0 = non-finite).
- bt:
  Number of line-search backtracking steps associated with this row.
  Convention: row r stores transition info from iterate r-1 -> r.
  Therefore row 1 is NA by construction.
  For fixed-step GD, this is structurally 0 from row 2 onward (not an active diagnostic).
- eta_used:
  Step size used for the transition into this row (same row convention as bt).
  Row 1 is NA by construction.
  For fixed-step GD, this is always eta_init from row 2 onward.
- ls_ok:
  Line-search acceptance indicator after backtracking (1 = accepted, 0 = failed),
  using the same row convention as bt/eta_used. Row 1 is NA.
  For fixed-step GD, this is structurally 1 from row 2 onward (line search not used).
- seed:
  Random seed used to generate this one-seed trajectory dataset.
- n:
  Number of nodes.
- k:
  Latent dimension.
- family:
  Edge distribution family: bernoulli / poisson / gaussian.
- method:
  GD method: fixed or linesearch.
- eta0:
  Baseline step size computed from RA-SVT initialization.
- eta_mult:
  Step-size multiplier (here in {10, 5} for this figure script).
- eta_label:
  String label of multiplier, e.g., "10*eta0", "5*eta0".
- eta_init:
  Initial step size used by this case, eta_init = eta_mult * eta0.
- eta_times_n:
  n-scaled step size, eta_init * n.
- case:
  Case label combining method and eta_label, e.g., "linesearch | 5*eta0".

Notes
- For fixed-step runs, bt is 0 from row 2 onward, eta_used = eta_init, ls_ok = 1.
- For line-search runs, bt/eta_used/ls_ok reflect actual backtracking outcomes.

================================================================================
summary.csv: one row per (family, method, eta_mult) case
================================================================================
Columns and meanings
- seed:
  Random seed used in this one-seed experiment.
- n:
  Number of nodes.
- k:
  Latent dimension.
- family:
  Edge distribution family.
- method:
  GD method: fixed or linesearch.
- eta0:
  Baseline step size from RA-SVT initialization.
- eta_mult:
  Step-size multiplier in this figure script (10 or 5).
- eta_label:
  String label for eta_mult.
- eta_init:
  Initial step size, eta_mult * eta0.
- eta_times_n:
  n-scaled step size.
- eps_grad:
  Stopping threshold on max(abs(grad_Y)).
- max_gd_steps:
  Maximum GD iteration cap.
- beta_ls:
  Line-search backtracking shrink factor.
  Not meaningful for fixed-step GD (only relevant to line-search runs).
- c_ls:
  Line-search condition constant.
  Not meaningful for fixed-step GD (only relevant to line-search runs).
- Rprime:
  Backtracking cap per GD iteration.
  0 for fixed-step cases; positive cap for line-search cases.
  Not meaningful for fixed-step GD beyond being recorded as 0.
- num_steps:
  Number of GD iterations executed in this case.
- converged:
  Convergence indicator (1/0), based on eps_grad threshold.
- exploded:
  Explosion indicator (1/0), non-finite score encountered.
- error_message:
  Error text if any runtime issue occurs; empty otherwise.
- bt_mean:
  Mean backtracking count across GD iterations (NA for fixed-step).
  Not meaningful for fixed-step GD.
- bt_max:
  Maximum backtracking count over GD iterations (NA for fixed-step).
  Not meaningful for fixed-step GD.
- bt_frac0:
  Fraction of GD iterations with zero backtracking (NA for fixed-step).
  Not meaningful for fixed-step GD.
- ls_fail_frac:
  Fraction of GD iterations where line-search conditions still fail after
  reaching the backtracking cap Rprime (NA for fixed-step).
  Not meaningful for fixed-step GD.
