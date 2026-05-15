This folder stores outputs from:
  simulation/linesearch_eta0.R

Working directory requirement
- Set working directory to BTAP before running the producing script.

Simulation scope
- Method: line-search GD only.
- Initial step size: eta_init = eta_0 (computed from RA-SVT initialization).
- Families: bernoulli, poisson, gaussian.
- Network sizes: n in {500, 1000, 2000, 4000}.

File naming pattern (one file per family x n):
  eta0_only_{family}_n{n}_asymptotic_records.rds

================================================================================
asymptotic_records.rds: list of per-rep records
================================================================================
Object structure
- Top level: list with one element per successful replication.
- Each element is a named list with fields:
  - rep:
    Replication index.
  - seed:
    Random seed used in that replication.
  - method:
    GD method; here fixed as "linesearch".
  - eta_mult:
    Multiplier relative to eta0; here fixed as 1.
  - eta_init:
    Initial step size used in this replication (equals eta0).
  - family:
    Edge distribution family: bernoulli / poisson / gaussian.
  - n:
    Number of nodes.
  - k:
    Latent dimension.
  - err_vec:
    Numeric vector storing true/estimated quantities and diagnostics
    (defined below).
  - record_id:
    Sequential record id within this file (1-based).

--------------------------------------------------------------------------------
err_vec definition (per replication)
--------------------------------------------------------------------------------
For each record, err_vec is:
  c(
    c(t(Z_true)),                     # nk
    c(alpha_true),                    # n
    c(t(Z_hat)),                      # nk
    c(alpha_hat),                     # n
    c(t(Z_hat_aligned - Z_true)),     # nk
    c(alpha_hat - alpha_true),        # n
    score_Y_max_abs,                  # 1
    score_Y_max_which,                # 2
    score_Y_abs_mean,                 # 1
    num_of_steps,                     # 1
    c(O),                             # k^2
    t11_hat,                          # 1
    t11_star                          # 1
  )

Length
- Let Np = n*k + n.
- Then length(err_vec) = 3*Np + 1 + 2 + 1 + 1 + k^2 + 2.
- In this project with k = 2, length(err_vec) = 7*n + 11.

Column-group meanings
- Z_true, alpha_true:
  Ground-truth latent positions and node effects used to generate data.
- Z_hat, alpha_hat:
  Final GD estimates before alignment.
- Z_hat_aligned - Z_true, alpha_hat - alpha_true:
  Estimation errors after Procrustes alignment of Z_hat to Z_true.
- score_Y_max_abs:
  Terminal max absolute entry of grad_Y.
- score_Y_max_which:
  2-entry index pair locating where |grad_Y| attains its maximum.
- score_Y_abs_mean:
  Mean absolute entry of grad_Y at termination.
- num_of_steps:
  Number of GD iterations executed.
- O:
  k x k Procrustes orthogonal matrix, vectorized by c(O).
- t11_hat:
  Standardized statistic using [DLii(hatY)^{-1}]_{1,1}.
- t11_star:
  Standardized statistic using [DLii(Y*)^{-1}]_{1,1}.

Notes
- This folder intentionally stores only RDS records (no run.csv/summary.csv).
- If a replication returns non-finite estimates, portions of err_vec can be NA.
