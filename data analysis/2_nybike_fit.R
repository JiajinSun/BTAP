### Working directory requirement:
### Set working directory to BTAP before running this script.
### Example: getwd() should return the BTAP folder path.

### Script overview:
### (1) Load functions for initialization + GD fitting
### (2) Read preprocessed NY bike data from `data analysis/results/nybike_data.rds`.
### (3) For each selected hour (08:00 and 18:00): run the initialization + GD fitting (RASVT_PGDwBLS) with line search on the symmetric adjacency.
### (4) Store full fitted objects for both hours in `results`.
### (5) Build a compact `summary_df` with convergence, runtime, score, backtracking,
###     and network-size/count diagnostics for quick inspection.
### (6) Save the full fit results as `data analysis/results/nybike_fit.rds`.

source("functions_RASVT_PGDwBLS.R")

# -------------------------
# User-configurable inputs
# -------------------------
family <- "poisson"        # Edge distribution family used by the fitting routine.
k <- 2L                     # Latent dimension.
eps_grad <- 1e-2            # Stopping threshold on gradient score.
max_gd_steps <- 3000L       # Maximum GD iterations.
beta_ls <- 0.5              # Backtracking shrink factor.
c_ls <- 1                   # Line-search condition constant.
eta_mult <- 20              # eta_init multiplier relative to eta0.

# Input/output file paths.
filtered_data_path <- file.path("data analysis/results", "nybike_data.rds")
results_out <- file.path("data analysis/results", "nybike_fit.rds")

# -------------------------
# Load preprocessed network data
# -------------------------
# `filtered` is produced by 1_preprocess.R and contains:
# selected hours, retained node indices/labels, and directed/symmetric matrices.
filtered <- readRDS(filtered_data_path)

# Number of retained nodes (same for both hours after preprocessing filter).
n <- nrow(filtered$A_symmetric_hour1)
Rprime_line <- max(1L, as.integer(ceiling(log(n))))   ## we set this here because we need "n"

# Pretty hour label for summary display.
hour_label <- function(hour){
  sprintf("%02d:00-%02d:00", hour, (hour + 1L) %% 24L)
}

# -------------------------
# Helper: fit one hour network
# -------------------------
# Input:
# - A: symmetric adjacency used for model fitting.
# - hour: integer hour index used for labels.
# Output:
# - list containing metadata, input networks, and `fit` returned by RASVT_PGDwBLS.
fit_one_hour <- function(A, hour){
  # Compute baseline step size from initialization, then scale it.
  eta0 <- as.numeric(eta_cornodiv6_from_init(A, k, family))
  eta_init <- eta_mult * eta0

  # Main optimization call: line-search GD with RASVT initialization.
  fit <- RASVT_PGDwBLS(
    A = A, k = k, family = family, eta_init = eta_init, method = "linesearch",
    eps_grad = eps_grad, max_gd_steps = max_gd_steps,
    beta_ls = beta_ls, c_ls = c_ls, Rprime_line = Rprime_line
  )

  list(
    hour = hour,
    hour_label = hour_label(hour),
    eta0 = eta0,
    eta_init = eta_init,
    A_symmetric = A,
    fit = fit
  )
}

# -------------------------
# Fit both target hours
# -------------------------
results <- list(
  hour_08 = fit_one_hour(filtered$A_symmetric_hour1, filtered$hour1),
  hour_18 = fit_one_hour(filtered$A_symmetric_hour2, filtered$hour2)
)

# -------------------------
# Build compact summary table
# -------------------------
# One row per hour with key diagnostics pulled from each fit object.
summary_df <- do.call(rbind, lapply(results, function(result){
  fit <- result$fit
  data.frame(
    hour = result$hour,
    hour_label = result$hour_label,
    eta0 = result$eta0,
    eta_multiplier = eta_mult,
    eta_init = result$eta_init,
    n_nodes = nrow(result$A_symmetric),
    converged = fit$converged,
    error = as.integer(fit$error),
    error_message = if(is.null(fit$message)) "" else fit$message,
    elapsed_sec = fit$elapsed,
    num_of_steps = fit$num_of_steps,
    score_Y_max_abs = fit$score_Y_max_abs,
    bt_mean = fit$bt_mean,
    bt_max = fit$bt_max,
    bt_frac0 = fit$bt_frac0,
    bt_frac_hit_cap = fit$bt_frac_hit_cap,
    ls_fail_frac = fit$ls_fail_frac,
    total_symmetric_count = sum(result$A_symmetric),
    stringsAsFactors = FALSE
  )
}))

# -------------------------
# Save full outputs and print summary
# -------------------------
# `results` keeps all model objects (for downstream plotting/analysis).
saveRDS(results, results_out)

# Quick console check.
print(summary_df)
