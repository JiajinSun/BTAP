### Script overview:
### (1) This file runs the 8-case study:
###     method in {fixed, linesearch} x eta_mult in {10, 5, 1, 0.2},
###     with eta_init = eta_mult * eta0 (eta0 from RA-SVT initialization).
### (2) Main worker is run_one_rep(rep_id):
###     it generates one dataset (Z, alpha, A), runs all 8 cases on that dataset,
###     and returns per-case algorithm outputs.
### (3) Parallel execution is launched by:
###     results_list <- foreach(rep_id = seq_len(reps), ...) %dopar% run_one_rep(rep_id)
### (4) User inputs are set in the section starting with:
###     "### User-configurable inputs" (family, n, reps, cores, and other defaults below).
### (5) Output files (saved to `simulation/results for GD behavior`) for each (family, n):
###     - fs_vs_ls_4etas_{family}_n{n}_run.csv:
###       one row per successful (rep, method, eta_mult) case with diagnostics.
###     - fs_vs_ls_4etas_{family}_n{n}_summary.csv:
###       aggregated means/medians/rates grouped by method and eta_mult.
###
### load function file and packages

### Working directory requirement:
### Set working directory to BTAP before running this script.
### Example: getwd() should return the BTAP folder path.
current_dir <- normalizePath("simulation")
out_dir <- normalizePath(file.path(current_dir, "results for GD behavior"))
if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

source("functions_RASVT_PGDwBLS.R")

library(foreach)
library(doParallel)
library(dplyr)
library(Rmpfr)

### pass arguments with Rscript (for running on cluster with .sh files)
args <- commandArgs(trailingOnly = TRUE)
parse_arg <- function(key, default){
  hit <- grep(paste0("^--", key, "="), args, value = TRUE)
  if(length(hit) == 0) return(default)
  sub(paste0("^--", key, "="), "", hit[1])
}


### User-configurable inputs
### family, n, reps, cores are passed through parse_arg (from, e.g., .sh files);
### if not given, the default is "poisson", "500", "100", "10".\
#### for testing purpose (on a local computer), one may change `reps` to 10
### to run this with R or RStudio, just change the defaults below to user specified choices.
family <- tolower(parse_arg("family", "poisson"))
    ## family must be one of: bernoulli, poisson, gaussian
n <- as.integer(parse_arg("n", 500L))
reps <- as.integer(parse_arg("reps", 100L))  #### for testing purpose, one may change `reps` to 10
cores <- as.integer(parse_arg("cores", 10L)) ## number of cores for parallel computing

k <- 2L   ## latent dimension
eps_grad <- 1e-2   ## GD stop criteria: score_Y max < eps_grad
max_gd_steps <- 2000L   ## GD stop criteria: # of iteration > max_gd_steps
beta_ls <- 0.5   ## line search backtrack shrink parameter
c_ls <- 1   ## line search condition parameter
Rprime_line <- max(1L, as.integer(ceiling(log(n))))   ## line search backtrack max times parameter

eta_mults <- c(10, 5, 1, 0.2)  ## we study eta_init = 10 eta0, 5 eta0, eta0, eta0/5 

# Reproducible seeds: seed_i = 1000 * (i + 2) + 123, i = 1,...,reps
seeds <- 1000L * (seq_len(reps) + 2L) + 123L

### 3 helper functions for data generation: for Z, for alpha, and for A
### generate Z parameter
make_Z <- function(n, k, z_scale = 0.5){
  Z_keep <- matrix(NA_real_, nrow = 0, ncol = k)
  while(nrow(Z_keep) < n){
    Z_raw <- matrix(rnorm(n * k * 5), n * 5, k)
    keep <- apply(Z_raw, 1, function(x) sum(x > 2) + sum(x < -2)) == 0
    Z_keep <- rbind(Z_keep, Z_raw[keep, , drop = FALSE])
  }
  Z <- Z_keep[1:n, , drop = FALSE]
  Z <- Z - rep(1, n) %*% (t(rep(1, n)) %*% Z) / n
  z_scale * sqrt(n) * svd(Z)$u
}

### generate alpha parameter
make_alpha <- function(n){
  alpha_tilde <- runif(n, 1, 3)
  -alpha_tilde / sum(alpha_tilde)
}

### generate A from Theta
simulate_A <- function(Theta, family){
  n <- nrow(Theta)
  if(family == "bernoulli"){
    P <- 1 / (1 + exp(-Theta))
    A <- matrix(rbinom(n * n, size = 1, prob = as.vector(P)), n, n)
  } else if(family == "poisson"){
    Lambda <- exp(Theta)
    A <- matrix(rpois(n * n, lambda = as.vector(Lambda)), n, n)
  } else if(family == "gaussian"){
    A <- matrix(rnorm(n * n, mean = as.vector(Theta), sd = 1), n, n)
  } else {
    stop("Unsupported family in simulate_A.")
  }
  A[upper.tri(A)] <- t(A)[upper.tri(A)]
  diag(A) <- 0
  A
}




cat(sprintf("Running family=%s, n=%d, reps=%d, GDmax=%d, cores=%d\n", family, n, reps, max_gd_steps, cores))
cat(sprintf("Seeds: %s\n", paste(seeds, collapse = ", ")))

### function of running one replication (one seed), 
### containing 8 cases: (fixed step, line search) * (4 choices of eta_init)
run_one_rep <- function(rep_id){
  seed <- seeds[rep_id]
  set.seed(seed)
  ### data generation
  Z <- make_Z(n, k, z_scale = 0.5)
  alpha <- make_alpha(n)
  Theta <- outer(alpha, alpha, "+") + Z %*% t(Z)
  A <- simulate_A(Theta, family)
  
  ### calculating eta0 from RA-SVT initializer
  eta0 <- as.numeric(eta_cornodiv6_from_init(A, k, family))
  
  run_rows_rep <- list()   ## for saving algorithmic results
  cc <- 1L   ### 4 etas * (fixed, ls) = 8 cases
  for(mult in eta_mults){
    eta_init <- mult * eta0
    for(method in c("fixed", "linesearch")){
      ### running the algorithm of Range Adaptive-Singular Value Thresholding + Projected Gradient Descent with Backtrakcing Line Search
      out <- RASVT_PGDwBLS(
        A = A, k = k, family = family, eta_init = eta_init, method = method,
        eps_grad = eps_grad, max_gd_steps = max_gd_steps,
        beta_ls = beta_ls, c_ls = c_ls, Rprime_line = Rprime_line
      )
      
      ## Output row saved to fs_vs_ls_4etas_{family}_n{n}_run.csv:
      ## - rep, seed: replication index and RNG seed.
      ## - n, k, family: simulation setup.
      ## - method: "fixed" or "linesearch".
      ## - eta0: baseline eta from RA-SVT initializer; eta_mult in {10, 5, 1, 0.2}.
      ## - eta_init, eta_times_n: actual initial step and n-scaled step size.
      ## - eps_grad, max_gd_steps, beta_ls, c_ls: stopping/line-search hyperparameters.
      ## - Rprime in fixed-step GD: set to 0 because no backtracking line search is used.
      ## - Rprime in line-search GD: maximum number of backtracking reductions per GD iteration.
      ## - error: 1 means the optimizer returned a run-level error (case failed); 0 means normal return.
      ## - error_message: text from out$message when error=1; empty string otherwise.
      ## - elapsed_sec: wall-clock runtime for this case.
      ## - num_of_steps, converged, exploded: GD iteration outcomes.
      ## - score_Y_max_abs: max absolute entry of final grad_Y.
      ## - bt_mean: mean # of backtracks per GD iteration.
      ## - bt_max: maximum # of backtracks observed in one GD iteration.
      ## - bt_frac0: fraction of GD iterations with zero backtracking.
      ## - bt_frac_hit_cap: fraction of GD iterations that used exactly Rprime backtracks.
      ## - ls_fail_frac: fraction of GD iterations where line-search conditions were still unsatisfied at Rprime backtracks
      run_rows_rep[[cc]] <- data.frame(  ## for saving algorithmic results; cc is 1-8
        rep = rep_id,
        seed = seed,
        n = n,
        k = k,
        family = family,
        method = method,
        eta0 = eta0,
        eta_mult = mult,
        eta_init = eta_init,
        eta_times_n = eta_init * n,
        eps_grad = eps_grad,
        max_gd_steps = max_gd_steps,
        beta_ls = beta_ls,
        c_ls = c_ls,
        Rprime = if(method == "linesearch") Rprime_line else 0L,
        error = as.integer(out$error),
        error_message = if(isTRUE(out$error)) out$message else "",
        elapsed_sec = out$elapsed,
        num_of_steps = if(isTRUE(out$error)) NA_integer_ else out$num_of_steps,
        converged = if(isTRUE(out$error)) 0L else out$converged,
        exploded = if(isTRUE(out$error)) NA_integer_ else out$exploded,
        score_Y_max_abs = if(isTRUE(out$error)) NA_real_ else out$score_Y_max_abs,
        bt_mean = if(isTRUE(out$error)) NA_real_ else out$bt_mean,
        bt_max = if(isTRUE(out$error)) NA_real_ else out$bt_max,
        bt_frac0 = if(isTRUE(out$error)) NA_real_ else out$bt_frac0,
        bt_frac_hit_cap = if(isTRUE(out$error)) NA_real_ else out$bt_frac_hit_cap,
        ls_fail_frac = if(isTRUE(out$error)) NA_real_ else out$ls_fail_frac,
        stringsAsFactors = FALSE
      )
      
      cc <- cc + 1L
    }
  }
  
  do.call(rbind, run_rows_rep)
}


### parallel cluster setup
cl <- makeCluster(cores)
registerDoParallel(cl,cores = cores)

### run the simulation with foreach parallel
results_list <- foreach(rep_id = seq_len(reps), .inorder = TRUE, .errorhandling = "pass") %dopar% run_one_rep(rep_id)

stopCluster(cl)

failed_ids <- which(vapply(results_list, inherits, logical(1), what = "error"))
if(length(failed_ids) > 0){
  warning(sprintf("Parallel worker errors at rep IDs: %s", paste(failed_ids, collapse = ", ")))
}
ok_ids <- setdiff(seq_along(results_list), failed_ids)


# run_df: one row per successful (replication, method, eta_mult) case,
# containing algorithm diagnostics.
run_df <- do.call(rbind, results_list[ok_ids])

# summary_df: aggregate of run_df by method and eta_mult
# (means/medians/rates across successful reps for each case).
summary_df <- run_df |>
  dplyr::group_by(method, eta_mult) |>
  dplyr::summarise(
    reps = dplyr::n(),                                      # number of rows in each (method, eta_mult) cell
    error_rate = mean(error),                               # fraction with run-level error==1
    conv_rate = mean(converged, na.rm = TRUE),              # fraction GD converged by eps_grad before R steps
    explode_rate = mean(exploded, na.rm = TRUE),            # fraction marked GD explosion
    mean_elapsed = mean(elapsed_sec, na.rm = TRUE),         # mean runtime (seconds)
    median_elapsed = stats::median(elapsed_sec, na.rm = TRUE), # median runtime (seconds)
    mean_steps = mean(num_of_steps, na.rm = TRUE),          # mean GD iterations used
    mean_score_Y = mean(score_Y_max_abs, na.rm = TRUE),     # mean terminal max(abs(grad_Y))
    mean_bt_frac0 = mean(bt_frac0, na.rm = TRUE),           # mean fraction of GD steps with 0 backtracks
    mean_bt_max = mean(bt_max, na.rm = TRUE),               # mean maximum backtracks used within a run
    mean_bt_hit_cap = mean(bt_frac_hit_cap, na.rm = TRUE),  # mean fraction of steps hitting Rprime cap
    mean_ls_fail_frac = mean(ls_fail_frac, na.rm = TRUE),   # mean fraction of steps failing line search at Rprime cap
    .groups = "drop"
  )

run_out <- file.path(out_dir, sprintf("fs_vs_ls_4etas_%s_n%d_run.csv", family, n))
sum_out <- file.path(out_dir, sprintf("fs_vs_ls_4etas_%s_n%d_summary.csv", family, n))

write.csv(run_df, run_out, row.names = FALSE)
write.csv(summary_df, sum_out, row.names = FALSE)
