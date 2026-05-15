### Script overview:
### (1) This file runs the eta0-only simulation case:
###     one method ("linesearch") with eta_init = eta0 from RA-SVT initialization.
### (2) Main worker is run_one_rep(rep_id):
###     it generates one dataset (Z, alpha, A), runs RASVT_PGDwBLS once,
###     and returns per-rep algorithm outputs and records.
### (3) Parallel execution is launched by:
###     results_list <- foreach(rep_id = seq_len(reps), ...) %dopar% run_one_rep(rep_id)
### (4) User inputs are set in the section starting with:
###     "### User-configurable inputs" (family, n, reps, cores, and other defaults below).
### (5) Output files (saved to `simulation/results for asymptotic distribution`) for each (family, n):
###     - eta0_only_{family}_n{n}_asymptotic_records.rds:
###       list of per-rep asymptotic records (metadata + err_vec).
###
### load function file and packages
### Working directory requirement:
### Set working directory to BTAP before running this script.
### Example: getwd() should return the BTAP folder path.
current_dir <- normalizePath("simulation")
out_dir <- normalizePath(file.path(current_dir, "results for asymptotic distribution"))
if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

source("functions_RASVT_PGDwBLS.R")

library(foreach)
library(doParallel)
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
### if not given, the default is "poisson", "500", "200", "10".
#### for testing purpose (on a local computer), one may change `reps` to 10
### to run this with R or RStudio, just change the defaults below to user specified choices.
family <- tolower(parse_arg("family", "poisson"))
n <- as.integer(parse_arg("n", 500L))
reps <- as.integer(parse_arg("reps", 200L))  #### for testing purpose, one may change `reps` to 10
cores <- as.integer(parse_arg("cores", 10L)) ## number of cores for parallel computing

k <- 2L   ## latent dimension
eps_grad <- 1e-2   ## GD stop criteria: score_Y max < eps_grad
max_gd_steps <- 2000L   ## GD stop criteria: # of iteration > max_gd_steps
beta_ls <- 0.5   ## line search backtrack shrink parameter
c_ls <- 1   ## line search condition parameter
Rprime_line <- max(1L, as.integer(ceiling(log(n))))   ## line search backtrack max times parameter


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

Np <- n * k + n
err_vec_len <- 3 * Np + 1 + 2 + 1 + 1 + k * k + 2   ## the length of err_vec, with t11_hat/t11_star appended at the end
### function of running one replication (one seed) with line search GD with eta_init = eta0
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
  
  eta_init <- eta0
  method <- "linesearch"
  ### running the algorithm of Range Adaptive-Singular Value Thresholding + Projected Gradient Descent with Backtrakcing Line Search
  out <- RASVT_PGDwBLS(
    A = A, k = k, family = family, eta_init = eta_init, method = method,
    eps_grad = eps_grad, max_gd_steps = max_gd_steps,
    beta_ls = beta_ls, c_ls = c_ls, Rprime_line = Rprime_line
  )
  
  if(isTRUE(out$error)){
    rec <- rep(NA_real_, err_vec_len)
  } else {
    rec <- build_asymptotic_record(
      Z_true = Z,
      alpha_true = alpha,
      est = out,
      family = family
    )
  }
  ## Per-run asymptotic record saved into eta0_only_{family}_n{n}_asymptotic_records.rds:
  ## - rep, seed, method, eta_mult, eta_init, family, n, k: metadata.
  ## - err_vec: concatenated vector from build_asymptotic_record(), ordered as
  ##   c(
  ##     c(t(Z_true)),                     # nk
  ##     c(alpha_true),                    # n
  ##     c(t(Z_hat)),                      # nk
  ##     c(alpha_hat),                     # n
  ##     c(t(Z_hat_aligned - Z_true)),     # nk
  ##     c(alpha_hat - alpha_true),        # n
  ##     score_Y_max_abs,                  # 1
  ##     score_Y_max_which,                # 2 (max-gradient location index pair in Y)
  ##     score_Y_abs_mean,                 # 1
  ##     num_of_steps,                     # 1
  ##     c(O),                             # k^2, Procrustes rotation matrix
  ##     t11_hat,                          # 1
  ##     t11_star                          # 1
  ##   ).
  rec_out <- list(  ## for saving statistical results
    rep = rep_id,
    seed = seed,
    method = method,
    eta_mult = 1,
    eta_init = eta_init,
    family = family,
    n = n,
    k = k,
    err_vec = rec
  )
  
  rec_out
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

# combine the results from all reps
# asymp_records: list (length = number of successful reps) of per-rep asymptotic records;
# each element stores rep/seed metadata and err_vec (whose last two entries are t11_hat and t11_star).
asymp_records <- results_list[ok_ids]
for(ii in seq_along(asymp_records)) {
asymp_records[[ii]]$record_id <- ii
}
asymp_out <- file.path(out_dir, sprintf("eta0_only_%s_n%d_asymptotic_records.rds", family, n))

saveRDS(asymp_records, asymp_out)
