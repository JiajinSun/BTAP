### Script overview:
### (1) This script generates one-seed, per-iteration trajectories for 4 cases per family:
###     method in {fixed, linesearch} x eta_mult in {10, 5}.
### (2) It mirrors the main optimization logic in RASVT_PGDwBLS while recording
###     per-iteration diagnostics (log-likelihood, score max, backtracking).
### (3) It is intended for Fig 1 style trajectory plots.
### (4) User inputs are set in the section "User settings" below.
### (5) Output files (in out_dir) are:
###     - fs_vs_ls_{family}_n{n}_trajectory.csv:
###         one row per GD iteration per case (4 cases per family), with diagnostics
###         such as loglik, score_max_abs, bt, eta_used, ls_ok,
###         plus case metadata (family/method/eta_label/eta_init).
###     - fs_vs_ls_{family}_n{n}_summary.csv:
###         one row per case (4 rows per family), aggregating convergence/explosion,
###         num_steps, and line-search backtracking summaries (bt_mean/bt_max/bt_frac0/ls_fail_frac).

### Working directory requirement:
### Set working directory to BTAP before running this script.
### Example: getwd() should return the BTAP folder path.
current_dir <- normalizePath("simulation")
out_dir <- normalizePath(file.path(current_dir, "Fig 1"))  # Output directory for trajectory/summary CSVs.

source("functions_RASVT_PGDwBLS.R")

library(Rmpfr)

### User settings: edit these values directly before running.
n <- 1000L  # Number of nodes
k <- 2L  # Latent dimension (number of columns in Z).
seed <- 123L 
families <- "bernoulli"
# families <- c("gaussian","bernoulli","poisson")
eta_mults <- c(10, 5)  # Fixed to two eta multipliers for this script.
eps_grad <- 1e-2  # GD stopping threshold on max(abs(grad_Y))
max_gd_steps <- 2000L  # Maximum number of GD iterations
beta_ls <- 0.5  # Backtracking shrink factor
c_ls <- 1  # Line-search constant in global/individual acceptance inequalities.
Rprime_line <- max(1L, as.integer(ceiling(log(n))))  # Backtracking cap R' for line-search runs.
write_summary_csv <- TRUE  # Whether to also write the one-row-per-case summary CSV.


make_Z <- function(n, k, z_scale = 0.5) {
  Z_keep <- matrix(NA_real_, nrow = 0, ncol = k)
  while (nrow(Z_keep) < n) {
    Z_raw <- matrix(rnorm(n * k * 5), n * 5, k)
    keep <- apply(Z_raw, 1, function(x) sum(x > 2) + sum(x < -2)) == 0
    Z_keep <- base::rbind(Z_keep, Z_raw[keep, , drop = FALSE])
  }
  Z <- Z_keep[1:n, , drop = FALSE]
  Z <- Z - rep(1, n) %*% (t(rep(1, n)) %*% Z) / n
  z_scale * sqrt(n) * svd(Z)$u
}

make_alpha <- function(n) {
  alpha_tilde <- runif(n, -3, -1)
  -alpha_tilde / sum(alpha_tilde)
}

simulate_A <- function(Theta, family) {
  n <- nrow(Theta)
  if (family == "bernoulli") {
    P <- 1 / (1 + exp(-Theta))
    A <- matrix(rbinom(n * n, size = 1, prob = as.vector(P)), n, n)
  } else if (family == "poisson") {
    Lambda <- exp(Theta)
    A <- matrix(rpois(n * n, lambda = as.vector(Lambda)), n, n)
  } else if (family == "gaussian") {
    A <- matrix(rnorm(n * n, mean = as.vector(Theta), sd = 1), n, n)
  } else {
    stop("Unsupported family in simulate_A.")
  }
  A[upper.tri(A)] <- t(A)[upper.tri(A)]
  diag(A) <- 0
  A
}

### Function: run_with_trace_GD
### Purpose:
###   Run one optimization case (fixed or linesearch GD) and record trajectory by iteration.
### Inputs:
###   A, family, k, eta_init: model/data and step-size inputs.
###   method      : "fixed" or "linesearch".
###   eps_grad    : stopping threshold on max(abs(grad_Y)).
###   max_gd_steps: maximum GD iterations.
###   beta_ls, c_ls, Rprime_line: line-search parameters.
### Output:
###   list with:
###   - trace: per-iteration data.frame with columns
###       iter       : integer GD iteration index r (1,2,...).
###       loglik     : log-likelihood at current iterate Y^r.
###       score_max_abs: max(abs(grad_Y)) at Y^r, where grad_Y combines grad_Z and grad_alpha.
###       finite_score : 1 if score_max_abs is finite at Y^r, 0 otherwise.
###       bt         : number of backtracking reductions used to obtain the step from Y^{r-1} to Y^r.
###                    Hence bt is NA at r=1 (no previous step).
###       eta_used   : accepted step size used for the update Y^{r-1} -> Y^r; NA at r=1.
###       ls_ok      : line-search acceptance indicator for Y^{r-1} -> Y^r (1=accepted, 0=failed after cap);
###                    NA at r=1.
###            one row = diagnostics at iterate Y^r (before applying update from r to r+1).
###   - converged/exploded/error_message/num_steps.
run_with_trace_GD <- function(A, family, k, eta_init,
                                 method = c("fixed", "linesearch"),
                                 eps_grad = 1e-2, max_gd_steps = 2000L,
                                 beta_ls = 0.5, c_ls = 1, Rprime_line = NULL) {
  method <- match.arg(method)
  n <- nrow(A)
  one <- rep(1, n)
  J <- diag(n) - one %*% t(one) / n
  g <- get_link_functions(family)$g
  max_backtrack <- if (method == "linesearch") {
    if (is.null(Rprime_line)) max(1L, as.integer(ceiling(log(n)))) else as.integer(Rprime_line)
  } else {
    0L
  }

  init <- RASVT(A, k, family = family)
  Z <- J %*% init$Z0
  alpha <- as.vector(init$alpha0)

  # Pre-allocate per-iteration trajectory storage.
  loglik_col <- rep(NA_real_, max_gd_steps)
  score_col <- rep(NA_real_, max_gd_steps)
  finite_col <- rep(NA_integer_, max_gd_steps)
  # bt/eta/ls_ok on row r store transition info from iterate r-1 -> r (row 1 is NA).
  bt_col <- rep(NA_integer_, max_gd_steps)
  eta_col <- rep(NA_real_, max_gd_steps)
  ls_col <- rep(NA_integer_, max_gd_steps)
  converged <- FALSE
  exploded <- FALSE
  err_msg <- ""

  t <- 1L
  repeat {   ### the main GD loop, with explosion or convergence or max iteration checks for stop
    Theta <- alpha %*% t(one) + one %*% t(alpha) + Z %*% t(Z)
    P <- g(Theta)
    M <- (A - P) - diag(diag(A - P))
    gradZ <- M %*% Z
    gradA <- M %*% one

    score_z_max <- suppressWarnings(max(abs(gradZ)))
    score_alpha_max <- suppressWarnings(max(abs(gradA)))
    score_max <- suppressWarnings(max(score_z_max, score_alpha_max))
    finite_score <- is.finite(score_max)

    # negloglik_term() returns entrywise negative log-likelihood at current Theta0.
    nll <- negloglik_term(Theta, A, family)
    diag(nll) <- 0
    loglik <- -sum(nll[upper.tri(nll)])

    loglik_col[t] <- loglik                                 # upper-triangular log-likelihood
    score_col[t] <- score_max                               # max(abs(grad_Y))
    finite_col[t] <- as.integer(finite_score)

    if (!finite_score) {
      exploded <- TRUE
      err_msg <- "non-finite score"
      break
    }
    if (score_max <= eps_grad) {
      converged <- TRUE
      break
    }
    if (t >= max_gd_steps) break

    dZ <- J %*% gradZ
    dA <- gradA
    eta_ls <- eta_init

    if (max_backtrack > 0L) {
      Theta0 <- Theta
      NLL0 <- nll
      # Line-search uses two expansions:
      # (a) Global update along d(Y) = (d.Z, d.a):
      #     Theta_full(eta) = Theta0 + eta * B_full + eta * t(B_full) + eta^2 * C_full,
      #     B_full = d.a %*% t(1_n) + d.Z %*% t(Z0),   C_full = d.Z %*% t(d.Z).
      #     Global check:
      #       [L(Theta_full(eta)) - L(Theta0)] - c_ls * n * eta^2 * d_inner_global <= 0.
      # (b) Individual row-wise update (for all rows i):
      #     Theta_row(eta) = Theta0 + eta * B_row,
      #     where row i of Theta_row(eta) updates Theta0_i for evaluating L_i.
      #       B_row = grad0.a %*% t(1_n) + grad0.Z %*% t(Z0).
      #     Individual checks:
      #       [L_i(Theta_row(eta)) - L_i(Theta0)] - c_ls * n * eta^2 * d_inner_i[i] <= 0,
      #       for each i = 1,...,n.
      B_full <- dA %*% t(one) + dZ %*% t(Z)
      C_full <- dZ %*% t(dZ)
      B_row <- gradA %*% t(one) + gradZ %*% t(Z)
      d_inner_global <- sum((-gradZ) * dZ) + sum((-gradA) * dA)
      d_inner_i <- -rowSums(gradZ^2) - gradA^2

      # Evaluate global + rowwise line-search acceptance conditions.
      ls <- check_linesearch_conditions(
        eta_ls = eta_ls,
        Theta0 = Theta0,
        NLL0 = NLL0,
        B_full = B_full,
        C_full = C_full,
        B_row = B_row,
        A = A,
        family = family,
        c_ls = c_ls,
        n = n,
        d_inner_global = d_inner_global,
        d_inner_i = d_inner_i
      )
      line_search_ok <- ls$ok
      bt <- 0L
      while ((!line_search_ok) && (bt < max_backtrack)) {
        eta_ls <- beta_ls * eta_ls
        ls <- check_linesearch_conditions(
          eta_ls = eta_ls,
          Theta0 = Theta0,
          NLL0 = NLL0,
          B_full = B_full,
          C_full = C_full,
          B_row = B_row,
          A = A,
          family = family,
          c_ls = c_ls,
          n = n,
          d_inner_global = d_inner_global,
          d_inner_i = d_inner_i
        )
        line_search_ok <- ls$ok
        bt <- bt + 1L
      }
      # Store transition diagnostics on the next iterate row (r+1).
      if ((t + 1L) <= max_gd_steps) {
        bt_col[t + 1L] <- bt
        eta_col[t + 1L] <- eta_ls
        ls_col[t + 1L] <- as.integer(line_search_ok)
      }
    } else { ## max_backtrack = 0L
      if ((t + 1L) <= max_gd_steps) {
        bt_col[t + 1L] <- 0L
        eta_col[t + 1L] <- eta_ls
        ls_col[t + 1L] <- 1L
      }
    }

    # GD update
    Z <- Z + eta_ls * dZ
    alpha <- alpha + eta_ls * dA
    t <- t + 1L
  }

  tr <- data.frame(
    iter = seq_len(t),
    loglik = loglik_col[seq_len(t)],
    score_max_abs = score_col[seq_len(t)],
    finite_score = finite_col[seq_len(t)],
    bt = bt_col[seq_len(t)],
    eta_used = eta_col[seq_len(t)],
    ls_ok = ls_col[seq_len(t)],
    stringsAsFactors = FALSE
  )

  list(
    trace = tr,
    converged = converged,
    exploded = exploded,
    error_message = err_msg,
    num_steps = t
  )
}

# Generate one common truth (Z_true, alpha_true) and corresponding Theta_true.
set.seed(seed)
Z_true <- make_Z(n, k, z_scale = 0.5)
alpha_true <- make_alpha(n)
Theta_true <- alpha_true %*% t(rep(1, n)) + rep(1, n) %*% t(alpha_true) + Z_true %*% t(Z_true)

family_seed_offsets <- c(
  gaussian = 40L,
  bernoulli = 41L,
  poisson = 42L
)

# For each family:
# (1) generate one A from the same (Z_true, alpha_true),
# (2) compute eta0 from that A,
# (3) run 4 cases = 2 methods x 2 eta multipliers,
# (4) append per-case outputs to family-specific CSV files (streaming write).
for (family in families) {

  file_prefix <- paste0(
    "fs_vs_ls_", family,
    "_n", n
  )
  traj_out <- file.path(out_dir, paste0(file_prefix, "_trajectory.csv"))
  sum_out <- file.path(out_dir, paste0(file_prefix, "_summary.csv"))
  if (file.exists(traj_out)) file.remove(traj_out)
  if (write_summary_csv && file.exists(sum_out)) file.remove(sum_out)
  traj_header_written <- FALSE
  sum_header_written <- FALSE

  set.seed(seed + family_seed_offsets[[family]])
  A <- simulate_A(Theta_true, family)
  eta0 <- eta_cornodiv6_from_init(A, k, family)

  for (mult in eta_mults) {
    eta_init <- as.numeric(mult * eta0)
    eta_label <- paste0(mult, "*eta0")

    for (method in c("fixed","linesearch")) {
      cat(sprintf("[%s] family=%s method=%s eta=%s\n",
                  format(Sys.time(), "%Y-%m-%d %H:%M:%S"), family, method, eta_label))

      out <- run_with_trace_GD(
        A = A,
        family = family,
        k = k,
        eta_init = eta_init,
        method = method,
        eps_grad = eps_grad,
        max_gd_steps = max_gd_steps,
        beta_ls = beta_ls,
        c_ls = c_ls,
        Rprime_line = Rprime_line
      )

      tr <- out$trace
      # Add case metadata to each iteration row for downstream plotting/filtering.
      tr$seed <- seed
      tr$n <- n
      tr$k <- k
      tr$family <- family
      tr$method <- method
      tr$eta0 <- eta0
      tr$eta_mult <- mult
      tr$eta_label <- eta_label
      tr$eta_init <- eta_init
      tr$eta_times_n <- eta_init * n
      tr$case <- paste0(method, " | ", eta_label)
      if (!traj_header_written) {
        write.csv(tr, traj_out, row.names = FALSE)
        traj_header_written <- TRUE
      } else {
        write.table(tr, traj_out, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
      }

      bt_vals <- tr$bt[!is.na(tr$bt)]
      ls_ok_vals <- tr$ls_ok[!is.na(tr$ls_ok)]

      # One-row case summary for this (family, method, eta_mult).
      sum_row <- data.frame(
        seed = seed,
        n = n,
        k = k,
        family = family,
        method = method,
        eta0 = eta0,
        eta_mult = mult,  # Multiplier in {10, 5, 1, 0.2} applied to eta0.
        eta_label = eta_label,
        eta_init = eta_init,
        eta_times_n = eta_init * n,
        eps_grad = eps_grad,
        max_gd_steps = max_gd_steps,
        beta_ls = beta_ls,
        c_ls = c_ls,
        Rprime = if (method == "linesearch") Rprime_line else 0L,
        num_steps = out$num_steps,  # Total iterations executed for this case.
        converged = as.integer(out$converged),  # 1 if score criterion reached; 0 otherwise.
        exploded = as.integer(out$exploded),  # 1 if non-finite score encountered; 0 otherwise.
        error_message = out$error_message,
        bt_mean = if (method == "linesearch" && length(bt_vals) > 0) mean(bt_vals) else NA_real_,  # Average number of backtracks per GD step (NA for fixed-step).
        bt_max = if (method == "linesearch" && length(bt_vals) > 0) max(bt_vals) else NA_real_,  # Maximum backtracks taken on any GD step (NA for fixed-step).
        bt_frac0 = if (method == "linesearch" && length(bt_vals) > 0) mean(bt_vals == 0) else NA_real_,  # Fraction of steps accepted with zero backtracking (NA for fixed-step).
        ls_fail_frac = if (method == "linesearch" && length(ls_ok_vals) > 0) mean(ls_ok_vals == 0) else NA_real_,  # Fraction of steps still failing LS after cap R' (NA for fixed-step).
        stringsAsFactors = FALSE
      )
      if (write_summary_csv) {
        if (!sum_header_written) {
          write.csv(sum_row, sum_out, row.names = FALSE)
          sum_header_written <- TRUE
        } else {
          write.table(sum_row, sum_out, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
        }
      }

      rm(out, tr, sum_row)
      invisible(gc(FALSE))
    }
  }
  cat(sprintf("Saved trajectory: %s\n", traj_out))
  if (write_summary_csv) {
    cat(sprintf("Saved summary: %s\n", sum_out))
  }
  rm(A, eta0)
  invisible(gc(FALSE))
}
