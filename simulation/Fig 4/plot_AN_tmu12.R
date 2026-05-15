### Script overview:
### (1) Locate eta0-only asymptotic-record RDS files in
###     `../results for asymptotic distribution`.
### (2) Load records and, for each replication, recompute
###     `t(mu_12) = [mu(hat Theta_12) - mu(Theta_12^*)] / sqrt(var_hat(mu_12))`
###     using the hat-variance form:
###     `var_hat(mu_12) = {mu'(hat Theta_12)}^2
###      [w_2^T Sigma_1(hat Y)^{-1} w_2 + w_1^T Sigma_2(hat Y)^{-1} w_1]`.
### (3) Assemble one long table with identifiers
###     (`family`, `n`, `file_base`, `rep`, `seed`) and `tmu12_hat`.
### (4) For each `(family, n)`, keep finite `tmu12_hat` values and create:
###     - QQ plot versus N(0,1)
###     - Histogram with N(0,1) density overlay
### (5) Output files (saved in `out_dir`):
###     - `qq_t_muTheta12_{family}_n{n}.pdf`
###     - `hist_t_muTheta12_{family}_n{n}.pdf`

## some global plotting parameters
textsize <- 27

library(ggplot2)

# directories
### Working directory requirement:
### Set working directory to BTAP before running this script.
### Example: getwd() should return the BTAP folder path.
current_dir <- normalizePath("simulation/Fig 4")
# Input folder
data_dir <- normalizePath("simulation/results for asymptotic distribution")
# Output folder for generated figures.
out_dir <- current_dir

# Collect all eta0-only asymptotic-record files.
rds_files <- sort(Sys.glob(file.path(data_dir, "eta0_only_*_n*_asymptotic_records.rds")))
if(length(rds_files) == 0){
  stop("No eta0-only asymptotic record files found in: ", data_dir)
}


### Reading all the RDS files being used
### (1) Parse each RDS filename into metadata `meta_list`, where each element is
###     list(family, n, file, file_base):
###       - file     : full file path, used by readRDS(file),
###       - file_base: basename(file), saved for source-tracking in output rows.
### (2) recs_by_file is a list aligned with meta_list, Each element recs_by_file[[j]] is the object read from one RDS file.
###     recs_by_file[[j]] is the per-replication record list read from meta_list[[j]]$file.

# Function: parse_meta
# Purpose:
#   Parse family and n from one RDS filename.
# Input:
#   path: one RDS filepath.
# Output:
#   list(family, n, file, file_base) or NULL if pattern mismatch.
parse_meta <- function(path){
  bn <- basename(path)
  m <- regexec("^eta0_only_(bernoulli|poisson|gaussian)_n([0-9]+)_asymptotic_records\\.rds$", bn)
  g <- regmatches(bn, m)[[1]]
  if(length(g) == 0) return(NULL)
  list(family = g[2], n = as.integer(g[3]), file = path, file_base = bn)
}

meta_list <- Filter(Negate(is.null), lapply(rds_files, parse_meta))

recs_by_file <- lapply(meta_list, function(meta) readRDS(meta$file))


### compute the t(mu12) stat for each replication
# t(mu_12) computation roadmap (used inside one_stat):
# (A) Parse one replication's err_vec into (Z_true, alpha_true, Z_est, alpha_est)
#     via parse_errvec_basic().
# (B) Build Theta_12^* and hat Theta_12, then compute mu(Theta) and mu'(hat Theta)
#     via mu_fun() and mu_prime_fun().
#     mu(Theta) for star and hat are used in the numerator of the t(mu12). 
#     mu'(hat Theta) is used in the variance formula, for the denominator of t(mu12)
# (C) Build Sigma_1(hat Y) and Sigma_2(hat Y) with Sigma_i_family(), where
#     neg_ell_second_fun() provides the family-specific -ell''(theta).
# (D) Form var_hat(mu_12) from Sigma_1(hat Y)^{-1}, Sigma_2(hat Y)^{-1},
#     w_1, w_2, and mu'(hat Theta_12), then return
#     [mu(hat Theta_12) - mu(Theta_12^*)] / sqrt(var_hat(mu_12)).
# The helper functions below exist to keep these steps explicit and readable.

# Function: parse_errvec_basic
# Purpose:
#   Extract true/estimated (Z, alpha) blocks from err_vec prefix.
### these (Z,alpha) hat and star will be used for calculating the t(mu12).
parse_errvec_basic <- function(err_vec, n, k){
  # Expected prefix:
  # c(t(Z_true), alpha_true, t(Z_est), alpha_est, ...)
  need_len <- 2 * (n * k + n)
  if(length(err_vec) < need_len){
    return(NULL)
  }
  idx <- 1L
  zt <- err_vec[idx:(idx + n * k - 1L)]; idx <- idx + n * k
  at <- err_vec[idx:(idx + n - 1L)];     idx <- idx + n
  ze <- err_vec[idx:(idx + n * k - 1L)]; idx <- idx + n * k
  ae <- err_vec[idx:(idx + n - 1L)]
  list(
    Z_true = t(matrix(zt, nrow = k, ncol = n)),
    alpha_true = as.vector(at),
    Z_est = t(matrix(ze, nrow = k, ncol = n)),
    alpha_est = as.vector(ae)
  )
}

# Function: mu_fun
# Purpose:
#   Family-specific mean/link mapping mu(theta).
mu_fun <- function(theta, family){
  if(family == "poisson") return(exp(theta))
  if(family == "bernoulli") return(plogis(theta))
  if(family == "gaussian") return(theta)
  stop("Unsupported family: ", family)
}

# Function: mu_prime_fun
# Purpose:
#   Family-specific derivative mu'(theta).
mu_prime_fun <- function(theta, family){
  if(family == "poisson") return(exp(theta))
  if(family == "bernoulli"){
    p <- plogis(theta)
    return(p * (1 - p))
  }
  if(family == "gaussian") return(rep(1, length(theta)))
  stop("Unsupported family: ", family)
}

# Function: neg_ell_second_fun
# Purpose:
#   Return -ell''(theta; x) evaluated at theta (family-specific form).
### this is used in calculating Sigma_i
neg_ell_second_fun <- function(theta, family){
  if(family == "poisson") return(exp(theta))
  if(family == "bernoulli"){
    p <- plogis(theta)
    return(p * (1 - p))
  }
  if(family == "gaussian") return(rep(1, length(theta)))
  stop("Unsupported family: ", family)
}

# Function: Sigma_i_family
# Purpose:
#   Compute Sigma_i(Y) from estimated parameters for one node i.
Sigma_i_family <- function(Z, alpha, i, family){
  # Sigma_i(hat Y) = -sum_{j != i} ell''(Theta_ij) w_j w_j^T
  # where w_j = [z_j^T, 1]^T.
  n <- nrow(Z)
  W <- cbind(Z, rep(1, n))  # n x (k+1)
  theta_i <- alpha + alpha[i] + as.vector(Z %*% Z[i, ])
  w <- neg_ell_second_fun(theta_i, family)  # length n
  S <- crossprod(W, W * w)
  # no-self-loop version (match existing simulation code path)
  S <- S - w[i] * tcrossprod(W[i, ])
  S
}

# Function: one_stat
# Purpose:
#   Compute one replication's t(mu_12) statistic from one record entry.
one_stat <- function(rec, family, n, k){
  ev <- rec$err_vec
  if(is.null(ev) || !is.numeric(ev)) return(NA_real_)
  parsed <- parse_errvec_basic(ev, n, k)
  if(is.null(parsed)) return(NA_real_)

  Zt <- parsed$Z_true; at <- parsed$alpha_true
  Ze <- parsed$Z_est;  ae <- parsed$alpha_est
  if(any(!is.finite(Zt)) || any(!is.finite(at)) || any(!is.finite(Ze)) || any(!is.finite(ae))){
    return(NA_real_)
  }

  i <- 1L; j <- 2L
  th_true_12 <- at[i] + at[j] + sum(Zt[i, ] * Zt[j, ])
  th_hat_12  <- ae[i] + ae[j] + sum(Ze[i, ] * Ze[j, ])
  mu_true_12 <- mu_fun(th_true_12, family)
  mu_hat_12  <- mu_fun(th_hat_12, family)
  mu_prime_hat_12 <- mu_prime_fun(th_hat_12, family)

  S1 <- Sigma_i_family(Ze, ae, i, family)
  S2 <- Sigma_i_family(Ze, ae, j, family)
  S1_inv <- tryCatch(solve(S1), error = function(e) NULL)
  S2_inv <- tryCatch(solve(S2), error = function(e) NULL)
  if(is.null(S1_inv) || is.null(S2_inv)) return(NA_real_)

  w1 <- c(Ze[i, ], 1)
  w2 <- c(Ze[j, ], 1)
  middle <- as.numeric(t(w2) %*% S1_inv %*% w2 + t(w1) %*% S2_inv %*% w1)
  var_hat <- as.numeric(mu_prime_hat_12^2) * middle
  if(!is.finite(var_hat) || var_hat <= 0) return(NA_real_)

  as.numeric(mu_hat_12 - mu_true_12) / sqrt(var_hat)
}



### Data preprocessing (before plotting):
### (1) Loop through all replications in all RDS files and append the calculated t(mu12) per replication
###     into `all_rows`:
###     family, n, file_base, rep, seed, t(mu12).
### (2) Bind the list into one long data frame `df` via do.call(rbind, all_rows).
all_rows <- list()
for(j in seq_along(meta_list)){
  meta <- meta_list[[j]]
  recs <- recs_by_file[[j]]
  if(length(recs) == 0L) next
  for(i in seq_along(recs)){
    rec <- recs[[i]]
    n <- if(!is.null(rec$n)) as.integer(rec$n) else meta$n
    k <- if(!is.null(rec$k)) as.integer(rec$k) else NA_integer_
    tmu <- if(!is.finite(k) || k <= 0) NA_real_ else one_stat(rec, meta$family, n, k)
    all_rows[[length(all_rows) + 1L]] <- data.frame(
      family = meta$family,
      n = meta$n,
      file_base = meta$file_base,
      rep = if(!is.null(rec$rep)) as.integer(rec$rep) else i,
      seed = if(!is.null(rec$seed)) as.integer(rec$seed) else NA_integer_,
      tmu12_hat = as.numeric(tmu),
      stringsAsFactors = FALSE
    )
  }
}
if(length(all_rows) == 0L){
  stop("All matched RDS files are empty.")
}

df <- do.call(rbind, all_rows)
df <- df[is.finite(df$tmu12_hat), , drop = FALSE]
if(nrow(df) == 0){
  stop("No finite tmu12_hat values were computed.")
}

# Function: make_qq
# Purpose:
#   Construct QQ plot against N(0,1) for one numeric vector.
make_qq <- function(x){
  ggplot(data.frame(x = x), aes(sample = x)) +
    stat_qq(distribution = qnorm) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    coord_equal(xlim = c(-3.5, 3.5), ylim = c(-3.5, 3.5), expand = FALSE) +
    labs(
      title = NULL,
      x = "Theoretical Quantiles",
      y = "Sample Quantiles"
    ) +
    theme_minimal(base_size = textsize) +
    theme(
      plot.title = element_text(size = textsize, face = "bold", hjust = 0.5),
      axis.title = element_text(size = textsize),
      axis.text = element_text(size = textsize),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "gray92", linewidth = 0.25),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7),
      plot.margin = margin(2,2,0,2)
    )
}

# Function: make_hist
# Purpose:
#   Construct histogram with N(0,1) density overlay for one numeric vector.
make_hist <- function(x){
  ggplot(data.frame(x = x), aes(x = x)) +
    geom_histogram(aes(y = after_stat(density)), bins = 30, color = "white") +
    stat_function(fun = dnorm, args = list(mean = 0, sd = 1), linewidth = 1) +
    coord_cartesian(xlim = c(-3.5, 3.5)) +
    labs(
      title = NULL,
      x = expression(hat(t)(mu[12])),
      y = "Density"
    ) +
    theme_minimal(base_size = textsize) +
    theme(
      plot.title = element_text(size = textsize, face = "bold", hjust = 0.5),
      axis.title = element_text(size = textsize),
      axis.text = element_text(size = textsize),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "gray92", linewidth = 0.25),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7),
      plot.margin = margin(2, 2, 2, 2)
    )
}

# Generate all plots for all families
families <- c("bernoulli", "poisson", "gaussian")
for(fam in families){
  sub <- df[df$family == fam, , drop = FALSE]
  ns <- sort(unique(sub$n))
  if(length(ns) == 0) next
  for(nv in ns){
    s <- sub[sub$n == nv, , drop = FALSE]
    x <- s$tmu12_hat
    seed_vec <- stats::na.omit(s$seed)
    n_seed <- length(unique(seed_vec))
    cat(sprintf("[t_muTheta12] family=%s n=%d: n_rep=%d, unique_seed=%d\\n",
                fam, nv, length(x), n_seed))
    qq_out <- file.path(out_dir, sprintf("qq_t_muTheta12_%s_n%d.pdf", fam, nv))
    hist_out <- file.path(out_dir, sprintf("hist_t_muTheta12_%s_n%d.pdf", fam, nv))
    ggsave(
      qq_out, make_qq(x),
      width = 5.6, height = 5.5, units = "in", device = "pdf", useDingbats = FALSE
    )
    ggsave(
      hist_out, make_hist(x),
      width = 5.6, height = 5.6, units = "in", device = "pdf", useDingbats = FALSE
    )
  }
}

cat("Saved t(mu12) asymptotic-normality outputs to:\n")
cat(out_dir, "\n")
