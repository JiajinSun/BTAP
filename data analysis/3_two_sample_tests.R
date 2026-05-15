### Working directory requirement:
### Set working directory to BTAP before running this script.
### Example: getwd() should return the BTAP folder path.

### Script overview:
### (1) Load fitted model outputs for two peak hours (08:00 and 18:00).
### (2) Build node-wise curvature caches (inverse local Hessian blocks) from each fit.
### (3) Conduct node-wise two-sample tests for alpha_i differences across hours.
### (4) Conduct pairwise two-sample tests for inner-product differences z_i^T z_j.
### (5) Apply BH correction and compute per-node rejection counts.
### (6) Save alpha tests, inner-product tests, and rejection-count CSV files.

### Output files produced by this script:
### (a) data analysis/results/nybike_08_vs_18_alpha.csv
###     - Node-wise test results for alpha_i differences across hours
###       (estimate difference, SE, z-stat, p-value, BH q-value, node metadata).
### (b) data analysis/results/nybike_08_vs_18_innerprod.csv
###     - Pairwise test results for <z_i, z_j> differences across hours
###       (estimate difference, SE, z-stat, p-value, BH q-value, pair metadata).
### (c) data analysis/results/nybike_innerprod_rejections_by_node.csv
###     - Node-level summary of pairwise rejections:
###       rejected_count, total_tests, rejected_frac for each retained node.

# -------------------------
# Paths, settings, Load input objects and helper functions
# -------------------------
results_path <- file.path("data analysis/results", "nybike_fit.rds")
filtered_data_path <- file.path("data analysis/results", "nybike_data.rds")
family <- "poisson"
output_prefix <- file.path("data analysis/results", "nybike_08_vs_18")
filtered <- readRDS(filtered_data_path)
borough <- filtered$borough
source("functions_RASVT_PGDwBLS.R")

# -------------------------
# Helper: node-wise curvature cache for one fitted hour
# -------------------------
# For each node i, cache inverse DLii block and extract:
# - alpha variance proxy: [inv_i]_{k+1,k+1}
# - z-block covariance proxy: [inv_i]_{1:k,1:k}
# Output (list of length 4):
# - inv: list of length n; inv[[i]] is the full inverse DLii matrix for node i.
# - alpha_var: numeric length n; alpha_var[i] = inv[[i]][k+1, k+1].
# - Vzz: list of length n; Vzz[[i]] is the k x k z-block covariance matrix.
# - ok: logical length n; TRUE if inversion/block extraction for node i succeeded.
build_curvature_cache <- function(Z, alpha, family){
  n <- nrow(Z)
  k <- ncol(Z)
  inv_list <- vector("list", n)
  alpha_var <- rep(NA_real_, n)
  Vzz_list <- vector("list", n)
  ok <- rep(FALSE, n)

  for(i in seq_len(n)){
    inv_i <- solve(DLii(Z, alpha, family, i))
    inv_list[[i]] <- inv_i
    if(!is.null(inv_i) && all(dim(inv_i) == c(k + 1L, k + 1L))){
      alpha_var[i] <- inv_i[k + 1L, k + 1L]
      Vzz_list[[i]] <- inv_i[seq_len(k), seq_len(k), drop = FALSE]
      ok[i] <- TRUE
    } else {
      Vzz_list[[i]] <- matrix(NA_real_, k, k)
    }
  }

  list(inv = inv_list, alpha_var = alpha_var, Vzz = Vzz_list, ok = ok)
}

# -------------------------
# Load fitted results (hour 08 vs hour 18)
# -------------------------
results_obj <- readRDS(results_path)
fit1 <- results_obj[[1]]$fit
fit2 <- results_obj[[2]]$fit
lab1 <- results_obj[[1]]$hour_label
lab2 <- results_obj[[2]]$hour_label
# cache inverse DLii block for each node i
cache1 <- build_curvature_cache(fit1$Z, as.vector(fit1$alpha), family)
cache2 <- build_curvature_cache(fit2$Z, as.vector(fit2$alpha), family)
n <- nrow(fit1$Z)

# -------------------------
# (A) Node-wise two-sample tests for alpha_i differences
## give the z-test and BH correction for alpha_i
# -------------------------
alpha_diff <- as.vector(fit1$alpha) - as.vector(fit2$alpha)
alpha_se2 <- cache1$alpha_var + cache2$alpha_var
alpha_se <- sqrt(alpha_se2)
alpha_z <- ifelse(is.finite(alpha_se) & alpha_se > 0, alpha_diff / alpha_se, NA_real_)  ## z-statistic
alpha_p <- ifelse(is.finite(alpha_z), 2 * pnorm(-abs(alpha_z)), NA_real_)  ## p-value

alpha_tests <- data.frame(
  filtered_node = seq_len(n),             # Node index in the filtered network (1..n).
  original_node = filtered$keep_idx,      # Corresponding node index in the original unfiltered network.
  borough = borough,                       # Borough label for the node.
  hour1 = lab1,                            # First hour label (e.g., "08:00-09:00").
  hour2 = lab2,                            # Second hour label (e.g., "18:00-19:00").
  alpha_hat_1 = as.vector(fit1$alpha),     # Estimated alpha_i at hour 1.
  alpha_hat_2 = as.vector(fit2$alpha),     # Estimated alpha_i at hour 2.
  diff = alpha_diff,                       # Difference alpha_hat_1 - alpha_hat_2.
  se = alpha_se,                           # Standard error of diff.
  z_stat = alpha_z,                        # Z statistic for testing diff = 0.
  p_value = alpha_p,                       # Two-sided p-value from z_stat.
  stringsAsFactors = FALSE
)
alpha_tests$q_value_bh <- p.adjust(alpha_tests$p_value, method = "BH")

# -------------------------
# (B) Pairwise tests for <z_i, z_j> differences
## give the z-test and BH correction for <z_i, z_j>
## output table meaning:
## - one row = one unordered node pair (i, j) with i < j;
## - columns store pair IDs/metadata plus estimates at both hours:
##     <z_i, z_j> at hour 1, <z_i, z_j> at hour 2, their difference,
##     standard error, z-statistic, p-value, and BH-adjusted q-value.
# -------------------------
pair_count <- n * (n - 1L) / 2L          # Total number of unordered node pairs (i, j), i < j.
pair_i <- integer(pair_count)            # Pair index i for each row in the pairwise output table.
pair_j <- integer(pair_count)            # Pair index j for each row in the pairwise output table.
prod_hat_1 <- numeric(pair_count)        # Estimated inner product <z_i, z_j> at hour 1 for each pair.
prod_hat_2 <- numeric(pair_count)        # Estimated inner product <z_i, z_j> at hour 2 for each pair.
prod_diff <- numeric(pair_count)         # Difference: prod_hat_1 - prod_hat_2 for each pair.
prod_se <- rep(NA_real_, pair_count)     # Standard error of prod_diff for each pair (NA if unavailable).
prod_z <- rep(NA_real_, pair_count)      # Z statistic for testing prod_diff = 0 (NA if SE unavailable).
prod_p <- rep(NA_real_, pair_count)      # Two-sided p-value corresponding to prod_z.

idx <- 1L
for(i in seq_len(n - 1L)){
  ## 1 is the hour1, 2 is the hour2
  zi1 <- fit1$Z[i, ]
  zi2 <- fit2$Z[i, ]
  Vii1 <- cache1$Vzz[[i]]  # Vzz[[i]] is the k x k z-block covariance matrix.
  Vii2 <- cache2$Vzz[[i]] 

  for(j in (i + 1L):n){
    zj1 <- fit1$Z[j, ]
    zj2 <- fit2$Z[j, ]
    Vjj1 <- cache1$Vzz[[j]]
    Vjj2 <- cache2$Vzz[[j]]

    pair_i[idx] <- i
    pair_j[idx] <- j

    g1 <- sum(zi1 * zj1)
    g2 <- sum(zi2 * zj2)
    prod_hat_1[idx] <- g1
    prod_hat_2[idx] <- g2
    prod_diff[idx] <- g1 - g2

    se2_now <- NA_real_
    if(cache1$ok[i] && cache1$ok[j] && cache2$ok[i] && cache2$ok[j]){
      ## variance formula for diff = <z_i,z_j>_1 - <z_i,z_j>_2: by delta method,
      ## var(<z_i,z_j>_1) ≈ z_j1^T Vii1 z_j1 + z_i1^T Vjj1 z_i1
      ## var(<z_i,z_j>_1) ≈ z_j2^T Vii2 z_j2 + z_i2^T Vjj2 z_i2
      ## var(diff) ≈ z_j1^T Vii1 z_j1 + z_i1^T Vjj1 z_i1 + z_j2^T Vii2 z_j2 + z_i2^T Vjj2 z_i2.
      se2_now <- as.numeric(
        t(zj1) %*% Vii1 %*% zj1 +
        t(zi1) %*% Vjj1 %*% zi1 +
        t(zj2) %*% Vii2 %*% zj2 +
        t(zi2) %*% Vjj2 %*% zi2
      )
    }

    if(is.finite(se2_now) && se2_now > 0){
      prod_se[idx] <- sqrt(se2_now)
      prod_z[idx] <- prod_diff[idx] / prod_se[idx]  ## z-stat for pair (i,j) z product 
      prod_p[idx] <- 2 * pnorm(-abs(prod_z[idx]))  ## p-value for pair (i,j) z product
    }

    idx <- idx + 1L
  }
}

innerprod_tests <- data.frame(
  filtered_node_i = pair_i,                # First node index i in filtered network for pair (i, j).
  original_node_i = filtered$keep_idx[pair_i], # Original unfiltered node index of i.
  borough_i = borough[pair_i],             # Borough label of node i.
  filtered_node_j = pair_j,                # Second node index j in filtered network for pair (i, j).
  original_node_j = filtered$keep_idx[pair_j], # Original unfiltered node index of j.
  borough_j = borough[pair_j],             # Borough label of node j.
  hour1 = lab1,                            # First hour label (e.g., "08:00-09:00").
  hour2 = lab2,                            # Second hour label (e.g., "18:00-19:00").
  innerprod_hat_1 = prod_hat_1,            # Estimated <z_i, z_j> at hour 1.
  innerprod_hat_2 = prod_hat_2,            # Estimated <z_i, z_j> at hour 2.
  diff = prod_diff,                        # Difference innerprod_hat_1 - innerprod_hat_2.
  se = prod_se,                            # Standard error of diff for the pair.
  z_stat = prod_z,                         # Z statistic for testing pairwise diff = 0.
  p_value = prod_p,                        # Two-sided p-value from z_stat.
  stringsAsFactors = FALSE
)
innerprod_tests$q_value_bh <- p.adjust(innerprod_tests$p_value, method = "BH")

# -------------------------
# Build node-wise rejection-count summary from pairwise BH results
# -------------------------
rej <- is.finite(innerprod_tests$q_value_bh) & innerprod_tests$q_value_bh < 0.05
rejected_count <- integer(n)
total_tests <- integer(n)
for(r in seq_len(nrow(innerprod_tests))){
  i <- innerprod_tests$filtered_node_i[r]
  j <- innerprod_tests$filtered_node_j[r]
  total_tests[i] <- total_tests[i] + 1L
  total_tests[j] <- total_tests[j] + 1L
  if(rej[r]){
    rejected_count[i] <- rejected_count[i] + 1L
    rejected_count[j] <- rejected_count[j] + 1L
  }
}

rejection_counts <- data.frame(
  filtered_node = seq_len(n),             # Node index in the filtered network (1..n).
  original_node = filtered$keep_idx,      # Corresponding node index in the original unfiltered network.
  borough = borough,                       # Borough label for the node.
  rejected_count = rejected_count,         # Number of pairwise tests involving this node rejected at BH q < 0.05.
  total_tests = total_tests,               # Total number of pairwise tests involving this node.
  rejected_frac = rejected_count / total_tests, # Fraction of rejected pairwise tests for this node.
  stringsAsFactors = FALSE
)

# -------------------------
# Save outputs
# -------------------------
alpha_csv <- paste0(output_prefix, "_alpha.csv")
innerprod_csv <- paste0(output_prefix, "_innerprod.csv")
rejection_counts_csv <- file.path(dirname(output_prefix), "nybike_innerprod_rejections_by_node.csv")

# alpha_csv stores node-wise two-sample tests for alpha_i.
# innerprod_csv stores pairwise two-sample tests for z_i^T z_j.
# rejection_counts_csv stores per-node rejected-pair counts and fractions.
write.csv(alpha_tests, alpha_csv, row.names = FALSE)
write.csv(innerprod_tests, innerprod_csv, row.names = FALSE)
write.csv(rejection_counts, rejection_counts_csv, row.names = FALSE)

cat(sprintf("Saved alpha tests: %s\n", alpha_csv))
cat(sprintf("Saved inner-product tests: %s\n", innerprod_csv))
cat(sprintf("Saved nodewise rejection counts: %s\n", rejection_counts_csv))
