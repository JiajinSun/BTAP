### Working directory requirement:
### Set working directory to BTAP before running this script.
### Example: getwd() should return the BTAP folder path.

### this r code produces Figure 5(b)(c) and Figure 6(a)(b) in the paper.

### Script overview:
### (1) Load fitted latent positions for 08:00 and 18:00 plus pairwise test results.
### (2) Build Gram matrices G1 = Z1 Z1^T, G2 = Z2 Z2^T, and difference D = G1 - G2.
### (3) Reorder nodes by borough and latent coordinates for block-structured visualization.
### (4) Create heatmaps for hour-1 inner products, hour-2 inner products, and differences.
### (5) Create a rejected-only heatmap (entries with BH q < 0.05).
### (6) Save outputs:
###     - data analysis/plots/nybike_innerprod_hour1.pdf            (heatmap of Z_08 Z_08^T),
###     - data analysis/plots/nybike_innerprod_hour2.pdf            (heatmap of Z_18 Z_18^T),
###     - data analysis/plots/nybike_innerprod_difference.pdf       (heatmap of Z_08 Z_08^T - Z_18 Z_18^T),
###     - data analysis/plots/nybike_innerprod_rejected_only.pdf    (difference heatmap masked by BH-significant pairs),
###     - data analysis/results/nybike_innerprod_heatmap_order.csv  (node ordering used by all heatmaps).

# -------------------------
# paths and files
# -------------------------
results_path <- file.path("data analysis/results", "nybike_fit.rds")  # filtered data
filtered_data_path <- file.path("data analysis/results", "nybike_data.rds")  # fit of proposed algorithm
innerprod_tests_path <- file.path("data analysis/results", "nybike_08_vs_18_innerprod.csv")  # innerprod test results
results <- readRDS(results_path)
filtered <- readRDS(filtered_data_path)
innerprod_tests <- read.csv(innerprod_tests_path, stringsAsFactors = FALSE)
## output paths
plot_prefix <- file.path("data analysis/plots", "nybike")
order_csv_path <- file.path("data analysis/results", "nybike_innerprod_heatmap_order.csv")
hour1_path <- paste0(plot_prefix, "_innerprod_hour1.pdf")
hour2_path <- paste0(plot_prefix, "_innerprod_hour2.pdf")
difference_path <- paste0(plot_prefix, "_innerprod_difference.pdf")
rejected_only_path <- paste0(plot_prefix, "_innerprod_rejected_only.pdf")


fit1 <- results[[1]]$fit
fit2 <- results[[2]]$fit
Z1 <- fit1$Z
Z2 <- fit2$Z
n <- nrow(Z1)  ## n=703
borough <- filtered$borough  ## length(borough) = 703

# -------------------------
# Node ordering for block visualization
# -------------------------
order_idx <- order(borough, Z1[, 1], Z1[, 2], decreasing = FALSE)
inv_order <- integer(n)
inv_order[order_idx] <- seq_len(n)

# -------------------------
# Inner-product matrices and difference
# -------------------------
G1 <- Z1 %*% t(Z1)
G2 <- Z2 %*% t(Z2)
D <- G1 - G2

# -------------------------
# Build symmetric rejection mask from BH-significant pairwise tests
# -------------------------
rej_mat <- matrix(FALSE, n, n)
rej_rows <- which(is.finite(innerprod_tests$q_value_bh) & innerprod_tests$q_value_bh < 0.05)
ii <- innerprod_tests$filtered_node_i[rej_rows]
jj <- innerprod_tests$filtered_node_j[rej_rows]
rej_mat[cbind(ii, jj)] <- TRUE
rej_mat[cbind(jj, ii)] <- TRUE
diag(rej_mat) <- FALSE

# Reordered matrices for plotting.
G1o <- G1[order_idx, order_idx]
G2o <- G2[order_idx, order_idx]
Do <- D[order_idx, order_idx]
rejo <- rej_mat[order_idx, order_idx]

borough_ord <- borough[order_idx]
block_end <- cumsum(as.integer(table(borough_ord)))
block_start <- c(1, head(block_end + 1L, -1L))
block_mid <- (block_start + block_end) / 2
borough_levels <- names(table(borough_ord))



# -------------------------
# Plotting helpers
# -------------------------
lim_G <- max(abs(c(G1o, G2o)))
pal_G <- colorRampPalette(c("#2166ac", "#92c5de", "#e6e6e6", "#f4a582", "#ca0020"))(256)

# Draw one heatmap panel from matrix M with borough block boundaries.
draw_heat <- function(M, pal, zlim){
  image(
    x = seq_len(n), y = seq_len(n),
    z = t(M[n:1, ]),
    col = pal, zlim = c(-zlim, zlim),
    xaxt = "n", yaxt = "n", xlab = "", ylab = "", main = "", useRaster = TRUE
  )
  abline(v = block_end + 0.5, h = n - block_end + 0.5, col = adjustcolor("#000000", 0.18), lwd = 0.7)
}

# Draw the horizontal color key used by heatmaps, with endpoint labels for the z-scale.
draw_color_key <- function(pal){
  ncol_key <- length(pal)
  plot(c(0, ncol_key), c(0, 1), type = "n", axes = FALSE, xlab = "", ylab = "",
       main = "", xaxs = "i", yaxs = "i")
  x_lo <- 0
  x_hi <- ncol_key
  y_lo <- 0.46
  y_hi <- 0.72
  x_seq <- seq(x_lo, x_hi, length.out = ncol_key + 1L)
  for(ii in seq_len(ncol_key)){
    rect(x_seq[ii], y_lo, x_seq[ii + 1L], y_hi, col = pal[ii], border = pal[ii])
  }
  rect(x_lo, y_lo, x_hi, y_hi, border = "black")
  tick_x <- c(x_lo, x_hi)
  tick_y0 <- y_lo
  tick_y1 <- y_lo - 0.10
  segments(tick_x, tick_y0, tick_x, tick_y1, lwd = 1)
  text(tick_x, tick_y1 - 0.22, labels = c("-22", "22"), cex = 1.4, xpd = TRUE)
}

# Add borough labels on both axes at block centers defined by block_mid.
draw_axes <- function(){
  axis(1, at = block_mid, labels = FALSE)
  axis(2, at = n - block_mid + 1, labels = borough_levels, las = 2, cex.axis = 1.5)
  text(
    x = block_mid,
    y = par("usr")[3] - 0.06 * diff(par("usr")[3:4]),
    labels = borough_levels,
    srt = 45,
    adj = 1,
    xpd = TRUE,
    cex = 1.5
  )
}

# Create one PDF figure with a top color key and a bottom heatmap panel.
# This wrapper standardizes layout/margins so all output panels share the same style.
draw_single_panel <- function(out_path, M, pal, zlim){
  pdf(file = out_path, width = 6, height = 5.8)
  old_par <- par(no.readonly = TRUE)
  layout(matrix(c(1, 2), nrow = 2), heights = c(0.15, 1))
  par(mar = c(1.0, 7.3, 0.1, 1.0))
  draw_color_key(pal)
  par(mar = c(6.0, 7.3, 0.3, 1.0))
  draw_heat(M, pal, zlim)
  draw_axes()
  par(old_par)
  dev.off()
}

# -------------------------
# Draw three full heatmaps: hour1, hour2, difference
# -------------------------
draw_single_panel(
  hour1_path,
  G1o,
  pal_G,
  lim_G
)

draw_single_panel(
  hour2_path,
  G2o,
  pal_G,
  lim_G
)

draw_single_panel(
  difference_path,
  Do,
  pal_G,
  lim_G
)

# -------------------------
# Draw rejected-only difference heatmap
# -------------------------
rej_only <- matrix(NA_real_, n, n)
rej_only[rejo] <- Do[rejo]

pdf(file = rejected_only_path, width = 6, height = 5.8)
old_par2 <- par(no.readonly = TRUE)
layout(matrix(c(1, 2), nrow = 2), heights = c(0.15, 1))
par(mar = c(1.0, 7.3, 0.1, 1.0))
draw_color_key(pal_G)
par(mar = c(6.0, 7.3, 0.3, 1.0))
image(
  x = seq_len(n), y = seq_len(n),
  z = t(rej_only[n:1, ]),
  col = pal_G, zlim = c(-lim_G, lim_G),
  xaxt = "n", yaxt = "n", xlab = "", ylab = "",
  main = "",
  useRaster = TRUE
)
abline(v = block_end + 0.5, h = n - block_end + 0.5, col = adjustcolor("#000000", 0.18), lwd = 0.7)
draw_axes()
par(old_par2)
dev.off()

# -------------------------
# Save ordering table used by heatmaps
# -------------------------
order_df <- data.frame(
  heatmap_rank = seq_len(n),
  filtered_node = order_idx,
  original_node = filtered$keep_idx[order_idx],
  borough = borough_ord,
  z1_hour1 = Z1[order_idx, 1],
  z2_hour1 = Z1[order_idx, 2],
  stringsAsFactors = FALSE
)
write.csv(order_df, order_csv_path, row.names = FALSE)

cat("Wrote:", hour1_path, "\n")
cat("Wrote:", hour2_path, "\n")
cat("Wrote:", difference_path, "\n")
cat("Wrote:", order_csv_path, "\n")
cat("Wrote:", rejected_only_path, "\n")
