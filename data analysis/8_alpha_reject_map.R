### Working directory requirement:
### Set working directory to BTAP before running this script.
### Example: getwd() should return the BTAP folder path.

### this r code produces Figure S10(c) in the paper.

### Script overview:
### (1) Load alpha two-sample test results and retained-node station locations.
### (2) Keep BH-significant nodes (q < bh_level) and compute p-value intensity score.
### (3) Map stronger evidence (smaller p) to larger/darker points.
### (4) Save output:
###     - data analysis/plots/nybike_station_locations_alpha_pvalue.pdf
###       (all stations in gray, BH-significant alpha-difference nodes highlighted
###       by point size/color according to p-value strength).

# -------------------------
# Paths and input files
# -------------------------
alpha_tests_path <- file.path("data analysis/results", "nybike_08_vs_18_alpha.csv") # alpha two-sample test table
station_latlon_csv <- file.path("data analysis/results", "nybike_station_locations_latlon.csv") # node station_id/lat/lon table
alpha_tests <- read.csv(alpha_tests_path, stringsAsFactors = FALSE)
metadata <- read.csv(station_latlon_csv, stringsAsFactors = FALSE)
# Keep station rows ordered by filtered-node index before matching.
metadata <- metadata[order(metadata$filtered_node), ]
# Align station metadata rows to alpha test row order by filtered node id.
meta_idx <- match(alpha_tests$filtered_node, metadata$filtered_node)
plot_df <- cbind(
  alpha_tests,
  metadata[meta_idx, c("lat", "lon")]
)
# output paths
plot_path <- file.path("data analysis/plots", "nybike_station_locations_alpha_pvalue.pdf")

# BH significance level used to highlight alpha-difference nodes.
bh_level <- 0.2

# -------------------------
# Derive significant subset and visual scale from p-values
# -------------------------
# Logic flow:
# (1) p_to_unit(p_vec, score_rng, p_eps) maps transformed p-value score to [0,1].
# (2) scale_p(...) and p_color_index(...) both call p_to_unit(...):
#     - scale_p maps normalized transformed p-value scores to point sizes (cex),
#     - p_color_index maps normalized transformed p-value scores to palette bin indices.
# (3) The resulting sizes/colors are used to draw highlighted significant nodes
#     and their legend entries in the final station map.
pvals <- as.numeric(plot_df$p_value)
qvals <- as.numeric(plot_df$q_value_bh)
sig_idx <- which(is.finite(pvals) & is.finite(qvals) & qvals < bh_level)

pvals_sig <- pvals[sig_idx]
p_eps <- 1e-12
score <- -log10(pmax(pvals_sig, p_eps))
score_rng <- range(score, na.rm = TRUE)
p_ticks <- c(min(pvals_sig, na.rm = TRUE), median(pvals_sig, na.rm = TRUE), max(pvals_sig, na.rm = TRUE))

# Map p-values to [0,1] by -log10(p), so smaller p-values get larger visual weight.
p_to_unit <- function(p_vec, score_rng, p_eps){
  score_vec <- -log10(pmax(p_vec, p_eps))
  if(diff(score_rng) < 1e-12){
    return(rep(0.5, length(score_vec)))
  }
  ((score_vec - score_rng[1]) / diff(score_rng))^(0.5)
}

# Convert p-values into plotting point sizes (cex) using the normalized p-score.
scale_p <- function(p_vec, score_rng, p_eps, cex_min = 0.75, cex_max = 2.35){
  p_unit <- p_to_unit(p_vec, score_rng, p_eps)
  cex_min + (cex_max - cex_min) * p_unit
}

# Convert p-values into color-bin indices (1..100) for the purple palette.
p_color_index <- function(p_vec, score_rng, p_eps){
  p_unit <- p_to_unit(p_vec, score_rng, p_eps)
  pmax(1L, pmin(100L, ceiling(99 * p_unit) + 1L))
}

pal <- colorRampPalette(c("#f2f0f7", "#cbc9e2", "#9e9ac8", "#6a51a3"))(100)
cols <- pal[p_color_index(pvals_sig, score_rng, p_eps)]
sizes <- scale_p(pvals_sig, score_rng, p_eps)
legend_cols <- pal[p_color_index(p_ticks, score_rng, p_eps)]
legend_sizes <- scale_p(p_ticks, score_rng, p_eps)

x_geo_rng <- range(plot_df$lon, na.rm = TRUE)
y_geo_rng <- range(plot_df$lat, na.rm = TRUE)
x_geo_pad_left <- 0.08 * max(1e-12, diff(x_geo_rng))
x_geo_pad_right <- 0.02 * max(1e-12, diff(x_geo_rng))
y_geo_pad <- 0.08 * max(1e-12, diff(y_geo_rng))

# -------------------------
# Draw p-value significance map
# -------------------------
pdf(file = plot_path, width = 6, height = 5.8)
old_par <- par(no.readonly = TRUE)
par(mar = c(4.7, 5.8, 3.2, 1.5))
plot(
  plot_df$lon, plot_df$lat,
  col = "#d9d9d9", pch = 16, cex = 0.45,
  xlim = c(x_geo_rng[1] - x_geo_pad_left, x_geo_rng[2] + x_geo_pad_right),
  ylim = y_geo_rng + c(-y_geo_pad, y_geo_pad),
  xlab = "Longitude", ylab = "Latitude",
  main = "",
  cex.lab = 1.8, cex.axis = 1.6,
  xaxt = "n"
)
points(plot_df$lon[sig_idx], plot_df$lat[sig_idx], col = cols, pch = 16, cex = sizes)
tick_x <- pretty(par("usr")[1:2])
tick_x <- tick_x[tick_x < (par("usr")[2] - 0.01)]
axis(1, at = tick_x, labels = formatC(tick_x, format = "f", digits = 2), cex.axis = 1.6)
usr <- par("usr")
legend(
  x = mean(usr[1:2]),
  y = usr[4] + 0.025 * diff(usr[3:4]),
  legend = c(
    formatC(p_ticks[1], format = "f", digits = 2),
    formatC(p_ticks[2], format = "f", digits = 2)
  ),
  pch = 16,
  col = legend_cols[1:2],
  pt.cex = legend_sizes[1:2],
  bty = "n",
  cex = 1.7,
  horiz = TRUE,
  x.intersp = 0.6,
  xpd = NA,
  xjust = 0.5,
  yjust = 0.1
)
par(old_par)
dev.off()
