### Working directory requirement:
### Set working directory to BTAP before running this script.
### Example: getwd() should return the BTAP folder path.

### this r code produces Figure S10(a)(b) in the paper.

### Script overview:
### (1) Load fitted alpha estimates for both peak hours and retained-node station locations.
### (2) Transform alpha by exp(alpha) to map intensity scale.
### (3) Build a common size/color scale shared by both hour maps.
### (4) Save outputs:
###     - data analysis/plots/nybike_station_locations_alpha_08.pdf
###       (station map with point size/color from exp(alpha) at 08:00),
###     - data analysis/plots/nybike_station_locations_alpha_18.pdf
###       (station map with point size/color from exp(alpha) at 18:00),
###     using the same scale for comparability across hours.

# -------------------------
# Paths and input files
# -------------------------
results_path <- file.path("data analysis/results", "nybike_fit.rds") # fitted model outputs for 08:00 and 18:00
station_latlon_csv <- file.path("data analysis/results", "nybike_station_locations_latlon.csv") # node station_id/lat/lon table
results <- readRDS(results_path)
metadata <- read.csv(station_latlon_csv, stringsAsFactors = FALSE)
# Keep station rows ordered by filtered-node index for consistent plotting.
metadata <- metadata[order(metadata$filtered_node), ]
## output paths
plot_prefix <- file.path("data analysis/plots", "nybike_station_locations_alpha")


# Convert alpha to baseline exp(alpha)
alpha1 <- exp(as.numeric(results[[1]]$fit$alpha))
alpha2 <- exp(as.numeric(results[[2]]$fit$alpha))
alpha_all <- c(alpha1, alpha2)
alpha_rng <- range(alpha_all, na.rm = TRUE)
alpha_ticks <- c(min(alpha_all), median(alpha_all), max(alpha_all))
alpha_pal <- colorRampPalette(c("#dfe9f3", "#9ecae1", "#3182bd", "#08519c"))(100)

# -------------------------
# Helpers: alpha -> normalized scale / point size / color bin
# -------------------------
# Logic flow:
# (1) alpha_to_unit(alpha_vec, alpha_rng) standardizes alpha to [0,1] on a shared range.
# (2) scale_alpha(...) and alpha_color_index(...) both call alpha_to_unit(...):
#     - scale_alpha maps normalized alpha to point sizes (cex),
#     - alpha_color_index maps normalized alpha to palette bin indices.
# (3) plot_alpha_map(...) calls scale_alpha and alpha_color_index to generate
#     per-node sizes/colors (and legend sizes/colors), then draws the map.

# Map alpha values to [0,1] using the shared global alpha range across both hours.
alpha_to_unit <- function(alpha_vec, alpha_rng){
  if(diff(alpha_rng) < 1e-12){
    return(rep(0.5, length(alpha_vec)))
  }
  ((alpha_vec - alpha_rng[1]) / diff(alpha_rng))^(0.5)
}

# Convert alpha values into plotting point sizes (cex) using the normalized scale.
scale_alpha <- function(alpha_vec, alpha_rng, cex_min = 0.75, cex_max = 2.35){
  alpha_unit <- alpha_to_unit(alpha_vec, alpha_rng)
  if(length(alpha_unit) == 0){
    return(rep((cex_min + cex_max) / 2, length(alpha_vec)))
  }
  cex_min + (cex_max - cex_min) * alpha_unit
}

# Convert alpha values into color-bin indices (1..100) on the shared palette.
alpha_color_index <- function(alpha_vec, alpha_rng){
  alpha_unit <- alpha_to_unit(alpha_vec, alpha_rng)
  if(length(alpha_unit) == 0){
    return(rep(50L, length(alpha_vec)))
  }
  pmax(1L, pmin(100L, ceiling(99 * alpha_unit) + 1L))
}

x_geo_rng <- range(metadata$lon, na.rm = TRUE)
y_geo_rng <- range(metadata$lat, na.rm = TRUE)
x_geo_pad_left <- 0.08 * max(1e-12, diff(x_geo_rng))
x_geo_pad_right <- 0.02 * max(1e-12, diff(x_geo_rng))
y_geo_pad <- 0.08 * max(1e-12, diff(y_geo_rng))

# -------------------------
# function: draw alpha map
# -------------------------
# Draw a station map for a given alpha vector using shared size/color scales
# and a common legend anchored to min/median/max(alpha_all).
plot_alpha_map <- function(alpha_vec, alpha_rng, out_path){
  sizes <- scale_alpha(alpha_vec, alpha_rng)
  color_idx <- alpha_color_index(alpha_vec, alpha_rng)
  cols <- alpha_pal[color_idx]
  legend_sizes <- scale_alpha(alpha_ticks, alpha_rng)
  legend_cols <- alpha_pal[alpha_color_index(alpha_ticks, alpha_rng)]

  pdf(file = out_path, width = 6, height = 5.8)
  old_par <- par(no.readonly = TRUE)
  par(mar = c(4.7, 5.8, 3.2, 1.5))
  plot(
    metadata$lon, metadata$lat,
    col = cols, pch = 16, cex = sizes,
    xlim = c(x_geo_rng[1] - x_geo_pad_left, x_geo_rng[2] + x_geo_pad_right),
    ylim = y_geo_rng + c(-y_geo_pad, y_geo_pad),
    xlab = "Longitude", ylab = "Latitude",
    main = "",
    cex.lab = 1.8, cex.axis = 1.6,
    xaxt = "n"
  )
  tick_x <- pretty(par("usr")[1:2])
  tick_x <- tick_x[tick_x < (par("usr")[2] - 0.01)]
  axis(1, at = tick_x, labels = formatC(tick_x, format = "f", digits = 2), cex.axis = 1.6)
  usr <- par("usr")
  legend(
    x = mean(usr[1:2]),
    y = usr[4] + 0.025 * diff(usr[3:4]),
    legend = c(
      sprintf("Min: %s", formatC(alpha_ticks[1], format = "f", digits = 2)),
      sprintf("Med: %s", formatC(alpha_ticks[2], format = "f", digits = 2)),
      sprintf("Max: %s", formatC(alpha_ticks[3], format = "f", digits = 2))
    ),
    pch = 16,
    col = legend_cols,
    pt.cex = legend_sizes,
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
}

# -------------------------
# Draw maps for 08:00 and 18:00
# -------------------------
plot_alpha_map(alpha1, alpha_rng, sprintf("%s_08.pdf", plot_prefix))
plot_alpha_map(alpha2, alpha_rng, sprintf("%s_18.pdf", plot_prefix))
