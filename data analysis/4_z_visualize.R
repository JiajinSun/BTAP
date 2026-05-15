### Working directory requirement:
### Set working directory to BTAP before running this script.
### Example: getwd() should return the BTAP folder path.

### this r code produces Figure 5(a) and Figure S9 in the paper.
### Script overview:
### (1) Load preprocessed data and fitted NY bike model results.
### (2) Rebuild station geolocation (lat and lon) from raw trip records (Aug 1 subset), and match with the order in the preprocessed data and algorithm fit.
### (3) Rotate latent embeddings from each hour toward true geographic (lat and lon) orientation.
### (4) Save station location table (nybike_station_locations_latlon.csv) to be used by other r scripts (the ones starting with 6_,7_ and 8_).
### (5) Plot: latent map for hour 08, latent map for hour 18, and geographic station map.

# -------------------------
# Paths and files
# -------------------------
filtered_data_path <- file.path("data analysis/results", "nybike_data.rds") # filtered data
results_path <- file.path("data analysis/results", "nybike_fit.rds") # fit of proposed algorithm
csv_path <- "data analysis/201908-citibike-tripdata_aug01.csv" # for extracting the station location info
plot_prefix <- file.path("data analysis/plots", "nybike")
station_latlon_csv <- file.path("data analysis/results", "nybike_station_locations_latlon.csv")
results <- readRDS(results_path)
filtered <- readRDS(filtered_data_path)

# -------------------------
# Build station metadata from raw trip CSV
# -------------------------
# Output of this block (used later):
# - station_df (data.frame; one row per original node):
#     * original_node: integer node index in the NYbike1 ordering (1..N).
#     * station_id   : Citi Bike station ID mapped to that node.
#     * lat, lon     : station latitude/longitude.
# - geo_df (data.frame; one row per retained node):
#     * same columns as station_df, but rows are re-ordered/subset to match
#       filtered$keep_idx exactly.
# - geo_xy (numeric matrix; n_kept x 2):
#     * column 1 = lon, column 2 = lat, in the same row order as geo_df.
#     * used by latent-to-geo rotation and by the geographic location plot.
trip_cols <- c(
  "starttime", "tripduration",
  "start station id", "start station name", "start station latitude", "start station longitude",
  "end station id", "end station name", "end station latitude", "end station longitude"
)
trip <- read.csv(csv_path, check.names = FALSE)[, trip_cols]

trip_del <- trip[
  trip[["start station id"]] != trip[["end station id"]] &
    trip$tripduration > 60 &
    trip$tripduration < 3 * 60 * 60,
]

stations <- unique(sort(as.numeric(c(trip_del[["start station id"]], trip_del[["end station id"]]))))

start_meta <- unique(trip[, c("start station id", "start station latitude", "start station longitude")])
names(start_meta) <- c("station_id", "lat", "lon")
end_meta <- unique(trip[, c("end station id", "end station latitude", "end station longitude")])
names(end_meta) <- c("station_id", "lat", "lon")
station_meta <- unique(rbind(start_meta, end_meta))
station_meta <- station_meta[!duplicated(station_meta$station_id), ]
station_meta <- station_meta[station_meta$station_id != "NULL", ]
station_meta$station_id <- as.numeric(station_meta$station_id)

# Map original-node indexing to station coordinates.
station_df <- data.frame(
  original_node = seq_along(stations),
  station_id = stations,
  stringsAsFactors = FALSE
)
station_df <- merge(station_df, station_meta, by = "station_id", all.x = TRUE, sort = FALSE)
station_df <- station_df[order(station_df$original_node), ]
geo_df <- station_df[match(filtered$keep_idx, station_df$original_node), ]
geo_xy <- as.matrix(geo_df[, c("lon", "lat")])
station_latlon_df <- data.frame(
  filtered_node = seq_len(nrow(geo_df)),
  station_id = geo_df$station_id,
  lat = geo_df$lat,
  lon = geo_df$lon,
  stringsAsFactors = FALSE
)

write.csv(station_latlon_df, station_latlon_csv, row.names = FALSE)
### (4) Save station location table (nybike_station_locations_latlon.csv) to be used by other r scripts.


# -------------------------
# Extract fitted objects and plotting labels
# -------------------------
fit_08 <- results[[1]][["fit"]]
fit_18 <- results[[2]][["fit"]]
borough <- filtered$borough
# Borough color mapping for points.
borough_factor <- factor(borough)
borough_levels <- levels(borough_factor)
borough_cols <- c("Brooklyn" = "#1b9e77", "Manhattan" = "#d95f02", "Queens" = "#7570b3")
col_map <- borough_cols[borough_levels]
station_cols <- col_map[as.integer(borough_factor)]

true_location_plot_path <- paste0(plot_prefix, "_true_locations.pdf")

# -------------------------
# Helper: Procrustes-style rotation from latent to geographic frame
# -------------------------
compute_geo_rotation <- function(z, geo_xy){
  z_centered <- scale(z, center = TRUE, scale = FALSE)
  geo_centered <- scale(geo_xy, center = TRUE, scale = FALSE)
  sv <- svd(t(z_centered) %*% geo_centered)
  sv$u %*% t(sv$v)
}

z_08_rot <- scale(fit_08$Z, center = TRUE, scale = FALSE) %*% compute_geo_rotation(fit_08$Z, geo_xy)
z_18_rot <- scale(fit_18$Z, center = TRUE, scale = FALSE) %*% compute_geo_rotation(fit_18$Z, geo_xy)

# -------------------------
# Helper: latent-space scatter by hour
# -------------------------
plot_hour_latent <- function(z, out_file, hour_index){
  x_lab <- "1st component"
  y_lab <- "2nd component"
  x_rng <- range(z[, 1], na.rm = TRUE)
  y_rng <- range(z[, 2], na.rm = TRUE)
  x_pad <- 0.04 * max(1e-12, diff(x_rng))
  y_pad <- 0.06 * max(1e-12, diff(y_rng))

  pdf(file = out_file, width = 6, height = 5.8)
  old_par <- par(no.readonly = TRUE)
  on.exit({ par(old_par); dev.off() }, add = TRUE)
  par(mar = c(4.7, 5.8, 3.2, 1.5))

  plot(z[, 1], z[, 2], col = station_cols, pch = 16, cex = 0.85,
       xlim = x_rng + c(-x_pad, x_pad),
       ylim = y_rng + c(-y_pad, y_pad),
       xlab = "", ylab = y_lab,
       main = "",
       cex.lab = 1.8, cex.axis = 1.6)
  mtext(x_lab, side = 1, line = 3.4, cex = 1.8)
  usr <- par("usr")
  legend(x = mean(usr[1:2]), y = usr[4] + 0.025 * diff(usr[3:4]),
         legend = borough_levels, col = col_map, pch = 16,
         pt.cex = 1.5, cex = 1.6, horiz = TRUE, xpd = NA, bty = "n",
         xjust = 0.5, yjust = 0)
}

# -------------------------
# Plot latent-space maps for hour 08 and 18
# -------------------------
for(result in results){
  if(result$hour == 8){
    z_plot <- z_08_rot
  } else if(result$hour == 18){
    z_plot <- z_18_rot
  }
  hour_index <- if(result$hour == 8) 1 else 2
  plot_hour_latent(z_plot, sprintf("%s_hour_%02d_latent.pdf", plot_prefix, result$hour), hour_index)
}

# -------------------------
# Plot true station geographic locations
# -------------------------
x_geo_rng <- range(geo_xy[, 1], na.rm = TRUE)
y_geo_rng <- range(geo_xy[, 2], na.rm = TRUE)
x_geo_pad_left <- 0.08 * max(1e-12, diff(x_geo_rng))
x_geo_pad_right <- 0.02 * max(1e-12, diff(x_geo_rng))
y_geo_pad <- 0.08 * max(1e-12, diff(y_geo_rng))

pdf(file = true_location_plot_path, width = 6, height = 5.8)
old_par_geo <- par(no.readonly = TRUE)
par(mar = c(4.7, 5.8, 3.2, 1.5))
plot(
  geo_xy[, 1], geo_xy[, 2],
  col = station_cols,
  pch = 16,
  cex = 0.85,
  xlim = c(x_geo_rng[1] - x_geo_pad_left, x_geo_rng[2] + x_geo_pad_right),
  ylim = y_geo_rng + c(-y_geo_pad, y_geo_pad),
  xlab = "Longitude",
  ylab = "Latitude",
  main = "",
  cex.lab = 1.8, cex.axis = 1.6,
  xaxt = "n"
)
tick_x <- pretty(par("usr")[1:2])
tick_x <- tick_x[tick_x < (par("usr")[2] - 0.01)]
axis(1, at = tick_x, labels = formatC(tick_x, format = "f", digits = 2), cex.axis = 1.6)
usr <- par("usr")
legend(x = mean(usr[1:2]), y = usr[4] + 0.025 * diff(usr[3:4]),
       legend = borough_levels, col = col_map, pch = 16,
       pt.cex = 1.5, cex = 1.6, horiz = TRUE, xpd = NA, bty = "n",
       xjust = 0.5, yjust = 0)
par(old_par_geo)
dev.off()


