### Working directory requirement:
### Set working directory to BTAP before running this script.
### Example: getwd() should return the BTAP folder path.

### this r code produces Figure 6(c) in the paper.

### Script overview:
### (1) Load preprocessing output and node-level inner-product rejection summaries.
### (2) Load retained-node station locations (filtered_node, station_id, lat, lon).
### (3) Align filtered nodes to station coordinates.
### (4) Build a map where color/size encode rejected-pair fraction per node.
### (5) Save outputs:
###     - data analysis/plots/nybike_station_locations_innerprod_fraction.pdf
###       (station map with color/size encoding node-level rejected fraction).

# -------------------------
# Paths and files
# -------------------------
filtered_data_path <- file.path("data analysis/results", "nybike_data.rds") # filtered node set from preprocessing
rejection_counts_path <- file.path("data analysis/results", "nybike_innerprod_rejections_by_node.csv") # node-level rejection summary
station_latlon_csv <- file.path("data analysis/results", "nybike_station_locations_latlon.csv") # node station_id/lat/lon 
filtered <- readRDS(filtered_data_path)
rejection_counts <- read.csv(rejection_counts_path, stringsAsFactors = FALSE)
station_latlon <- read.csv(station_latlon_csv, stringsAsFactors = FALSE)
# Keep station rows ordered by filtered-node index for consistent merging/alignment.
station_latlon <- station_latlon[order(station_latlon$filtered_node), ]
# output path
plot_prefix <- file.path("data analysis/plots", "nybike_station_locations")

# -------------------------
# Merge node-level rejection metrics with station metadata
# -------------------------
map_df <- merge(rejection_counts, station_latlon, by = "filtered_node", all.x = TRUE, sort = FALSE)
map_df$is_innerprod_rejected <- map_df$rejected_count > 0

x_geo_rng <- range(map_df$lon, na.rm = TRUE)
y_geo_rng <- range(map_df$lat, na.rm = TRUE)
x_geo_pad_left <- 0.08 * max(1e-12, diff(x_geo_rng))
x_geo_pad_right <- 0.02 * max(1e-12, diff(x_geo_rng))
y_geo_pad <- 0.08 * max(1e-12, diff(y_geo_rng))

fraction_map_path <- paste0(plot_prefix, "_innerprod_fraction.pdf")

# -------------------------
# function to map rejected fraction to color + point size
# -------------------------
plot_fraction_map <- function(df, out_path){
  pal <- colorRampPalette(c("#fff5eb", "#fdcc8a", "#fc8d59", "#d7301f"))(100)
  frac <- df$rejected_frac
  idx <- pmax(1L, pmin(100L, ceiling(100 * (frac - min(frac)) / max(1e-12, diff(range(frac))))))
  cols <- pal[idx]
  sizes <- 0.6 + 1.8 * frac

  pdf(file = out_path, width = 6, height = 5.8)
  old_par <- par(no.readonly = TRUE)
  par(mar = c(4.7, 5.8, 3.2, 1.5))
  plot(
    df$lon, df$lat,
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
    x = mean(usr[1:2]), y = usr[4] + 0.025 * diff(usr[3:4]),
    legend = c(
      sprintf("Min: %.2f", min(frac)),
      sprintf("Med: %.2f", median(frac)),
      sprintf("Max: %.2f", max(frac))
    ),
    pch = 16,
    col = c(pal[1], pal[50], pal[100]),
    pt.cex = c(0.6 + 1.8 * min(frac), 0.6 + 1.8 * median(frac), 0.6 + 1.8 * max(frac)),
    bty = "n", cex = 1.7, horiz = TRUE, x.intersp = 0.6, xpd = NA, xjust = 0.5, yjust = 0.1
  )
  par(old_par)
  dev.off()
}

# -------------------------
# Draw map and report outputs
# -------------------------
plot_fraction_map(map_df, fraction_map_path)
