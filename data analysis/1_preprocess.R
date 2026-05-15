### Working directory requirement:
### Set working directory to BTAP before running this script.
### Example: getwd() should return the BTAP folder path.

### Script overview:
### (1) Load pre-aggregated NY bike network array data from `NYbike.rda`.
### (2) Extract two target hours (default: 08:00 and 18:00), and for each hour:
###     - keep the directed adjacency,
###     - build a symmetric (undirected) adjacency by A + t(A), with zero diagonal.
### (3) Compute node degrees in both symmetric networks.
### (4) Filter out nodes with low activity:
###     keep only nodes with degree > `degree_threshold` in BOTH hours.
### (5) Save the filtered directed/symmetric networks and metadata to
###     `results/nybike_data.rds` for downstream fitting/visualization scripts.




# Expected network data from NYbike.rda: 
### NYbike1 (3D array: 782 x 782 x 24 dim, i.e., node x node x hour-index).
### borough (782 dim char vector)
load("data analysis/NYbike.rda") 
hour1 <- 8L
hour2 <- 18L
degree_threshold <- 1
out_rds <- file.path("data analysis/results", "nybike_data.rds")

# -------------------------
# Helper: extract one hour network, making it symmetric and no-self-loop
# -------------------------
# NOTE: Hour index in NYbike1 is offset by +1 relative to clock hour. e.g., hour = 8 uses slice [, , 9].
get_symmetric_hour <- function(hour){
  A_dir <- NYbike1[, , hour + 1L]
  A_sym <- A_dir + t(A_dir)
  diag(A_sym) <- 0
  list(A_directed = A_dir, A_symmetric = A_sym)
}

# Build hour-specific directed/symmetric networks.
hour1_net <- get_symmetric_hour(hour1)
hour2_net <- get_symmetric_hour(hour2)

# -------------------------
# Degree-based node filtering
# -------------------------
# Degree is computed from symmetric adjacency so activity is undirected volume.
deg1 <- rowSums(hour1_net$A_symmetric)
deg2 <- rowSums(hour2_net$A_symmetric)

# Keep nodes active in BOTH hours under the threshold rule.
keep <- (deg1 > degree_threshold) & (deg2 > degree_threshold)
### 703 out of 782 nodes were kept

# -------------------------
# Package filtered output
# -------------------------
# This list is the standardized preprocessing output consumed by later scripts.
filtered <- list(
  hour1 = hour1,  # Integer: first selected hour (clock hour) used in analysis.
  hour2 = hour2,  # Integer: second selected hour (clock hour) used in analysis.
  keep_idx = which(keep),  # Integer vector: original node indices retained after filtering.
  borough = borough[keep],  # Character vector: borough labels for retained nodes only (aligned with filtered matrices).
  A_directed_hour1 = hour1_net$A_directed[keep, keep, drop = FALSE],  # Matrix (kept x kept): directed adjacency at hour1.
  A_symmetric_hour1 = hour1_net$A_symmetric[keep, keep, drop = FALSE],  # Matrix (kept x kept): symmetric adjacency at hour1.
  A_directed_hour2 = hour2_net$A_directed[keep, keep, drop = FALSE],  # Matrix (kept x kept): directed adjacency at hour2.
  A_symmetric_hour2 = hour2_net$A_symmetric[keep, keep, drop = FALSE]  # Matrix (kept x kept): symmetric adjacency at hour2.
)

# -------------------------
# Save and report
# -------------------------
saveRDS(filtered, out_rds)

cat(sprintf("Saved filtered data: %s\n", out_rds))
cat(sprintf("Kept %d of %d nodes after dropping nodes with total degree <= %.0f in either hour.\n",
            sum(keep), length(keep), degree_threshold))
