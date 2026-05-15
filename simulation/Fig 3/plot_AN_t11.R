### Script overview:
### (1) Locate eta0-only asymptotic-record RDS files in
###     `../results for asymptotic distribution`.
### (2) For each replication record, extract `t11_hat` from the tail of `err_vec`,
###     and assemble one long table with identifiers
###     (`family`, `n`, `file_base`, `rep`, `seed`).
### (3) For each `(family, n)`, keep finite `t11_hat` values and create:
###     - QQ plot versus N(0,1)
###     - Histogram with N(0,1) density overlay
### (4) Output files (saved in `out_dir`):
###     - `qq_t_z11_{family}_n{n}.pdf`
###     - `hist_t_z11_{family}_n{n}.pdf`

## some global plotting parameters
textsize <- 27

library(ggplot2)

# directories
### Working directory requirement:
### Set working directory to BTAP before running this script.
### Example: getwd() should return the BTAP folder path.
current_dir <- normalizePath("simulation/Fig 3")
# Input folder
data_dir <- normalizePath("simulation/results for asymptotic distribution")
# Output folder for generated figures.
out_dir <- current_dir

# Collect all eta0-only asymptotic-record files.
rds_files <- sort(Sys.glob(file.path(data_dir, "eta0_only_*_n*_asymptotic_records.rds")))
if(length(rds_files) == 0){
  stop("No eta0-only asymptotic record files found in: ", data_dir)
}


### Data preprocessing summary (before plotting):
### (1) Parse each RDS filename into metadata `meta_list`, where each element is
###     list(family, n, file, file_base):
###       - file     : full file path, used by readRDS(file),
###       - file_base: basename(file), saved for source-tracking in output rows.
### (2) Loop through all replications in all RDS files and append one row per replication
###     into `all_rows`:
###     family, n, file_base, rep, seed, t11_hat.
### (3) Bind the list into one long data frame `df` via do.call(rbind, all_rows).

# Function: parse_meta
# Purpose:
#   Parse family and n from one RDS filename.
# Input:
#   path: one RDS filepath.
# Output:
#   list(family, n, file, file_base) or NULL if pattern mismatch.
parse_meta <- function(path){
  bn <- basename(path)
  # Example:
  # eta0_only_poisson_n500_asymptotic_records.rds
  m <- regexec("^eta0_only_(bernoulli|poisson|gaussian)_n([0-9]+)_asymptotic_records\\.rds$", bn)
  g <- regmatches(bn, m)[[1]]
  if(length(g) == 0){
    return(NULL)
  }
  list(family = g[2], n = as.integer(g[3]), file = path, file_base = bn)
}

meta_list <- lapply(rds_files, parse_meta)
meta_list <- Filter(Negate(is.null), meta_list)

# `meta_list` is a list with one element per matched RDS file; each element is
# list(family, n, file, file_base), used to drive file loading and labeling.

# Build one combined long table over all files and all replications.
all_rows <- list()
for(meta in meta_list){
  recs <- readRDS(meta$file)
  if(length(recs) == 0L) next
  for(i in seq_along(recs)){
    rec <- recs[[i]]
    ev <- rec$err_vec
    all_rows[[length(all_rows) + 1L]] <- data.frame(
      family = meta$family,
      n = meta$n,
      file_base = meta$file_base,
      rep = if(!is.null(rec$rep)) as.integer(rec$rep) else i,
      seed = if(!is.null(rec$seed)) as.integer(rec$seed) else NA_integer_,
      t11_hat = ev[length(ev) - 1L],
      stringsAsFactors = FALSE
    )
  }
}

# One long table: one row = one replication.
df <- do.call(rbind, all_rows)


#### plotting: QQ plot and histogram for t11_hat
# Function: make_qq
# Purpose:
#   Construct QQ plot against N(0,1) for one numeric vector.
make_qq <- function(x, axis_label_expr){
  ggplot(data.frame(x = x), aes(sample = x)) +
    stat_qq(distribution = qnorm) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    coord_equal(xlim = c(-3.5, 3.5), ylim = c(-3.5, 3.5), expand = FALSE) +
    # labs(title = NULL, x = "Theoretical Quantiles (N(0,1))", y = axis_label_expr) +
    labs(title = NULL, x = "Theoretical Quantiles", y = axis_label_expr) +
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
make_hist <- function(x, axis_label_expr){
  ggplot(data.frame(x = x), aes(x = x)) +
    geom_histogram(aes(y = after_stat(density)), bins = 30, color = "white") +
    stat_function(fun = dnorm, args = list(mean = 0, sd = 1), linewidth = 1) +
    coord_cartesian(xlim = c(-3.5, 3.5)) +
    labs(title = NULL, x = axis_label_expr, y = "Density") +
    theme_minimal(base_size = textsize) +
    theme(
      plot.title = element_text(size = textsize, face = "bold", hjust = 0.5),
      axis.title = element_text(size = textsize),
      axis.text = element_text(size = textsize),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "gray92", linewidth = 0.25),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7),
      plot.margin = margin(2, 2, -0.5, 2)
    )
}

# Function: plot_family
# Purpose:
#   For one family, generate QQ/hist PDFs across all available n values.
plot_family <- function(family_name){
  sub <- df[df$family == family_name, , drop = FALSE]
  ns <- sort(unique(sub$n))
  if(length(ns) == 0) return(NULL)

  for(i in seq_along(ns)){
    n <- ns[i]    ### n
    s <- sub[sub$n == n, , drop = FALSE]
    s_hat <- s[is.finite(s$t11_hat), , drop = FALSE]
    xh <- s_hat$t11_hat
    seed_vec <- stats::na.omit(s_hat$seed)
    n_seed <- length(unique(seed_vec))
    cat(sprintf("[t_z11] family=%s n=%d: n_rep=%d, unique_seed=%d\n",
                family_name, n, length(xh), n_seed))

    # y_qq_hat <- expression(paste("Sample Quantiles (", hat(t)(z[11]), ")"))
    y_qq_hat <- "Sample Quantiles"
    x_hist_hat <- expression(hat(t)(z[11]))
    qq_out <- file.path(out_dir, sprintf("qq_t_z11_%s_n%d.pdf", family_name, n))
    hist_out <- file.path(out_dir, sprintf("hist_t_z11_%s_n%d.pdf", family_name, n))
    ggsave(
      qq_out, make_qq(xh, y_qq_hat),
      width = 5.6, height = 5.5, units = "in", device = "pdf", useDingbats = FALSE
    )
    ggsave(
      hist_out, make_hist(xh, x_hist_hat),
      width = 5.6, height = 5.6, units = "in", device = "pdf", useDingbats = FALSE
    )
  }
}

# Generate all plots in a fixed family order for reproducible file ordering.
for(fam in c("bernoulli", "poisson", "gaussian")){
  plot_family(fam)
}

cat("Saved merged direct-only t11 outputs to:\n")
cat(out_dir, "\n")
