### Script overview:
### (1) Read n=1000 run-level results from:
###     ../results for GD behavior
### (2) Keep line-search runs only (method == "linesearch").
### (3) For each replication, use bt_mean = average number of backtracks per GD step.
### (4) Plot boxplots of bt_mean vs eta_init / eta_0.
### (5) Create one PDF per family (bernoulli/poisson/gaussian).
### (6) Save:
###     boxplot_avg_backtracks_n1000_{family}.pdf

library(ggplot2)
# directories
### Working directory requirement:
### Set working directory to BTAP before running this script.
### Example: getwd() should return the BTAP folder path.
current_dir <- normalizePath("simulation/Fig 2 boxplot")
# Input folder
data_dir <- normalizePath("simulation/results for GD behavior")
# Output folder for generated figures.
out_dir <- current_dir

# One run CSV per family at n=1000.
run_files <- sort(Sys.glob(file.path(data_dir, "fs_vs_ls_4etas_*_n1000_run.csv")))

all_rows <- lapply(run_files, function(f){
  read.csv(f, stringsAsFactors = FALSE)
})
df <- do.call(rbind, all_rows)

# Keep line-search rows
ls_df <- df[df$method == "linesearch", , drop = FALSE]

# Use eta_mult when available (equal to eta_init / eta0 in these run files),
# avoiding floating-point label duplication.
if("eta_mult" %in% names(ls_df)){
  ls_df$eta_ratio <- as.numeric(ls_df$eta_mult)
} else {
  ls_df$eta_ratio <- ls_df$eta_init / ls_df$eta0
}
ratio_to_label <- function(x){
  out <- ifelse(abs(x - 0.2) < 1e-12, "1/5",
         ifelse(abs(x - 1.0) < 1e-12, "1",
         ifelse(abs(x - 5.0) < 1e-12, "5",
         ifelse(abs(x - 10.0) < 1e-12, "10",
                format(x, trim = TRUE, scientific = FALSE)))))
  out
}
ls_df$eta_ratio_label <- ratio_to_label(ls_df$eta_ratio)
desired_order <- c("1/5", "1", "5", "10")
ls_df$eta_ratio_f <- factor(ls_df$eta_ratio_label, levels = desired_order)
family_levels <- c("bernoulli", "poisson", "gaussian")
ls_df$family <- factor(ls_df$family, levels = family_levels)

cat(sprintf("Loaded %d line-search rows from %d files (n=1000).\n", nrow(ls_df), length(run_files)))

for(fam in family_levels){
  sub <- ls_df[ls_df$family == fam, , drop = FALSE]
  if(nrow(sub) == 0) next

  # Quick diagnostic to show whether there is replication-level variation.
  by_ratio <- split(sub$bt_mean, sub$eta_ratio_f)
  cat(sprintf("\nFamily: %s\n", fam))
  for(rn in names(by_ratio)){
    xv <- by_ratio[[rn]]
    cat(sprintf("  eta ratio %s -> unique bt_mean count = %d, value(s): %s\n",
                rn, length(unique(xv)), paste(sort(unique(xv)), collapse = ", ")))
  }

  p <- ggplot(sub, aes(x = eta_ratio_f, y = bt_mean)) +
    geom_boxplot(width = 0.62) +
    scale_y_continuous(
      breaks = function(lim){
        seq(from = floor(lim[1]), to = ceiling(lim[2]), by = 1)
      }
    ) +
    labs(
      x = expression(eta[init] / eta[0]),
      y = "Backtracking Steps",
      title = NULL
    ) +
    theme_bw(base_size = 26) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "gray92", linewidth = 0.25),
      axis.title = element_text(size = 26),
      axis.title.y = element_text(
        angle = 90,
        vjust = 0,
        margin = margin(r = 19.6)
      ),
      axis.text = element_text(size = 26),
      plot.title = element_blank(),
      plot.margin = margin(2, 2, 0, 2)
    )

  out_pdf <- file.path(out_dir, sprintf("boxplot_avg_backtracks_n1000_%s.pdf", fam))
  ggsave(out_pdf, p, width = 7.2, height = 4.5, units = "in", device = "pdf", useDingbats = FALSE)
  cat("Saved:\n", out_pdf, "\n", sep = "")
}
