### Script overview:
### (1) Read n=1000 run-level results from:
###     ../results for GD behavior
### (2) Keep line-search runs only (method == "linesearch").
### (3) Plot boxplots of GD steps vs eta_init / eta_0.
### (4) Create one PDF per family (bernoulli/poisson/gaussian).
### (5) Save:
###     boxplot_gd_steps_n1000_{family}.pdf


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

# Use eta_mult when available (it is eta_init / eta0 in these run files),
# avoiding floating-point label duplication from eta_init/eta0 arithmetic.
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

  p <- ggplot(sub, aes(x = eta_ratio_f, y = num_of_steps)) +
    geom_boxplot(width = 0.62) +
    labs(
      x = expression(eta[init] / eta[0]),
      y = expression(paste("Iterations ", R[conv])),
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
        margin = margin(r = 16)
      ),
      axis.text = element_text(size = 26),
      plot.title = element_blank(),
      plot.margin = margin(2, 2, 0, 2)
    )

  out_pdf <- file.path(out_dir, sprintf("boxplot_gd_steps_n1000_%s.pdf", fam))
  ggsave(out_pdf, p, width = 7.2, height = 4.5, units = "in", device = "pdf", useDingbats = FALSE)
  cat("Saved:\n", out_pdf, "\n", sep = "")
}
