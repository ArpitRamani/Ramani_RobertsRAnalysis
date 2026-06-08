# Volcano + heatmap for each DE_*.csv in a results folder.
ROOT <- local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
                normalizePath(if (length(f)) file.path(dirname(f), "..") else ".") })
source(file.path(ROOT, "R", "plotting_functions.R"))

args        <- commandArgs(trailingOnly = TRUE)
results_dir <- if (length(args) >= 1) args[1] else file.path(ROOT, "results/McEachin_DE_paper")
plot_dir    <- file.path(results_dir, "plots")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

mat     <- as.matrix(read.delim(file.path(results_dir, "cleanDat_normalized_imputed.tsv"),
                                row.names = 1, check.names = FALSE))
samples <- read.csv(file.path(results_dir, "samples.csv"), check.names = FALSE)

de_files <- list.files(results_dir, pattern = "^DE_.*\\.csv$", full.names = TRUE)
de_files <- de_files[!grepl("summary", de_files)]

for (f in de_files) {
  nm <- sub("^DE_", "", sub("\\.csv$", "", basename(f)))
  de <- read.csv(f, check.names = FALSE)
  v  <- plot_volcano(de, title = gsub("_", " ", nm))
  ggplot2::ggsave(file.path(plot_dir, sprintf("volcano_%s.pdf", nm)), v, width = 8, height = 6)
  ggplot2::ggsave(file.path(plot_dir, sprintf("volcano_%s.png", nm)), v, width = 8, height = 6, dpi = 150)
  sig <- de$ID[which(de$P.Value < 0.05)]
  plot_heatmap(mat, samples, sig, file = file.path(plot_dir, sprintf("heatmap_%s.pdf", nm)),
               title = sprintf("%s (p<0.05, n=%d)", gsub("_", " ", nm), length(sig)))
}
