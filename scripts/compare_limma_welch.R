# limma (moderated t) vs Welch t-test on the same sites: concordance table + plots.
# Reads DE_*.csv (both p-value columns) and the saved cleanDat to expose limma's
# empirical-Bayes variance moderation -- the actual statistical difference between the two.
suppressWarnings(suppressMessages(library(limma)))
ROOT <- local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
                normalizePath(if (length(f)) file.path(dirname(f), "..") else ".") })
source(file.path(ROOT, "R", "plotting_functions.R"))

args        <- commandArgs(trailingOnly = TRUE)
results_dir <- if (length(args) >= 1) args[1] else file.path(ROOT, "results/McEachin_DE_paper")
out_dir     <- file.path(results_dir, "comparison")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Refit limma on the saved matrix to recover raw vs moderated per-site variance.
cd  <- as.matrix(read.delim(file.path(results_dir, "cleanDat_normalized_imputed.tsv"),
                            check.names = FALSE, row.names = 1))
sm  <- read.csv(file.path(results_dir, "samples.csv"), check.names = FALSE)
cd  <- cd[, sm$sample, drop = FALSE]
grp <- factor(sm$group)
des <- model.matrix(~0 + grp); colnames(des) <- make.names(colnames(des))
fit <- limma::eBayes(limma::lmFit(cd, des))
vm  <- plot_variance_moderation(fit$sigma^2, fit$s2.post, fit$s2.prior)
ggplot2::ggsave(file.path(out_dir, "variance_moderation.pdf"), vm, width = 7, height = 6)

# Per-contrast p-value concordance + a summary row each.
de_files <- list.files(results_dir, pattern = "^DE_.*\\.csv$", full.names = TRUE)
de_files <- de_files[!grepl("summary", de_files)]

rows <- list()
for (f in de_files) {
  nm <- sub("^DE_", "", sub("\\.csv$", "", basename(f)))
  de <- read.csv(f, check.names = FALSE)
  if (!all(c("P.Value", "welch_P") %in% names(de))) next
  d  <- de[is.finite(de$P.Value) & is.finite(de$welch_P), ]

  ggplot2::ggsave(file.path(out_dir, sprintf("%s_p_concordance.pdf", nm)),
                  plot_p_concordance(d, title = gsub("_", " ", nm)), width = 7.5, height = 6.5)

  rows[[nm]] <- data.frame(
    contrast            = nm,
    n                   = nrow(d),
    limma_sig           = sum(d$P.Value < 0.05),
    welch_sig           = sum(d$welch_P < 0.05),
    both_sig            = sum(d$P.Value < 0.05 & d$welch_P < 0.05),
    limma_only          = sum(d$P.Value < 0.05 & d$welch_P >= 0.05),
    welch_only          = sum(d$P.Value >= 0.05 & d$welch_P < 0.05),
    spearman_p          = cor(d$P.Value, d$welch_P, method = "spearman"),
    pearson_neglog10p   = cor(-log10(d$P.Value), -log10(d$welch_P)),
    median_abs_dneglogp = median(abs(log10(d$welch_P) - log10(d$P.Value))),
    pct_limma_more_sig  = round(mean(d$P.Value < d$welch_P) * 100, 1))
}
summ <- do.call(rbind, rows)
write.csv(summ, file.path(out_dir, "limma_vs_welch_concordance.csv"), row.names = FALSE)
ggplot2::ggsave(file.path(out_dir, "method_significance_counts.pdf"),
                plot_method_counts(summ), width = 7, height = 5)
