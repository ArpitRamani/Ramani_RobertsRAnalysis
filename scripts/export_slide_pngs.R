# Render every analysis figure to PNG in <results_dir>/slides_png/ for slide assembly.
suppressWarnings(suppressMessages(library(limma)))
ROOT <- local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
                normalizePath(if (length(f)) file.path(dirname(f), "..") else ".") })
source(file.path(ROOT, "R", "plotting_functions.R"))

args        <- commandArgs(trailingOnly = TRUE)
results_dir <- if (length(args) >= 1) args[1] else file.path(ROOT, "results/McEachin_DE_paper")
png_dir     <- file.path(results_dir, "slides_png")
dir.create(png_dir, showWarnings = FALSE, recursive = TRUE)
P <- function(f) file.path(png_dir, f)

mat     <- as.matrix(read.delim(file.path(results_dir, "cleanDat_normalized_imputed.tsv"),
                                row.names = 1, check.names = FALSE))
samples <- read.csv(file.path(results_dir, "samples.csv"), check.names = FALSE)
de_files <- list.files(results_dir, pattern = "^DE_.*\\.csv$", full.names = TRUE)
de_files <- de_files[!grepl("summary", de_files)]

# Volcano + heatmap + p-concordance, per contrast.
for (f in de_files) {
  nm <- sub("^DE_", "", sub("\\.csv$", "", basename(f)))
  de <- read.csv(f, check.names = FALSE)
  ggplot2::ggsave(P(sprintf("volcano_%s.png", nm)),
                  plot_volcano(de, title = gsub("_", " ", nm)), width = 8, height = 6, dpi = 150)
  sig <- de$ID[which(de$P.Value < 0.05)]
  plot_heatmap(mat, samples, sig, file = P(sprintf("heatmap_%s.png", nm)),
               title = sprintf("%s (p<0.05, n=%d)", gsub("_", " ", nm), length(sig)))
  if (all(c("P.Value", "welch_P") %in% names(de)))
    ggplot2::ggsave(P(sprintf("concordance_%s.png", nm)),
                    plot_p_concordance(de, title = gsub("_", " ", nm)), width = 7.5, height = 6.5, dpi = 150)
}

# Combined STRING network + enrichment-bar figures, per contrast.
ppi_dir  <- file.path(results_dir, "PPI")
PATHWAYS <- c("KEGG", "RCTM", "WikiPathways")
for (ef in list.files(ppi_dir, pattern = "_enrichment\\.tsv$", full.names = TRUE)) {
  nm  <- sub("_enrichment\\.tsv$", "", basename(ef))
  enr <- tryCatch(read.delim(ef, check.names = FALSE), error = function(e) NULL)
  png_file <- file.path(ppi_dir, sprintf("%s_network.png", nm))
  if (is.null(enr) || !nrow(enr) || !file.exists(png_file)) next
  pbars <- plot_enrichment_bars(enr, PATHWAYS, top_n = 15,
                                title = "Enriched pathways (KEGG / Reactome / WikiPathways)")
  if (!is.null(pbars))
    combine_net_bars(png_file, pbars, P(sprintf("PPI_%s_pathways.png", nm)),
                     title = gsub("_", " ", nm), dpi = 150)
  gbars <- plot_enrichment_bars(enr, "Process", top_n = 15, title = "GO Biological Process")
  if (!is.null(gbars))
    combine_net_bars(png_file, gbars, P(sprintf("PPI_%s_GO-BP.png", nm)),
                     title = gsub("_", " ", nm), dpi = 150)
}

# limma vs Welch: variance moderation + method-count bar.
cd  <- mat[, samples$sample, drop = FALSE]
des <- model.matrix(~0 + factor(samples$group)); colnames(des) <- make.names(colnames(des))
fit <- limma::eBayes(limma::lmFit(cd, des))
ggplot2::ggsave(P("variance_moderation.png"),
                plot_variance_moderation(fit$sigma^2, fit$s2.post, fit$s2.prior),
                width = 7, height = 6, dpi = 150)
conc <- file.path(results_dir, "comparison", "limma_vs_welch_concordance.csv")
if (file.exists(conc))
  ggplot2::ggsave(P("method_significance_counts.png"),
                  plot_method_counts(read.csv(conc, check.names = FALSE)),
                  width = 7, height = 5, dpi = 150)
