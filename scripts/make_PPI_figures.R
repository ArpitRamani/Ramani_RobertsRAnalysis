# Combined STRING network + enrichment-bar figures, per contrast.
# Reads the PPI/*_enrichment.tsv + *_network.png written by run_STRING_PPI.R.
ROOT <- local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
                normalizePath(if (length(f)) file.path(dirname(f), "..") else ".") })
source(file.path(ROOT, "R", "plotting_functions.R"))

args        <- commandArgs(trailingOnly = TRUE)
results_dir <- if (length(args) >= 1) args[1] else file.path(ROOT, "results/McEachin_DE_paper")
ppi_dir     <- file.path(results_dir, "PPI")

PATHWAYS <- c("KEGG", "RCTM", "WikiPathways")
GO       <- list(Process = "GO Biological Process",
                 Function = "GO Molecular Function",
                 Component = "GO Cellular Component")

enr_files <- list.files(ppi_dir, pattern = "_enrichment\\.tsv$", full.names = TRUE)

for (ef in enr_files) {
  nm  <- sub("_enrichment\\.tsv$", "", basename(ef))
  enr <- tryCatch(read.delim(ef, check.names = FALSE), error = function(e) NULL)
  if (is.null(enr) || !nrow(enr)) next
  png_file <- file.path(ppi_dir, sprintf("%s_network.png", nm))

  pbars <- plot_enrichment_bars(enr, PATHWAYS, top_n = 15,
                                title = "Enriched pathways (KEGG / Reactome / WikiPathways)")
  if (!is.null(pbars) && file.exists(png_file))
    combine_net_bars(png_file, pbars, file.path(ppi_dir, sprintf("%s_network_pathways.pdf", nm)),
                     title = gsub("_", " ", nm))

  gbars <- plot_enrichment_bars(enr, "Process", top_n = 15, title = "GO Biological Process")
  if (!is.null(gbars) && file.exists(png_file))
    combine_net_bars(png_file, gbars, file.path(ppi_dir, sprintf("%s_network_GO-BP.pdf", nm)),
                     title = gsub("_", " ", nm))

  for (cat in names(GO)) {
    b <- plot_enrichment_bars(enr, cat, top_n = 15, title = paste(gsub("_", " ", nm), "-", GO[[cat]]))
    if (!is.null(b)) ggplot2::ggsave(file.path(ppi_dir, sprintf("%s_%s_bars.pdf", nm, cat)),
                                     b, width = 8, height = 6)
  }
}
