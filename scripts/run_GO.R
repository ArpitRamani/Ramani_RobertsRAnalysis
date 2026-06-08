# GO:BP enrichment per contrast (sig = p<0.05; background = all tested genes).
ROOT <- local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
                normalizePath(if (length(f)) file.path(dirname(f), "..") else ".") })
source(file.path(ROOT, "R", "enrichment_functions.R"))

args        <- commandArgs(trailingOnly = TRUE)
results_dir <- if (length(args) >= 1) args[1] else file.path(ROOT, "results/McEachin_DE_paper")
go_dir      <- file.path(results_dir, "GO")
dir.create(go_dir, showWarnings = FALSE, recursive = TRUE)

de_files <- list.files(results_dir, pattern = "^DE_.*\\.csv$", full.names = TRUE)
de_files <- de_files[!grepl("summary", de_files)]

for (f in de_files) {
  nm  <- sub("^DE_", "", sub("\\.csv$", "", basename(f)))
  de  <- read.csv(f, check.names = FALSE)
  ego <- tryCatch(run_go(de$Gene[which(de$P.Value < 0.05)], de$Gene, ont = "BP"),
                  error = function(e) NULL)
  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) next
  write.csv(as.data.frame(ego), file.path(go_dir, sprintf("GO_BP_%s.csv", nm)), row.names = FALSE)
  p <- enrichplot::dotplot(ego, showCategory = 15) + ggplot2::ggtitle(paste("GO:BP", gsub("_", " ", nm)))
  ggplot2::ggsave(file.path(go_dir, sprintf("GO_BP_%s_dotplot.pdf", nm)), p, width = 8, height = 7)
}
