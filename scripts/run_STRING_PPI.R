# STRING PPI per contrast (sig = p<0.05): network image, edge table, enrichment.
ROOT <- local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
                normalizePath(if (length(f)) file.path(dirname(f), "..") else ".") })
source(file.path(ROOT, "R", "enrichment_functions.R"))

args        <- commandArgs(trailingOnly = TRUE)
results_dir <- if (length(args) >= 1) args[1] else file.path(ROOT, "results/McEachin_DE_paper")
ppi_dir     <- file.path(results_dir, "PPI")
dir.create(ppi_dir, showWarnings = FALSE, recursive = TRUE)

de_files <- list.files(results_dir, pattern = "^DE_.*\\.csv$", full.names = TRUE)
de_files <- de_files[!grepl("summary", de_files)]

for (f in de_files) {
  nm <- sub("^DE_", "", sub("\\.csv$", "", basename(f)))
  de <- read.csv(f, check.names = FALSE)
  string_ppi(unique(de$Gene[which(de$P.Value < 0.05)]), ppi_dir, label = nm,
             species = 9606, required_score = 400)
}
