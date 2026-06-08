# TAMPOR cross-batch normalization (robust core; keeps sites present in both plexes).
suppressWarnings(suppressMessages(library(limma)))

ROOT <- local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
                normalizePath(if (length(f)) file.path(dirname(f), "..") else ".") })
source(file.path(ROOT, "R", "nitro_DE_functions.R"))

cfg <- list(
  mode        = "tampor",
  concordance = FALSE,
  site_file   = file.path(ROOT, "data/McEachin_single-site_nitro_snapshot.tsv"),
  meta_file   = file.path(ROOT, "data/metadata/original_metaData.csv"),
  tampor_src  = file.path(ROOT, "R/vendor/TAMPOR.R"),
  out_dir     = file.path(ROOT, "results/McEachin_DE_TAMPOR"),
  min_present = 0.50,
  impute      = TRUE,
  ref_channel = "126",
  contrasts   = list(
    c9ALS_vs_Control = c("c9ALS", "Control"),
    ASO_vs_Control   = c("ASO",   "Control"),
    c9ALS_vs_ASO     = c("c9ALS", "ASO")
  )
)
run_pipeline(cfg)
