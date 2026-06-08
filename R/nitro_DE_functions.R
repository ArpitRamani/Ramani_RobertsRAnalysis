# DE pipeline: load -> normalize -> missingness filter -> impute -> limma (+ optional Welch).
# Modes: "paper" (log2 only), "tampor" (cross-batch median polish), "fallback" (log2 + median center).

normalize_paper <- function(abund) {
  m <- log2(abund); m[!is.finite(m)] <- NA; m
}

normalize_fallback <- function(abund) {
  m <- log2(abund); m[!is.finite(m)] <- NA
  cm <- apply(m, 2, median, na.rm = TRUE)
  sweep(m, 2, cm - median(cm, na.rm = TRUE))
}

normalize_tampor <- function(abund, meta, src, out_dir) {
  source(src)
  traits <- data.frame(Batch = meta$Batch, row.names = meta$sample, stringsAsFactors = FALSE)
  res <- TAMPOR(dat = as.data.frame(abund), traits = traits, noGIS = TRUE,
                batchPrefixInSampleNames = FALSE, parallelThreads = 2,
                outputSuffix = "McEachin_nitro", path = out_dir)
  as.matrix(res$cleanDat)
}

welch_contrast <- function(mat, meta, case, ref) {
  gi <- which(meta$group == case); ri <- which(meta$group == ref)
  r <- t(apply(mat, 1, function(x) {
    fc <- mean(x[gi], na.rm = TRUE) - mean(x[ri], na.rm = TRUE)
    p  <- tryCatch(t.test(x[gi], x[ri])$p.value, error = function(e) NA_real_)
    c(welch_logFC = fc, welch_P = p)
  }))
  data.frame(ID = rownames(mat), r, welch_adjP = p.adjust(r[, "welch_P"], "BH"),
             row.names = NULL, check.names = FALSE)
}

run_pipeline <- function(cfg) {
  stopifnot(cfg$mode %in% c("paper", "tampor", "fallback"))
  dir.create(cfg$out_dir, showWarnings = FALSE, recursive = TRUE)

  raw   <- read.delim(cfg$site_file, check.names = FALSE, stringsAsFactors = FALSE)
  chan  <- grep("^b[0-9]+_", names(raw), value = TRUE)
  abund <- as.matrix(raw[, chan]); storage.mode(abund) <- "double"
  rownames(abund) <- raw$Index
  annot <- raw[, intersect(c("Index","Gene","ProteinID","Peptide"), names(raw)), drop = FALSE]
  rownames(annot) <- raw$Index

  meta <- read.csv(cfg$meta_file, check.names = FALSE, stringsAsFactors = FALSE)
  meta$sample <- paste0(meta$Batch, "_", meta$Channel)
  meta$group  <- meta$DiseaseGroup
  meta <- meta[meta$Channel != cfg$ref_channel, ]
  keep <- intersect(meta$sample, colnames(abund))
  meta <- meta[match(keep, meta$sample), ]
  abund <- abund[, keep, drop = FALSE]

  if (cfg$mode == "paper") {
    cleanDat <- normalize_paper(abund);    batch_in_model <- FALSE
  } else if (cfg$mode == "fallback") {
    cleanDat <- normalize_fallback(abund); batch_in_model <- TRUE
  } else {
    cleanDat <- normalize_tampor(abund, meta, cfg$tampor_src, cfg$out_dir)
    batch_in_model <- FALSE
  }
  cleanDat <- cleanDat[, meta$sample, drop = FALSE]
  cleanDat <- cleanDat[rowMeans(!is.na(cleanDat)) >= cfg$min_present, , drop = FALSE]

  if (isTRUE(cfg$impute) && anyNA(cleanDat)) {
    set.seed(42)
    fit <- missForest::missForest(as.data.frame(t(cleanDat)), maxiter = 5, ntree = 100)
    cleanDat <- t(as.matrix(fit$ximp))
  }

  write.table(data.frame(ID = rownames(cleanDat), cleanDat, check.names = FALSE),
              file.path(cfg$out_dir, "cleanDat_normalized_imputed.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  write.csv(meta[, intersect(c("sample","group","Batch","Sex"), names(meta))],
            file.path(cfg$out_dir, "samples.csv"), row.names = FALSE)

  grp <- factor(meta$group); dat_terms <- list(grp = grp); form <- "~0 + grp"
  if (batch_in_model) { dat_terms$batch <- factor(meta$Batch); form <- "~0 + grp + batch" }
  design <- model.matrix(as.formula(form), data = dat_terms)
  colnames(design) <- make.names(colnames(design))
  fit <- limma::lmFit(cleanDat, design)

  limma_contrast <- function(case, ref) {
    cn <- make.names(paste0("grp", c(case, ref)))
    cm <- limma::makeContrasts(contrasts = paste(cn[1], "-", cn[2]), levels = design)
    f2 <- limma::eBayes(limma::contrasts.fit(fit, cm))
    tt <- limma::topTable(f2, number = Inf, sort.by = "none")
    tt$ID <- rownames(tt); tt$Gene <- annot[tt$ID, "Gene"]
    tt$minus_log10_p <- -log10(tt$P.Value)
    tt[, c("ID","Gene","logFC","AveExpr","t","P.Value","adj.P.Val","minus_log10_p")]
  }

  summ <- list()
  for (nm in names(cfg$contrasts)) {
    cc  <- cfg$contrasts[[nm]]
    out <- limma_contrast(cc[1], cc[2])
    if (isTRUE(cfg$concordance)) {
      w   <- welch_contrast(cleanDat, meta, cc[1], cc[2])
      out <- merge(out, w[, c("ID","welch_P","welch_adjP")], by = "ID", all.x = TRUE, sort = FALSE)
    }
    out <- out[order(out$P.Value), ]
    write.csv(out, file.path(cfg$out_dir, sprintf("DE_%s.csv", nm)), row.names = FALSE)
    summ[[nm]] <- data.frame(contrast = nm, sites = nrow(out),
                             limma_p05 = sum(out$P.Value < 0.05, na.rm = TRUE),
                             limma_FDR05 = sum(out$adj.P.Val < 0.05, na.rm = TRUE))
  }
  summ <- do.call(rbind, summ)
  write.csv(summ, file.path(cfg$out_dir, "DE_summary.csv"), row.names = FALSE)
  invisible(summ)
}
