# Volcano + heatmap + STRING-enrichment bars + network/enrichment combiner.
suppressWarnings(suppressMessages({ library(ggplot2); library(ggrepel) }))

plot_volcano <- function(de, title = "", p_thresh = 0.05, fc_thresh = 0, top_n = 25) {
  d <- de[is.finite(de$logFC) & is.finite(de$P.Value), ]
  d$status <- with(d, ifelse(P.Value < p_thresh & logFC >  fc_thresh, "Increased",
                      ifelse(P.Value < p_thresh & logFC < -fc_thresh, "Decreased", "Unchanged")))
  d$status <- factor(d$status, levels = c("Decreased", "Unchanged", "Increased"))
  lab <- head(d[order(d$P.Value), ], top_n)
  ggplot(d, aes(logFC, -log10(P.Value), color = status)) +
    geom_hline(yintercept = -log10(p_thresh), linetype = 2, linewidth = 0.3) +
    geom_vline(xintercept = c(-fc_thresh, fc_thresh), linetype = 3, linewidth = 0.3) +
    geom_point(alpha = 0.85, size = 1.8) +
    ggrepel::geom_text_repel(data = lab, aes(label = Gene), size = 3, max.overlaps = Inf,
                             box.padding = 0.3, segment.size = 0.2, show.legend = FALSE) +
    scale_color_manual(name = "Change",
                       values = c(Decreased = "#2166AC", Unchanged = "grey75", Increased = "#B2182B"),
                       breaks = c("Increased", "Decreased", "Unchanged")) +
    labs(title = title, x = "log2 fold change", y = expression(-log[10](p))) +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))
}

plot_heatmap <- function(mat, samples, sites, file, title = "") {
  sites <- intersect(sites, rownames(mat))
  if (length(sites) < 2) return(invisible())
  m   <- t(scale(t(mat[sites, , drop = FALSE])))
  ann <- data.frame(group = samples$group, row.names = samples$sample)[colnames(m), , drop = FALSE]
  ord <- order(samples$group[match(colnames(m), samples$sample)])
  pheatmap::pheatmap(m[, ord], annotation_col = ann, cluster_cols = FALSE,
                     show_rownames = length(sites) <= 60, show_colnames = FALSE,
                     main = title, fontsize_row = 6, filename = file, width = 7, height = 8)
}

# Horizontal -log10(FDR) bar chart of the top STRING-enriched terms in given categories.
# STRING categories: Process / Function / Component (GO BP/MF/CC), KEGG, RCTM (Reactome),
# WikiPathways, HPO, TISSUES, COMPARTMENTS, ...
plot_enrichment_bars <- function(enr, categories, top_n = 15, title = "") {
  d <- enr[enr$category %in% categories & is.finite(enr$fdr), , drop = FALSE]
  if (!nrow(d)) return(NULL)
  d <- d[order(d$fdr), , drop = FALSE]
  d <- d[seq_len(min(top_n, nrow(d))), ]
  wrap <- function(x) vapply(x, function(s) paste(strwrap(s, 40), collapse = "\n"), character(1))
  d$lab <- factor(wrap(d$description), levels = rev(wrap(d$description)))
  d$nlp <- -log10(d$fdr)
  ggplot(d, aes(nlp, lab, fill = nlp)) +
    geom_col(width = 0.72) +
    scale_fill_gradient(low = "#a8ddb5", high = "#2b8cbe", guide = "none") +
    labs(x = expression(-log[10](FDR)), y = NULL, title = title) +
    theme_bw(base_size = 11) +
    theme(axis.text.y = element_text(size = 8), plot.title = element_text(face = "bold", size = 11))
}

# Pad a raster array to a square with white margins (no crop, no node loss),
# so it can fill a square panel without aspect-ratio distortion.
pad_to_square <- function(img) {
  d <- dim(img); h <- d[1]; w <- d[2]
  if (h == w) return(img)
  s  <- max(h, w)
  ch <- if (length(d) == 3) d[3] else 1L
  canvas <- array(1, dim = c(s, s, ch))          # opaque white
  y0 <- (s - h) %/% 2; x0 <- (s - w) %/% 2
  if (length(d) == 3) canvas[(y0 + 1):(y0 + h), (x0 + 1):(x0 + w), ] <- img
  else                canvas[(y0 + 1):(y0 + h), (x0 + 1):(x0 + w)]   <- img
  canvas
}

# Network PNG (left, padded to square) + enrichment bar chart (right), as one figure.
combine_net_bars <- function(png_file, bar_gg, out_file, title = "") {
  img   <- pad_to_square(png::readPNG(png_file))
  g_img <- grid::rasterGrob(img, interpolate = TRUE,
                            width = grid::unit(1, "npc"), height = grid::unit(1, "npc"))
  p <- patchwork::wrap_elements(full = g_img) + bar_gg +
       patchwork::plot_layout(widths = c(1, 1)) +
       patchwork::plot_annotation(title = title,
                                  theme = theme(plot.title = element_text(face = "bold")))
  ggsave(out_file, p, width = 16, height = 8)
}
