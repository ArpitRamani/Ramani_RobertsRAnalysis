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

# Pad a raster array to a target width/height aspect with white margins (no crop,
# no node loss), so it can fill a panel of that aspect without distortion.
pad_to_aspect <- function(img, target) {            # target = width / height
  d <- dim(img); h <- d[1]; w <- d[2]; cur <- w / h
  ch <- if (length(d) == 3) d[3] else 1L
  if (abs(cur - target) < 1e-3) return(img)
  if (cur < target) {                               # too narrow -> add width
    nw <- round(h * target); x0 <- (nw - w) %/% 2
    canvas <- array(1, dim = c(h, nw, ch))
    if (length(d) == 3) canvas[, (x0 + 1):(x0 + w), ] <- img else canvas[, (x0 + 1):(x0 + w)] <- img
  } else {                                          # too wide -> add height
    nh <- round(w / target); y0 <- (nh - h) %/% 2
    canvas <- array(1, dim = c(nh, w, ch))
    if (length(d) == 3) canvas[(y0 + 1):(y0 + h), , ] <- img else canvas[(y0 + 1):(y0 + h), ] <- img
  }
  canvas
}

# Legend for STRING evidence-flavor edge colors (what the lines in the network mean).
string_edge_legend <- function() {
  key <- data.frame(
    label = c("Curated database", "Experimental", "Co-expression", "Text-mining",
              "Neighborhood", "Gene fusion", "Co-occurrence", "Homology"),
    col   = c("#00B7EB", "#D63CC8", "#000000", "#A6CE39",
              "#37B34A", "#ED1C24", "#2E3192", "#B6B6E8"),
    stringsAsFactors = FALSE)
  key$x <- rep(c(0, 1, 2, 3), each = 2)
  key$y <- rep(c(1, 0), times = 4)
  ggplot(key) +
    geom_segment(aes(x = x, xend = x + 0.16, y = y, yend = y, colour = col), linewidth = 1.7) +
    geom_text(aes(x = x + 0.20, y = y, label = label), hjust = 0, size = 3) +
    scale_colour_identity() +
    scale_x_continuous(limits = c(0, 4)) + scale_y_continuous(limits = c(-0.6, 1.6)) +
    labs(title = "Edge evidence (STRING)") +
    theme_void(base_size = 10) +
    theme(plot.title = element_text(face = "bold", size = 9.5),
          plot.margin = margin(2, 6, 2, 6))
}

# Network PNG (left, sized to its own aspect) + enrichment bars (right) + edge legend (bottom).
combine_net_bars <- function(png_file, bar_gg, out_file, title = "") {
  img <- png::readPNG(png_file)
  ih  <- dim(img)[1]; iw <- dim(img)[2]

  net_w <- 10; bars_w <- 6; leg_h <- 1.2                 # network gets the larger column
  net_h <- max(5, min(net_w * ih / iw, 13))              # clamp portrait/landscape extremes
  img   <- pad_to_aspect(img, net_w / net_h)             # pad to the panel aspect -> no stretch
  g_img <- grid::rasterGrob(img, interpolate = TRUE,
                            width = grid::unit(1, "npc"), height = grid::unit(1, "npc"))

  top <- patchwork::wrap_elements(full = g_img) + bar_gg + patchwork::plot_layout(widths = c(net_w, bars_w))
  p   <- top / string_edge_legend() + patchwork::plot_layout(heights = c(net_h, leg_h)) +
         patchwork::plot_annotation(title = title,
                                    theme = theme(plot.title = element_text(face = "bold")))
  ggsave(out_file, p, width = net_w + bars_w, height = net_h + leg_h, limitsize = FALSE)
}
