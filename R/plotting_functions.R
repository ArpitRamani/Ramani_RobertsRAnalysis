# Volcano + heatmap helpers.
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
