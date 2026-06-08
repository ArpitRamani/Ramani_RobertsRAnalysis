# STRING protein-protein interactions + functional enrichment via the REST API.
# enrichment covers GO (Process/Function/Component), pathways (KEGG/Reactome/WikiPathways),
# HPO, and more -- filter by the `category` column.

string_ppi <- function(genes, out_dir, label, species = 9606,
                       required_score = 400, caller = "nitration_pipeline") {
  genes <- unique(na.omit(genes))
  if (length(genes) < 3) return(NULL)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  ids  <- paste(genes, collapse = "%0d")               # STRING wants CR-separated identifiers
  base <- "https://string-db.org/api"
  q    <- sprintf("identifiers=%s&species=%d&required_score=%d&caller_identity=%s",
                  ids, species, required_score, caller)
  get  <- function(endpoint, dest, binary = FALSE)
    tryCatch(utils::download.file(sprintf("%s/%s?%s", base, endpoint, q), dest,
                                  mode = if (binary) "wb" else "w", quiet = TRUE),
             error = function(e) NULL)
  get("highres_image/network", file.path(out_dir, sprintf("%s_network.png", label)), binary = TRUE)
  get("tsv/network",           file.path(out_dir, sprintf("%s_interactions.tsv", label)))
  ok <- get("tsv/enrichment",  file.path(out_dir, sprintf("%s_enrichment.tsv", label)))
  if (is.null(ok)) return(NULL)
  tryCatch(read.delim(file.path(out_dir, sprintf("%s_enrichment.tsv", label)), check.names = FALSE),
           error = function(e) NULL)
}
