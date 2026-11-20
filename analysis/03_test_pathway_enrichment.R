# Chapters 2–3 / Figs. 2D, S8–S9: one-sided ORA with the frozen mouse gene sets.
ora <- function(foreground, background, gene_sets) {
  foreground <- intersect(foreground, background)
  pieces <- lapply(names(gene_sets), function(term) {
    genes <- intersect(gene_sets[[term]], background)
    overlap <- intersect(foreground, genes)
    a <- length(overlap)
    g <- length(genes)
    f <- length(foreground)
    if (a < 2L || g < 5L) return(NULL)
    fit <- fisher.test(matrix(c(a, f - a, g - a, length(background) - g - f + a), 2L,
                              byrow = TRUE), alternative = "greater")
    data.table(Term = term, SetSize = g, ForegroundSize = f, Overlap = a,
               OddsRatio = unname(fit$estimate), P.Value = fit$p.value)
  })
  ans <- rbindlist(pieces)
  if (nrow(ans)) ans[, FDR := p.adjust(P.Value, "BH")]
  ans
}
resources <- readRDS(file.path(inputs, "enrichment_inputs.rds"))
expected <- read_ref("Tissue_and_ageing_enrichment.csv")
groups <- unique(expected[, .(Analysis, Collection, Tissue, Omics, Direction)])
all_results <- vector("list", nrow(groups))
gene_set_list <- lapply(resources$sets, function(sets) lapply(sets, function(genes) unique(toupper(genes))))
for (i in seq_len(nrow(groups))) {
  g <- groups[i]
  subset <- resources$contrasts[[g$Analysis]][
    Tissue == sub("^Skeletal muscle$|^Muscle$", "Mus", g$Tissue) & Omics == g$Omics & Testable]
  genes <- toupper(subset$GeneUpper)
  keep <- !is.na(genes) & nzchar(genes)
  direction <- if (g$Direction == "Up") 1 else -1
  foreground <- which(keep & subset$P.Value < .05 & direction * subset$logFC >= .5)
  result <- ora(unique(genes[foreground]), unique(genes[keep]), gene_set_list[[g$Collection]])
  if (nrow(result)) all_results[[i]] <- cbind(g[rep(1L, nrow(result))], result)
}
result <- rbindlist(all_results)
keys <- c(names(groups), "Term")
compare_by("Pathway enrichment", result, expected, keys,
           c("SetSize", "ForegroundSize", "Overlap", "P.Value", "FDR", "OddsRatio"))
stopifnot(nrow(result) == nrow(expected))
write_table(result, "03_functional_context/over_representation.tsv")

# S7: use the induced STRING graphs and replay the original Louvain restarts.
graphs <- c(outer(tissues, c("corrected_SNO", "corrected_SSG"),
                  function(tissue, layer) paste(layer, tissue, sep = "__")))
for (i in seq_along(graphs)) {
  base <- file.path("string", paste0("functional__", graphs[i]))
  vertices <- read_input(paste0(base, "_nodes.tsv"))
  edges <- read_input(paste0(base, "_edges.tsv"))
  g <- igraph::graph_from_data_frame(edges[, .(from, to, weight)], directed = FALSE,
                                    vertices = data.frame(name = sort(vertices$stringId)))
  set.seed(20260926L)
  seed <- 20260926L + (16L + i) * 30L
  fits <- lapply(seq_len(20L), function(start) {
    set.seed(seed + start)
    igraph::cluster_louvain(g, weights = igraph::edge_attr(g, "weight"), resolution = 1)
  })
  best <- fits[[which.max(vapply(fits, igraph::modularity, numeric(1)))]]
  membership <- igraph::membership(best)[vertices$stringId]
  check_numbers(paste("S7 Louvain adjusted Rand index", graphs[i]),
                igraph::compare(membership, vertices$community, method = "adjusted.rand"), 1)
}

# S7 PTM composition of the eight displayed modules.
network_features <- read_ref("Ageing_displayed_module_features.csv")
network_ref <- read_ref("Ageing_displayed_module_summary.csv")
module_counts <- network_features[, .(N_selected_features = .N, N_gain = sum(logFC > 0),
                                       N_loss = sum(logFC < 0)),
                                  by = .(Figure_panel, Module, Tissue, Layer)]
compare_by("S7 displayed modules", module_counts, network_ref,
           c("Figure_panel", "Module", "Tissue", "Layer"),
           c("N_selected_features", "N_gain", "N_loss"), 0)
write_table(module_counts, "03_functional_context/displayed_module_counts.tsv")
