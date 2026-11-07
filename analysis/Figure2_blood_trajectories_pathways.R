# Figure 2 and Tables S5: trajectory clusters, pathway enrichment and protein interactions of the blood omics layers.
# Input : results/preprocessed/<layer>.rds, results/time_model/<layer>.rds
# Output: results/navmix/<layer>.rds, results/enrichment/clusters_<layer>.rds, results/all_visit_summary_blood.rds,
#         results/ppi/serum_proteome.rds
source("analysis/utils.R")

# ---------------- Part 1: settings ----------------
STRING_VERSION   <- "12.0"
STRING_MIN_SCORE <- 400   # combined score (medium confidence)

# ---------------- Part 2: NAvMix trajectory clustering (Figure 2A-C; Tables S5.1-S5.3) ----------------
for (layer in BLOOD_LAYERS) {
  nav <- run_navmix(load_layer(layer), time_responsive_ids(layer), NAVMIX_K[[layer]], SEED + match(layer, LAYERS))
  saveRDS(nav, result_path("navmix", paste0(layer, ".rds")))
}

# ---------------- Part 3: features changed at all three post-baseline visits (Figure 2D) ----------------
saveRDS(rbind_list(lapply(BLOOD_LAYERS, all_visit_summary)), result_path("all_visit_summary_blood.rds"))

# ---------------- Part 4: KEGG over-representation analysis of trajectory clusters (Figure 2E; Table S5.4) ----------------
for (layer in BLOOD_LAYERS) {
  prep <- load_layer(layer)
  a <- read_result("navmix", paste0(layer, ".rds"))$assignments
  a <- a[a$cluster != "Noise", ]
  enrichment <- kegg_ora_layer(layer, split(a$feature_id, a$cluster), colnames(prep$z), prep$annotation)
  saveRDS(enrichment, result_path("enrichment", paste0("clusters_", layer, ".rds")))
}

# ---------------- Part 5: STRING interactions and k-core decomposition of time-responsive proteins (Figure 2F; Tables S5.5-S5.6) ----------------
prep <- load_layer("serum_proteome")
tm <- load_time("serum_proteome")
proteins <- data.frame(feature_id = tm$feature_id[tm$time_responsive], beta_w24 = tm$beta_w24[tm$time_responsive])
proteins$gene_symbol <- prep$annotation$gene_symbol[match(proteins$feature_id, prep$annotation$feature_id)]

string_db <- STRINGdb::STRINGdb$new(version = STRING_VERSION, species = 9606, score_threshold = STRING_MIN_SCORE, input_directory = tempdir())
mapped <- string_db$map(proteins, "gene_symbol", removeUnmappedRows = TRUE, takeFirst = TRUE)
edges <- string_db$get_interactions(unique(mapped$STRING_id))
graph <- igraph::simplify(igraph::graph_from_data_frame(edges[, c("from", "to", "combined_score")], directed = FALSE,
                                                        vertices = unique(mapped$STRING_id)),
                          edge.attr.comb = list(combined_score = "max"))
nodes <- data.frame(string_id = igraph::V(graph)$name, degree = igraph::degree(graph), core = igraph::coreness(graph))
nodes <- merge(nodes, mapped[!duplicated(mapped$STRING_id), c("STRING_id", "feature_id", "gene_symbol", "beta_w24")],
               by.x = "string_id", by.y = "STRING_id")
nodes$highest_core <- nodes$core == max(nodes$core)
saveRDS(list(nodes = nodes, edges = igraph::as_data_frame(graph, what = "edges")), result_path("ppi", "serum_proteome.rds"))
