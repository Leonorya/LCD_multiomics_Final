# Figure 4 and Tables S7: cross-omics network of the most time-responsive features.
# Input : results/preprocessed/<layer>.rds, results/time_model/<layer>.rds, data/microbiome/species_kegg_pathways.csv
# Output: results/network/cross_omics.rds
source("analysis/utils.R")

# ---------------- Part 1: settings ----------------
TOP_N <- 50            # features with the smallest P for the overall time effect per layer
MIN_MODULE_SIZE <- 3   # modules reported and displayed

# ---------------- Part 2: feature selection and partial repeated-measures correlation (Table S7.1) ----------------
features <- rbind_list(lapply(LAYERS, function(layer) {
  tm <- load_time(layer)
  tm <- head(tm[tm$time_responsive, ][order(tm$p_global[tm$time_responsive]), ], TOP_N)
  data.frame(layer = layer, feature_id = tm$feature_id, feature_name = tm$feature_name, beta_w24 = tm$beta_w24,
             node = paste(layer, tm$feature_id, sep = "::"))
}))
wide <- wide_feature_table(features)
tests <- correlate_pairs(wide, features$node)   # one BH adjustment over all feature pairs
tests$layer1 <- features$layer[match(tests$node1, features$node)]
tests$layer2 <- features$layer[match(tests$node2, features$node)]
tests$cross_layer <- tests$layer1 != tests$layer2
edges <- tests[is_sig(tests$fdr), ]

# ---------------- Part 3: Louvain modules and cross-layer connectivity (Figure 4A-D; Tables S7.2-S7.3) ----------------
graph <- igraph::graph_from_data_frame(data.frame(from = edges$node1, to = edges$node2, weight = abs(edges$r)),
                                       directed = FALSE, vertices = features$node)
set.seed(SEED)
features$louvain <- as.integer(igraph::membership(igraph::cluster_louvain(graph, weights = igraph::E(graph)$weight))[features$node])

# Modules with at least MIN_MODULE_SIZE connected members, numbered by size.
connected <- features$node %in% c(edges$node1, edges$node2)
sizes <- sort(table(features$louvain[connected]), decreasing = TRUE)
kept <- as.integer(names(sizes)[sizes >= MIN_MODULE_SIZE])
features$module <- ifelse(connected & features$louvain %in% kept, paste0("Module", match(features$louvain, kept)), NA)

# Cross-layer connectivity: sum of |r| over a feature's edges with features from other layers.
cross <- edges[edges$cross_layer, ]
strength <- tapply(abs(c(cross$r, cross$r)), c(cross$node1, cross$node2), sum)
features$cross_layer_strength <- ifelse(is.na(strength[features$node]), 0, strength[features$node])

# ---------------- Part 4: KEGG analysis of cross-layer modules (Figure 4E; Table S7.4) ----------------
# Modules spanning at least two layers; the members of each layer form one feature set.
species_pathways <- read.csv(file.path(DATA_DIR, "microbiome", "species_kegg_pathways.csv"))
module_layers <- tapply(features$layer, features$module, function(x) length(unique(x)))
cross_modules <- names(module_layers)[module_layers >= 2]
module_enrichment <- rbind_list(lapply(setdiff(LAYERS, "fecal_metatranscriptome"), function(layer) {
  members <- features[features$layer == layer & features$module %in% cross_modules, ]
  if (nrow(members) < 2) return(NULL)
  prep <- load_layer(layer)
  if (layer == "fecal_metabolome" && is.null(prep$annotation$kegg_id)) prep$annotation$kegg_id <- match_kegg_by_name(prep$annotation$feature_name)
  kegg_ora_layer(layer, split(members$feature_id, members$module), colnames(prep$z), prep$annotation, species_pathways)
}))

saveRDS(list(features = features, tests = tests, edges = edges, module_enrichment = module_enrichment),
        result_path("network", "cross_omics.rds"))
