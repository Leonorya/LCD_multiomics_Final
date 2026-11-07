# Figure 3 and Tables S4.2, S5 and S6: gut microbial diversity, trajectory clusters and pathway mapping of the fecal
# omics layers, and effect modification of trajectory-cluster scores of all six layers.
# Input : results/preprocessed/<layer>.rds, results/time_model/<layer>.rds, results/navmix/<blood layer>.rds,
#         data/microbiome/species_relative_abundance.csv, data/microbiome/gene_counts.csv,
#         data/microbiome/species_kegg_pathways.csv, data/medication.csv
# Output: results/diversity.rds, results/navmix/<layer>.rds, results/enrichment/clusters_<layer>.rds,
#         results/all_visit_summary_stool.rds, results/effect_modification/clusters.rds
source("analysis/utils.R")

# ---------------- Part 1: settings ----------------
N_PERMUTATIONS <- 999
MODIFIERS <- c("diet", "sex", "metformin")

# ---------------- Part 2: alpha diversity and gene count (Figure 3A; Tables S6.1-S6.2) ----------------
meta <- load_layer("fecal_metagenome")$metadata
species <- read.csv(file.path(DATA_DIR, "microbiome", "species_relative_abundance.csv"), check.names = FALSE)
gene_counts <- read.csv(file.path(DATA_DIR, "microbiome", "gene_counts.csv"))
rel <- as.matrix(species[match(meta$sample_id, species$sample_id), -1])
rel <- rel / rowSums(rel)   # classified species, re-normalized within each sample
rownames(rel) <- meta$sample_id

genes <- gene_counts[match(meta$sample_id, gene_counts$sample_id), ]
meta$shannon <- vegan::diversity(rel, index = "shannon")
meta$simpson <- vegan::diversity(rel, index = "simpson")
meta$observed_species <- rowSums(rel > 0)
meta$gene_count <- genes$gene_count
meta$assigned_gene_count <- genes$assigned_gene_count
meta$unclassified_gene_count <- genes$gene_count - genes$assigned_gene_count
meta$assigned_gene_fraction <- genes$assigned_gene_count / genes$gene_count

# Table S6.1
DESCRIBED <- c("shannon", "simpson", "observed_species", "gene_count", "assigned_gene_count", "unclassified_gene_count", "assigned_gene_fraction")
alpha_by_visit <- rbind_list(lapply(DESCRIBED, function(m) rbind_list(lapply(VISITS, function(v) { x <- meta[[m]][meta$visit == v]
  data.frame(metric = m, visit = v, mean = mean(x, na.rm = TRUE), sd = stats::sd(x, na.rm = TRUE), median = stats::median(x, na.rm = TRUE),
             q1 = stats::quantile(x, 0.25, na.rm = TRUE), q3 = stats::quantile(x, 0.75, na.rm = TRUE)) }))))

# Table S6.2
meta$log10_observed_species <- log10(meta$observed_species)
meta$log10_gene_count <- log10(meta$gene_count)
MODELLED <- c("shannon", "simpson", "log10_observed_species", "log10_gene_count")
alpha_lmm <- rbind_list(lapply(MODELLED, function(m) data.frame(metric = m, t(fit_time_lmm(meta[[m]], meta)))))
for (p in c("p_w1", "p_w12", "p_w24", "p_global")) alpha_lmm[[sub("^p_", "fdr_", p)]] <- stats::p.adjust(alpha_lmm[[p]], "BH")

# ---------------- Part 3: community composition (Tables S6.3-S6.5) ----------------
# Aitchison distance: zeros replaced by half the smallest positive abundance of the sample, then CLR.
half_min_clr <- function(x) t(apply(x, 1, function(r) { r[r <= 0] <- min(r[r > 0]) / 2; r <- r / sum(r); log(r) - mean(log(r)) }))
# Distance of each sample from the same participant's baseline sample.
distance_to_baseline <- function(m) {
  ref <- stats::setNames(meta$sample_id[meta$visit == "0"], meta$participant_id[meta$visit == "0"])[as.character(meta$participant_id)]
  vapply(seq_along(ref), function(i) if (is.na(ref[i])) NA_real_ else m[meta$sample_id[i], ref[[i]]], numeric(1))
}
distances <- list(`Bray-Curtis` = vegan::vegdist(rel, method = "bray"),
                  Aitchison = stats::dist(half_min_clr(rel)))
design <- data.frame(participant_id = factor(meta$participant_id), visit = factor(meta$visit, levels = VISITS))
within_participant <- permute::how(blocks = design$participant_id, nperm = N_PERMUTATIONS)

composition <- lapply(names(distances), function(name) {
  d <- distances[[name]]
  set.seed(SEED)
  permanova <- vegan::adonis2(d ~ participant_id + visit, data = design, permutations = within_participant, by = "terms")
  dispersion <- vegan::betadisper(d, design$visit, type = "median", bias.adjust = TRUE)
  set.seed(SEED)
  dispersion_test <- vegan::permutest(dispersion, permutations = within_participant)
  to_baseline <- data.frame(meta[, c("participant_id", "visit")], distance = distance_to_baseline(as.matrix(d)))
  pcoa <- stats::cmdscale(d, k = 2, eig = TRUE, add = TRUE)
  list(distance = name, permanova = as.data.frame(permanova), dispersion = as.data.frame(dispersion_test$tab),
       to_baseline = to_baseline,
       mean_to_baseline = stats::aggregate(distance ~ visit, to_baseline[to_baseline$visit != "0", ], mean),
       pcoa = data.frame(meta[, c("participant_id", "visit")], axis1 = pcoa$points[, 1], axis2 = pcoa$points[, 2]),
       pcoa_variance = pcoa$eig[1:2] / sum(pcoa$eig[pcoa$eig > 0]))
})
saveRDS(list(samples = meta, alpha_by_visit = alpha_by_visit, alpha_lmm = alpha_lmm, composition = stats::setNames(composition, names(distances))),
        result_path("diversity.rds"))

# ---------------- Part 4: NAvMix trajectory clustering (Figure 3B-D; Tables S5.1-S5.3) ----------------
for (layer in STOOL_LAYERS) {
  nav <- run_navmix(load_layer(layer), time_responsive_ids(layer), NAVMIX_K[[layer]], SEED + match(layer, LAYERS))
  saveRDS(nav, result_path("navmix", paste0(layer, ".rds")))
}

# ---------------- Part 5: features changed at all three post-baseline visits (Figure 3E) ----------------
saveRDS(rbind_list(lapply(STOOL_LAYERS, all_visit_summary)), result_path("all_visit_summary_stool.rds"))

# ---------------- Part 6: species-to-pathway mapping and fecal metabolite KEGG analysis (Figure 3F; Table S5.4) ----------------
# Metatranscriptomic species have no functional annotation and are not tested.
species_pathways <- read.csv(file.path(DATA_DIR, "microbiome", "species_kegg_pathways.csv"))   # feature_id, pathway_id, pathway
for (layer in c("fecal_metagenome", "fecal_metabolome")) {
  prep <- load_layer(layer)
  if (layer == "fecal_metabolome" && is.null(prep$annotation$kegg_id)) prep$annotation$kegg_id <- match_kegg_by_name(prep$annotation$feature_name)
  a <- read_result("navmix", paste0(layer, ".rds"))$assignments
  a <- a[a$cluster != "Noise", ]
  enrichment <- kegg_ora_layer(layer, split(a$feature_id, a$cluster), colnames(prep$z), prep$annotation, species_pathways)
  saveRDS(enrichment, result_path("enrichment", paste0("clusters_", layer, ".rds")))
}

# ---------------- Part 7: effect modification of trajectory-cluster scores, all six layers (Table S4.2) ----------------
medication <- read.csv(file.path(DATA_DIR, "medication.csv"))
cluster_modification <- rbind_list(lapply(LAYERS, function(layer) {
  scores <- read_result("navmix", paste0(layer, ".rds"))$cluster_scores
  m <- load_layer(layer)$metadata
  m <- m[match(scores$sample_id, m$sample_id), ]
  m$metformin <- medication$metformin[match(m$participant_id, medication$participant_id)]
  clusters <- setdiff(names(scores), c("sample_id", "participant_id", "visit"))
  rbind_list(lapply(MODIFIERS, function(modifier) rbind_list(lapply(clusters, function(cl)
    data.frame(layer = layer, modifier = modifier, cluster = cl, t(fit_time_modifier(scores[[cl]], m, modifier)))))))
}))
cluster_modification$fdr_interaction <- stats::ave(cluster_modification$p_interaction, cluster_modification$layer,
                                                  cluster_modification$modifier,
                                                  FUN = function(p) stats::p.adjust(p, "BH"))
saveRDS(cluster_modification, result_path("effect_modification", "clusters.rds"))
