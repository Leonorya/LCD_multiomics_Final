# Shared settings and plotting helpers for the figure scripts.
# Run from the repository root after the analysis scripts, e.g. `Rscript plotting/Figure1.R`.
# The scripts draw the content of each panel; the final panel arrangement and styling were edited manually.
source("analysis/utils.R")
library(ggplot2)

# ---------------- Part 1: settings ----------------
FIGURE_DIR <- "figures"
LAYER_LABELS <- c(pbmc_transcriptome = "PBMC transcriptomics", serum_proteome = "Serum proteomics",
                  serum_metabolome = "Serum metabolomics", fecal_metagenome = "Metagenomics",
                  fecal_metatranscriptome = "Metatranscriptomics", fecal_metabolome = "Fecal metabolomics")
SHORT_LAYER_LABELS <- c(pbmc_transcriptome = "PBMC", serum_proteome = "Serum protein", serum_metabolome = "Serum metabolite",
                        fecal_metagenome = "Metagenome", fecal_metatranscriptome = "Metatranscriptome", fecal_metabolome = "Fecal metabolite")
LAYER_COLORS <- c(pbmc_transcriptome = "#4C78A8", serum_proteome = "#E45756", serum_metabolome = "#72B7B2",
                  fecal_metagenome = "#54A24B", fecal_metatranscriptome = "#B279A2", fecal_metabolome = "#E6A93D")
TRAIT_LABELS <- c(bmi = "BMI", liver_fat = "Liver fat", pancreatic_fat = "Pancreatic fat",
                  hba1c = "HbA1c", fpg = "FPG", homa2_ir = "HOMA2-IR")
TRAIT_COLORS <- c(bmi = "#1F6E8C", liver_fat = "#4A9BB8", pancreatic_fat = "#9CCBDB",
                  hba1c = "#B03A2E", fpg = "#D9735E", homa2_ir = "#EFB3A3")
VISIT_COLORS  <- c(`0` = "#8C8C8C", `1` = "#E3A33B", `12` = "#3FA796", `24` = "#3B6FB6")
MODEL_COLORS  <- c(Base = "#7FA6C9", Omics = "#D1495B")
MODEL_LABELS  <- c(Base = "Base model", Omics = "Omics model")
STATUS_COLORS <- c(Remission = "#2A9D8F", `No remission` = "#E76F51")
TOP5_COLORS   <- c("#D7301F", "#2C7FB8", "#31A354", "#F0A02C", "#756BB1")
DIVERGING     <- c(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B")
KEGG_COLORS   <- c(Metabolism = "#E07A9A", `Genetic Information Processing` = "#4DB89A", `Environmental Information Processing` = "#E8A04D",
                   `Cellular Processes` = "#C9A93A", `Organismal Systems` = "#A88BD0", `Human Diseases` = "#3D73B0")
SHORT_PATHWAY_NAMES <- c(`Biosynthesis of siderophore group nonribosomal peptides` = "Biosynthesis of siderophore NRP",
                         `Glucocorticoid and mineralocorticoid receptor agonists/antagonists` = "Corticosteroid receptor agonists/antagonists")

# ---------------- Part 2: theme and output ----------------
theme_paper <- function(base_size = 7) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(axis.line = element_line(linewidth = 0.35), axis.ticks = element_line(linewidth = 0.35),
          strip.background = element_blank(), strip.text = element_text(face = "bold"),
          plot.title = element_text(size = base_size + 0.5, face = "bold"), legend.position = "bottom")
}
theme_set(theme_paper())
save_figure <- function(plot, name, width_mm, height_mm) {
  dir.create(FIGURE_DIR, showWarnings = FALSE)
  ggsave(file.path(FIGURE_DIR, paste0(name, ".pdf")), plot, width = width_mm, height = height_mm, units = "mm", device = grDevices::cairo_pdf)
}
week_axis <- function() scale_x_continuous(breaks = c(0, 1, 12, 24), labels = c("0", "1", "12", "24"))
feature_label <- function(layer, ids) {
  a <- load_layer(layer)$annotation
  ifelse(is.na(a$feature_name[match(ids, a$feature_id)]), ids, a$feature_name[match(ids, a$feature_id)])
}
# Features putatively annotated as exogenous compounds are not labeled in the figures.
LABEL_EXCLUSIONS <- read.csv("annotation/label_exclusions.csv")
labelable <- function(layer, ids) !paste(layer, ids) %in% paste(LABEL_EXCLUSIONS$layer, LABEL_EXCLUSIONS$feature_id)

# ---------------- Part 3: shared panels ----------------
# Mean baseline-standardized trajectory of every feature in each cluster; the five labelable features with the smallest
# P for the overall time effect are coloured. The noise group is shown last.
plot_trajectories <- function(layer, ncol) {
  nav <- read_result("navmix", paste0(layer, ".rds"))
  tm <- load_time(layer)
  a <- nav$assignments
  a$p_global <- tm$p_global[match(a$feature_id, tm$feature_id)]
  a <- a[order(a$p_global), ]
  a$label <- feature_label(layer, a$feature_id)
  eligible <- labelable(layer, a$feature_id)
  eligible[eligible] <- !duplicated(paste(a$cluster, a$label)[eligible])
  a$rank <- NA_integer_
  a$rank[eligible] <- stats::ave(seq_len(sum(eligible)), a$cluster[eligible], FUN = seq_along)
  d <- data.frame(feature_id = rep(rownames(nav$trajectories), 4), week = rep(as.numeric(VISITS), each = nrow(nav$trajectories)),
                  mean_z = as.vector(nav$trajectories))
  d <- merge(d, a, by = "feature_id")
  n <- table(a$cluster)
  d$cluster <- factor(d$cluster, c(sort(setdiff(names(n), "Noise")), "Noise"))
  labels <- stats::setNames(paste0(sub("^C", "Cluster ", names(n)), " (n = ", n, ")"), names(n))
  top <- d[!is.na(d$rank) & d$rank <= 5, ]
  ggplot(d, aes(week, mean_z, group = feature_id)) +
    geom_hline(yintercept = 0, linetype = 2, linewidth = 0.3) +
    geom_line(colour = "grey80", linewidth = 0.2, alpha = 0.4) +
    geom_line(data = top, aes(colour = factor(rank)), linewidth = 0.5) +
    geom_text(data = top[top$week == 24, ], aes(x = 25, label = label, colour = factor(rank)), hjust = 0, size = 1.6) +
    facet_wrap(~cluster, ncol = ncol, scales = "free_y", labeller = as_labeller(labels)) +
    scale_colour_manual(values = TOP5_COLORS, guide = "none") + week_axis() +
    labs(x = "Time (week)", y = "Z-score", title = LAYER_LABELS[[layer]]) +
    theme(panel.border = element_rect(fill = NA, linewidth = 0.4))
}

# Week-1, week-12 and week-24 effects of features changed at all three post-baseline visits, ordered by the week-24 effect.
plot_all_visit_heatmap <- function(layers) {
  d <- rbind_list(lapply(layers, function(layer) {
    tm <- load_time(layer)
    tm <- tm[tm$selected_w1 & tm$selected_w12 & tm$selected_w24, ]
    tm$order <- rank(-tm$beta_w24, ties.method = "first")
    data.frame(layer = factor(LAYER_LABELS[[layer]], LAYER_LABELS[layers]), order = rep(tm$order, 3), visit = rep(c("W1", "W12", "W24"), each = nrow(tm)),
               beta = c(tm$beta_w1, tm$beta_w12, tm$beta_w24))
  }))
  limit <- max(abs(d$beta))
  ggplot(d, aes(visit, -order, fill = beta)) + geom_raster() +
    facet_wrap(~layer, nrow = 1, scales = "free_y") +
    scale_fill_gradient2(low = DIVERGING[["low"]], mid = DIVERGING[["mid"]], high = DIVERGING[["high"]], limits = c(-limit, limit), name = "β vs W0") +
    labs(x = NULL, y = NULL) + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line = element_blank())
}

# Top KEGG pathways per feature set: x = gene ratio, size = -log10(P), colour = KEGG category.
plot_kegg <- function(enrichment, group, n_top = 3) {
  enrichment$group <- group
  top <- rbind_list(lapply(split(enrichment, enrichment$group), function(d) head(d[order(d$p), ], n_top)))
  top$pathway <- ifelse(top$pathway %in% names(SHORT_PATHWAY_NAMES), SHORT_PATHWAY_NAMES[top$pathway], top$pathway)
  top$ratio <- vapply(strsplit(top$gene_ratio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]), numeric(1))
  top$pathway <- stats::reorder(paste(top$group, top$pathway, sep = "___"), top$ratio)
  ggplot(top, aes(ratio, pathway)) +
    geom_segment(aes(x = 0, xend = ratio, yend = pathway), colour = "grey60", linewidth = 0.25) +
    geom_point(aes(size = -log10(p), colour = pathway_class)) +
    facet_grid(group ~ ., scales = "free_y", space = "free_y") +
    scale_y_discrete(labels = function(x) sub("^.*___", "", x)) +
    scale_colour_manual(values = KEGG_COLORS, name = "KEGG category") +
    labs(x = "Gene ratio", y = NULL, size = "−log10(P)") +
    theme(strip.text.y = element_text(angle = 0, hjust = 0))
}

# Force-directed network layout; nodes and edges are data frames with node / node1, node2 columns.
layout_network <- function(nodes, edges, seed = SEED) {
  g <- igraph::graph_from_data_frame(data.frame(from = edges$node1, to = edges$node2, weight = abs(edges$r)), directed = FALSE, vertices = nodes$node)
  set.seed(seed)
  xy <- igraph::layout_with_fr(g, weights = igraph::E(g)$weight)
  nodes$x <- xy[, 1]; nodes$y <- xy[, 2]
  edges$x <- nodes$x[match(edges$node1, nodes$node)]; edges$y <- nodes$y[match(edges$node1, nodes$node)]
  edges$xend <- nodes$x[match(edges$node2, nodes$node)]; edges$yend <- nodes$y[match(edges$node2, nodes$node)]
  list(nodes = nodes, edges = edges)
}
