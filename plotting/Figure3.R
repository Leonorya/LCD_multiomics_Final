# Figure 3: gut microbial diversity, trajectories and pathway mapping of the fecal omics layers.
source("plotting/plot_utils.R")
library(patchwork)

# ---------------- Panel A: diversity, gene count and community composition ----------------
div <- read_result("diversity.rds")
alpha <- rbind(data.frame(div$samples[, c("participant_id", "visit")], metric = "Shannon index", value = div$samples$shannon),
               data.frame(div$samples[, c("participant_id", "visit")], metric = "Gene count", value = 10^div$samples$log10_gene_count))
alpha_p <- data.frame(metric = c("Shannon index", "Gene count"),
                      label = sprintf("Time P = %.3g", div$alpha_lmm$p_global[match(c("shannon", "log10_gene_count"), div$alpha_lmm$metric)]))
p_a1 <- ggplot(alpha, aes(factor(visit, VISITS), value, fill = factor(visit, VISITS))) +
  geom_boxplot(width = 0.55, outlier.shape = NA, linewidth = 0.3) +
  geom_line(aes(group = participant_id), colour = "grey85", linewidth = 0.15) + geom_point(size = 0.45, alpha = 0.55) +
  geom_text(data = alpha_p, aes(x = 0.6, y = Inf, label = label), inherit.aes = FALSE, hjust = 0, vjust = 1.3, size = 1.8) +
  facet_wrap(~metric, scales = "free_y") + scale_fill_manual(values = VISIT_COLORS, guide = "none") +
  scale_x_discrete(labels = paste0("W", VISITS)) + labs(x = NULL, y = NULL)

bray <- div$composition$`Bray-Curtis`
centroids <- stats::aggregate(cbind(axis1, axis2) ~ visit, bray$pcoa, mean)
permanova_time <- bray$permanova["visit", ]
p_a2 <- ggplot(bray$pcoa, aes(axis1, axis2, colour = factor(visit, VISITS))) +
  stat_ellipse(level = 0.68, linewidth = 0.35) + geom_point(size = 0.7) +
  geom_path(data = centroids[match(VISITS, centroids$visit), ], aes(group = 1), colour = "grey20", arrow = grid::arrow(length = grid::unit(1, "mm"))) +
  scale_colour_manual(values = VISIT_COLORS, labels = paste0("W", VISITS), name = NULL) +
  labs(x = sprintf("PCoA1 (%.1f%%)", 100 * bray$pcoa_variance[1]), y = sprintf("PCoA2 (%.1f%%)", 100 * bray$pcoa_variance[2]),
       subtitle = sprintf("PERMANOVA R² = %.1f%%, P = %.3f", 100 * permanova_time$R2, permanova_time$`Pr(>F)`))

to_base <- bray$to_baseline[bray$to_baseline$visit != "0", ]
p_a3 <- ggplot(to_base, aes(factor(visit, VISITS[-1]), distance, fill = factor(visit, VISITS[-1]))) +
  geom_boxplot(width = 0.55, outlier.shape = NA, linewidth = 0.3) +
  geom_line(aes(group = participant_id), colour = "grey85", linewidth = 0.15) + geom_point(size = 0.45, alpha = 0.55) +
  scale_fill_manual(values = VISIT_COLORS, guide = "none") + scale_x_discrete(labels = paste0("W", VISITS[-1])) +
  labs(x = NULL, y = "Bray–Curtis distance to W0")

# ---------------- Panels B-D: trajectory clusters ----------------
p_b <- plot_trajectories("fecal_metagenome", ncol = 4)
p_c <- plot_trajectories("fecal_metatranscriptome", ncol = 2)
p_d <- plot_trajectories("fecal_metabolome", ncol = 3)

# ---------------- Panel E: features changed at weeks 1, 12 and 24 ----------------
p_e <- plot_all_visit_heatmap(STOOL_LAYERS)

# ---------------- Panel F: species-to-pathway mapping and fecal metabolite KEGG pathways ----------------
enrichment <- rbind_list(lapply(c("fecal_metagenome", "fecal_metabolome"), function(layer) read_result("enrichment", paste0("clusters_", layer, ".rds"))))
p_f <- plot_kegg(enrichment, paste(SHORT_LAYER_LABELS[enrichment$layer], sub("^C", "Cluster ", enrichment$set), sep = " · "))

figure <- (p_b + labs(title = "(B) Metagenome") | p_c + labs(title = "(C) Metatranscriptome") | p_d + labs(title = "(D) Fecal metabolites")) /
  ((p_a1 + labs(title = "(A) Diversity, gene count, and community composition")) / (p_a2 | p_a3) |
     p_e + labs(title = "(E) Changed at W1, W12, and W24") | p_f + labs(title = "(F) KEGG pathways of trajectory clusters")) +
  plot_layout(heights = c(0.9, 1.1))
save_figure(figure, "Figure3", 297, 210)
