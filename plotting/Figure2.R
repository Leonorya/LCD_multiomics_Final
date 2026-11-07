# Figure 2: trajectories, pathways and protein interactions of the blood omics layers.
source("plotting/plot_utils.R")
library(patchwork)

# ---------------- Panels A-C: trajectory clusters ----------------
p_a <- plot_trajectories("pbmc_transcriptome", ncol = 3)
p_b <- plot_trajectories("serum_proteome", ncol = 2)
p_c <- plot_trajectories("serum_metabolome", ncol = 3)

# ---------------- Panel D: features changed at weeks 1, 12 and 24 ----------------
p_d <- plot_all_visit_heatmap(BLOOD_LAYERS)

# ---------------- Panel E: KEGG pathways of trajectory clusters ----------------
enrichment <- rbind_list(lapply(BLOOD_LAYERS, function(layer) read_result("enrichment", paste0("clusters_", layer, ".rds"))))
p_e <- plot_kegg(enrichment, paste(SHORT_LAYER_LABELS[enrichment$layer], sub("^C", "Cluster ", enrichment$set), sep = " · "))

# ---------------- Panel F: STRING network of time-responsive serum proteins ----------------
ppi <- read_result("ppi", "serum_proteome.rds")
nodes <- ppi$nodes[ppi$nodes$degree > 0, ]
names(nodes)[names(nodes) == "string_id"] <- "node"
edges <- data.frame(node1 = ppi$edges$from, node2 = ppi$edges$to, r = 1)
net <- layout_network(nodes, edges)
limit <- max(abs(nodes$beta_w24))
p_f <- ggplot() +
  geom_segment(data = net$edges, aes(x, y, xend = xend, yend = yend), colour = "grey45", linewidth = 0.3, linetype = 2) +
  geom_point(data = net$nodes[!net$nodes$highest_core, ], aes(x, y, fill = beta_w24, size = degree), shape = 21, stroke = 0.25) +
  geom_point(data = net$nodes[net$nodes$highest_core, ], aes(x, y, fill = beta_w24, size = degree), shape = 21, stroke = 0.9) +   # highest k-core
  ggrepel::geom_text_repel(data = net$nodes, aes(x, y, label = gene_symbol), size = 1.8, seed = SEED) +
  scale_fill_gradient2(low = DIVERGING[["low"]], mid = DIVERGING[["mid"]], high = DIVERGING[["high"]], limits = c(-limit, limit), name = "W24 β") +
  scale_size_continuous(range = c(1.4, 4), name = "Degree") + coord_equal() + theme_void(base_size = 7) + theme(legend.position = "bottom")

figure <- (p_a + labs(title = "(A) PBMC transcripts") | p_b + labs(title = "(B) Serum proteins") | p_c + labs(title = "(C) Serum metabolites")) /
  (p_d + labs(title = "(D) Changed at W1, W12, and W24") | p_e + labs(title = "(E) KEGG pathways of trajectory clusters") |
     p_f + labs(title = "(F) Protein–protein interaction network")) +
  plot_layout(heights = c(0.9, 1.1))
save_figure(figure, "Figure2", 297, 210)
