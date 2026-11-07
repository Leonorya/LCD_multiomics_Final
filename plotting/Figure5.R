# Figure 5: longitudinal associations of omics features with adiposity and glycemic traits.
source("plotting/plot_utils.R")
library(patchwork)
res <- read_result("clinical_association", "model2.rds")
sig <- res$model2[is_sig(res$model2$fdr), ]
sig$node <- paste(sig$layer, sig$feature_id, sep = "::")
sig$label <- NA_character_
for (layer in unique(sig$layer)) sig$label[sig$layer == layer] <- feature_label(layer, sig$feature_id[sig$layer == layer])
sig$labelable <- labelable(sig$layer, sig$feature_id)

# ---------------- Panel A: circular map of feature-trait associations ----------------
# Outer ring: omics layers with feature labels; middle ring: adiposity traits; inner ring: glycemic traits.
# Shape and colour identify the trait; filled = positive, open = negative beta.
features <- unique(sig[, c("node", "layer", "label", "labelable")])
features <- features[order(match(features$layer, LAYERS), features$node), ]
features$x <- stats::ave(seq_len(nrow(features)), features$layer, FUN = seq_along) - 0.5
sig$x <- features$x[match(sig$node, features$node)]
sig$ring <- ifelse(sig$trait %in% ADIPOSITY_TRAITS, "Adiposity", "Glycemia")
sig$row <- c(bmi = 2.5, liver_fat = 1.5, pancreatic_fat = 0.5, hba1c = 2.5, fpg = 1.5, homa2_ir = 0.5)[sig$trait]
shape <- c(bmi = 21, liver_fat = 24, pancreatic_fat = 22, hba1c = 21, fpg = 24, homa2_ir = 22)
draw_circle <- function() {
  present <- LAYERS[LAYERS %in% features$layer]
  circlize::circos.clear()
  circlize::circos.par(start.degree = 90, gap.after = c(rep(2, length(present) - 1), 8), cell.padding = c(0, 0, 0, 0))
  circlize::circos.initialize(sectors = factor(features$layer, present), xlim = cbind(0, pmax(as.numeric(table(factor(features$layer, present))), 5)))
  circlize::circos.labels(sectors = features$layer, x = features$x, labels = ifelse(features$labelable, features$label, ""), side = "outside", cex = 0.4,
                          col = LAYER_COLORS[features$layer])
  circlize::circos.track(ylim = c(0, 1), track.height = 0.05, bg.border = NA, panel.fun = function(x, y)
    circlize::circos.rect(circlize::CELL_META$xlim[1], 0, circlize::CELL_META$xlim[2], 1, col = LAYER_COLORS[circlize::CELL_META$sector.index], border = NA))
  for (ring in c("Adiposity", "Glycemia")) {
    circlize::circos.track(ylim = c(0, 3), track.height = 0.12, bg.border = NA, bg.col = if (ring == "Adiposity") "#F2F7FA" else "#FBF3F1",
      panel.fun = function(x, y) {
        d <- sig[sig$layer == circlize::CELL_META$sector.index & sig$ring == ring, ]
        if (nrow(d)) circlize::circos.points(d$x, d$row, pch = shape[d$trait], cex = 0.6, col = TRAIT_COLORS[d$trait],
                                             bg = ifelse(d$beta > 0, TRAIT_COLORS[d$trait], "white"))
      })
  }
  circlize::circos.clear()
}
p_a <- ggplotify::as.ggplot(draw_circle) + theme(aspect.ratio = 1)

# ---------------- Panel B: significant associations per omics layer and trait ----------------
counts <- as.data.frame(table(layer = factor(sig$layer, LAYERS), trait = factor(sig$trait, CLINICAL_TRAITS)))
p_b <- ggplot(counts, aes(Freq, factor(layer, rev(LAYERS)), fill = trait)) + geom_col(width = 0.68) +
  scale_fill_manual(values = TRAIT_COLORS, labels = TRAIT_LABELS, name = NULL) + scale_y_discrete(labels = LAYER_LABELS) +
  labs(x = "Number of significant associations", y = NULL) + theme(legend.position = "right")

# ---------------- Panel C: strongest associations per trait (beta with 95% CI) ----------------
top <- rbind_list(lapply(split(sig[sig$labelable, ], factor(sig$trait[sig$labelable], CLINICAL_TRAITS)), function(d) head(d[order(d$p), ], 3)))
top$low95 <- top$beta - stats::qt(0.975, top$df) * top$se
top$high95 <- top$beta + stats::qt(0.975, top$df) * top$se
p_c <- ggplot(top, aes(beta, paste(trait, label, sep = "___"))) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey60") +
  geom_errorbar(aes(xmin = low95, xmax = high95), width = 0, orientation = "y", colour = "grey35") +
  geom_point(aes(fill = layer), shape = 21, size = 1.9) +
  facet_grid(factor(trait, CLINICAL_TRAITS) ~ ., scales = "free_y", space = "free_y", labeller = as_labeller(TRAIT_LABELS)) +
  scale_y_discrete(labels = function(x) sub("^.*___", "", x)) + scale_fill_manual(values = LAYER_COLORS, labels = LAYER_LABELS, name = NULL) +
  labs(x = "Standardized β (95% CI)", y = NULL) + theme(strip.text.y = element_text(angle = 0))

# ---------------- Panel D: KEGG pathways of adiposity- and glycemia-associated features ----------------
e <- res$set_enrichment
p_d <- plot_kegg(e, paste(e$set, SHORT_LAYER_LABELS[e$layer], sep = " · "))

# ---------------- Panels E-F: correlation networks (cross-layer edges) ----------------
plot_set_network <- function(set, flag) {
  edges <- res$set_networks[[set]]$edges
  edges <- edges[edges$cross_layer, ]
  nodes <- res$feature_sets[res$feature_sets[[flag]] & res$feature_sets$node %in% c(edges$node1, edges$node2), ]
  nodes$label <- features$label[match(nodes$node, features$node)]
  nodes$degree <- vapply(nodes$node, function(n) sum(edges$node1 == n | edges$node2 == n), integer(1))
  net <- layout_network(nodes, edges)
  labelled <- net$nodes[net$nodes$shared | net$nodes$node %in% head(net$nodes$node[order(-net$nodes$degree)], 6), ]
  labelled <- labelled[labelable(labelled$layer, labelled$feature_id), ]
  ggplot() +
    geom_segment(data = net$edges, aes(x, y, xend = xend, yend = yend, colour = r), alpha = 0.25, linewidth = 0.22) +
    geom_point(data = net$nodes, aes(x, y, fill = layer), shape = 21, size = 1.8, colour = "white") +
    geom_point(data = net$nodes[net$nodes$shared, ], aes(x, y), shape = 21, size = 2.2, stroke = 0.6) +   # shared features
    ggrepel::geom_text_repel(data = labelled, aes(x, y, label = label), size = 1.8, seed = SEED) +
    scale_colour_gradient2(low = DIVERGING[["low"]], mid = "grey85", high = DIVERGING[["high"]], limits = c(-1, 1), name = "Partial r") +
    scale_fill_manual(values = LAYER_COLORS, guide = "none") + coord_equal() + theme_void(base_size = 7) + theme(legend.position = "bottom")
}
p_e <- plot_set_network("Adiposity", "adiposity")
p_f <- plot_set_network("Glycemia", "glycemia")

figure <- ((p_a + labs(title = "(A) Feature–trait associations")) / (p_b + labs(title = "(B) Associations per omics layer")) + plot_layout(heights = c(3.3, 1))) |
  (((p_c + labs(title = "(C) Top associations")) | (p_d + labs(title = "(D) KEGG pathways"))) /
     ((p_e + labs(title = "(E) Adiposity network")) | (p_f + labs(title = "(F) Glycemia network"))))
save_figure(figure, "Figure5", 297, 210)
