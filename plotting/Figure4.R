# Figure 4: cross-omics network of the most time-responsive features.
source("plotting/plot_utils.R")
library(patchwork)
net <- read_result("network", "cross_omics.rds")
nodes <- net$features[!is.na(net$features$module), ]
nodes$labelable <- labelable(nodes$layer, nodes$feature_id)
modules <- paste0("Module", seq_along(unique(nodes$module)))
nodes$module <- factor(nodes$module, modules)
edges <- net$edges[net$edges$node1 %in% nodes$node & net$edges$node2 %in% nodes$node, ]

# ---------------- Panel A: module map (each module laid out separately and placed on a circle) ----------------
angle <- pi / 2 - 2 * pi * (seq_along(modules) - 1) / length(modules)
radius <- sqrt(as.numeric(table(nodes$module))) / max(sqrt(as.numeric(table(nodes$module)))) * 0.42
placed <- rbind_list(lapply(seq_along(modules), function(i) {
  n <- nodes[nodes$module == modules[i], ]
  l <- layout_network(n, edges[edges$node1 %in% n$node & edges$node2 %in% n$node, ])$nodes
  xy <- scale(cbind(l$x, l$y), scale = FALSE); xy <- xy / max(sqrt(rowSums(xy^2)), 1e-9) * radius[i]
  l$x <- xy[, 1] + cos(angle[i]); l$y <- xy[, 2] + sin(angle[i])
  l
}))
hulls <- rbind_list(lapply(split(placed, placed$module, drop = TRUE), function(d) d[grDevices::chull(d$x, d$y), c("x", "y", "module")]))
hubs <- rbind_list(lapply(split(placed[placed$labelable, ], placed$module[placed$labelable], drop = TRUE), function(d) head(d[order(-d$cross_layer_strength), ], 3)))
cross <- edges[edges$cross_layer, ]
cross$x <- placed$x[match(cross$node1, placed$node)]; cross$y <- placed$y[match(cross$node1, placed$node)]
cross$xend <- placed$x[match(cross$node2, placed$node)]; cross$yend <- placed$y[match(cross$node2, placed$node)]
p_a <- ggplot() +
  geom_polygon(data = hulls, aes(x, y, group = module), fill = "grey96", colour = "grey80", linewidth = 0.3) +
  geom_segment(data = cross, aes(x, y, xend = xend, yend = yend), colour = "grey45", alpha = 0.08, linewidth = 0.18) +
  geom_point(data = placed, aes(x, y, fill = layer, size = cross_layer_strength), shape = 21, colour = "white", stroke = 0.2) +
  ggrepel::geom_text_repel(data = hubs, aes(x, y, label = feature_name), size = 1.8, min.segment.length = 0, seed = SEED) +
  scale_fill_manual(values = LAYER_COLORS, guide = "none") + scale_size_continuous(range = c(0.6, 3), name = "Cross-layer connectivity (Σ|r|)") +
  coord_equal() + theme_void(base_size = 7) + theme(legend.position = "bottom")

# ---------------- Panel B: module composition ----------------
p_b <- ggplot(nodes, aes(sub("Module", "", module), fill = layer)) + geom_bar(width = 0.7) +
  scale_fill_manual(values = LAYER_COLORS, labels = LAYER_LABELS, name = NULL) + labs(x = "Module", y = "Features") +
  theme(legend.position = "right")

# ---------------- Panel C: hub features per module ----------------
ranked <- rbind_list(lapply(split(nodes[nodes$labelable, ], nodes$module[nodes$labelable], drop = TRUE), function(d) head(d[order(-d$cross_layer_strength), ], 5)))
ranked$label <- paste(ranked$feature_name, ifelse(ranked$beta_w24 > 0, "↑", "↓"))
p_c <- ggplot(ranked, aes(cross_layer_strength, stats::reorder(paste(module, label, sep = "___"), cross_layer_strength))) +
  geom_segment(aes(x = 0, xend = cross_layer_strength, yend = stats::reorder(paste(module, label, sep = "___"), cross_layer_strength)), colour = "grey60") +
  geom_point(aes(colour = layer), size = 1.8) + facet_grid(module ~ ., scales = "free_y", space = "free_y") +
  scale_y_discrete(labels = function(x) sub("^.*___", "", x)) + scale_colour_manual(values = LAYER_COLORS, guide = "none") +
  labs(x = "Cross-layer connectivity (Σ|r|)", y = NULL) + theme(strip.text.y = element_text(angle = 0))

# ---------------- Panel D: strongest cross-layer partners of each module hub ----------------
stars <- rbind_list(lapply(split(nodes, nodes$module, drop = TRUE), function(d) {
  d <- d[d$labelable, ]
  hub <- d[which.max(d$cross_layer_strength), ]
  e <- net$edges[net$edges$cross_layer & (net$edges$node1 == hub$node | net$edges$node2 == hub$node), ]
  other <- ifelse(e$node1 == hub$node, e$node2, e$node1)
  e <- e[labelable(net$features$layer[match(other, net$features$node)], net$features$feature_id[match(other, net$features$node)]), ]
  e <- head(e[order(-abs(e$r)), ], 5)
  partner <- ifelse(e$node1 == hub$node, e$node2, e$node1)
  a <- pi / 2 - 2 * pi * (seq_along(partner) - 1) / max(length(partner), 1)
  data.frame(hub = paste0(hub$module, ": ", hub$feature_name), hub_layer = hub$layer, r = e$r, x = cos(a), y = sin(a),
             partner = net$features$feature_name[match(partner, net$features$node)], partner_layer = net$features$layer[match(partner, net$features$node)])
}))
p_d <- ggplot(stars) +
  geom_segment(aes(0, 0, xend = x, yend = y, linetype = r > 0), colour = "grey45", linewidth = 0.45) +
  geom_point(aes(x, y, fill = partner_layer), shape = 21, size = 1.9, colour = "white") +
  geom_point(data = unique(stars[, c("hub", "hub_layer")]), aes(0, 0, fill = hub_layer), shape = 21, size = 3) +
  geom_text(aes(1.25 * x, 1.25 * y, label = partner), size = 1.5) +
  facet_wrap(~hub, ncol = 2) +
  scale_linetype_manual(values = c(`TRUE` = "solid", `FALSE` = "22"), labels = c(`TRUE` = "Positive", `FALSE` = "Negative"), name = "Partial r") +
  scale_fill_manual(values = LAYER_COLORS, guide = "none") + coord_equal(xlim = c(-2.4, 2.4), ylim = c(-1.45, 1.45)) +
  theme_void(base_size = 7) + theme(legend.position = "bottom")

# ---------------- Panel E: KEGG pathways of cross-layer modules ----------------
e <- net$module_enrichment
p_e <- plot_kegg(e, paste(sub("Module", "Module ", e$set), SHORT_LAYER_LABELS[e$layer], sep = " · "))

figure <- ((p_a + labs(title = "(A) Cross-omics network")) / (p_b + labs(title = "(B) Module composition"))) |
  ((p_c + labs(title = "(C) Hub features per module")) / (p_d + labs(title = "(D) Strongest cross-layer partners of each module hub"))) |
  (p_e + labs(title = "(E) KEGG pathways of cross-layer modules"))
save_figure(figure, "Figure4", 297, 210)
