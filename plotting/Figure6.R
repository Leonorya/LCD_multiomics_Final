# Figure 6: molecular features associated with T2D remission.
source("plotting/plot_utils.R")
library(patchwork)
res <- read_result("remission", "remission.rds")

# ---------------- Panel A: remission associations per layer ----------------
r <- res$remission
r$log_fdr <- pmin(-log10(r$fdr), 20)
r$layer <- factor(r$layer, LAYERS)
r$class <- ifelse(is_sig(r$fdr), ifelse(r$beta > 0, "Higher odds of remission (FDR-adjusted P < 0.05)", "Lower odds of remission (FDR-adjusted P < 0.05)"), "Not significant")
r$label <- NA_character_
for (layer in unique(r$layer)) r$label[r$layer == layer] <- feature_label(as.character(layer), r$feature_id[r$layer == layer])
candidates <- r[r$class != "Not significant" & labelable(as.character(r$layer), r$feature_id), ]
labels <- rbind_list(lapply(split(candidates, candidates$layer, drop = TRUE),
                            function(d) head(d[order(-(abs(d$std_beta) + d$log_fdr / 20)), ], 6)))
p_a <- ggplot(r, aes(std_beta, log_fdr)) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, colour = "grey40", linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey40", linewidth = 0.3) +
  geom_point(data = r[r$class == "Not significant", ], colour = "grey80", size = 0.9) +
  geom_point(data = r[r$class != "Not significant", ], aes(fill = class), shape = 21, colour = "white", size = 1.5, stroke = 0.15) +
  ggrepel::geom_label_repel(data = labels, aes(label = label), size = 1.6, label.size = 0.15, min.segment.length = 0, seed = SEED) +
  facet_wrap(~layer, nrow = 2, labeller = as_labeller(LAYER_LABELS)) +
  scale_fill_manual(values = c(`Higher odds of remission (FDR-adjusted P < 0.05)` = DIVERGING[["high"]], `Lower odds of remission (FDR-adjusted P < 0.05)` = DIVERGING[["low"]]), name = NULL) +
  labs(x = "Standardized β", y = "−log10(FDR-adjusted P)") +
  theme(panel.border = element_rect(fill = NA, linewidth = 0.4))

# ---------------- Panel B: mean gut species composition by week-24 remission status ----------------
comp <- res$composition
species_cols <- setdiff(names(comp), c("remission_w24", "visit"))
long <- data.frame(status = rep(ifelse(comp$remission_w24 == 1, "Remission", "No remission"), length(species_cols)),
                   visit = rep(comp$visit, length(species_cols)), species = factor(rep(species_cols, each = nrow(comp)), rev(species_cols)),
                   abundance = unlist(comp[, species_cols]))
p_b <- ggplot(long, aes(factor(visit, VISITS), abundance, fill = species)) + geom_col(width = 0.42) +
  facet_wrap(~factor(status, names(STATUS_COLORS)), ncol = 1) +
  scale_x_discrete(labels = paste0("W", VISITS)) + scale_y_continuous(expand = c(0, 0)) +
  labs(x = "Time (week)", y = "Relative abundance", fill = NULL) +
  theme(legend.position = "right", legend.text = element_text(face = "italic"))

# ---------------- Panel C: species differing at week 24 ----------------
top <- head(res$species_w24[order(res$species_w24$p), ], 10)
top$p_label <- ifelse(top$p < 0.001, "P < 0.001", sprintf("P = %.3f", top$p))
p_c <- ggplot(top, aes(log2_fold_change, stats::reorder(species, -p), fill = log2_fold_change > 0)) + geom_col(width = 0.65) +
  geom_vline(xintercept = 0, linewidth = 0.3) +
  geom_text(aes(x = ifelse(log2_fold_change > 0, -0.15, 0.15), label = p_label, hjust = ifelse(log2_fold_change > 0, 1, 0)), size = 1.6) +
  scale_fill_manual(values = c(`TRUE` = STATUS_COLORS[["Remission"]], `FALSE` = STATUS_COLORS[["No remission"]]),
                    labels = c(`TRUE` = "Higher in remission", `FALSE` = "Higher in no remission"), name = NULL) +
  labs(x = "log2 fold change of mean abundance at W24", y = NULL) + theme(axis.text.y = element_text(face = "italic"))

# ---------------- Panel D: agreement with the longitudinal HbA1c association ----------------
con <- res$concordance
con$std_beta <- r$std_beta[match(paste(con$layer, con$feature_id), paste(r$layer, r$feature_id))]
p_d <- ggplot(con, aes(std_beta, beta_hba1c, fill = layer)) +
  annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = 0, fill = "grey95") +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = 0, ymax = Inf, fill = "grey95") +   # concordant quadrants
  geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.25) + geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.25) +
  geom_point(shape = 21, colour = "white", size = 1.1, stroke = 0.1, alpha = 0.8) +
  annotate("text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.4, size = 1.9,
           label = sprintf("Concordant (shaded): %.0f%% of %d", 100 * mean(con$concordant, na.rm = TRUE), nrow(con))) +
  scale_fill_manual(values = LAYER_COLORS, labels = LAYER_LABELS, name = NULL) +
  labs(x = "Remission standardized β", y = "Longitudinal HbA1c association β")

figure <- (p_a + labs(title = "(A) Multi-omic features associated with T2D remission")) /
  ((p_b + labs(title = "(B) Gut microbiome composition by remission status")) |
     ((p_c + labs(title = "(C) Species differing at W24")) / (p_d + labs(title = "(D) Agreement with the continuous HbA1c association")))) +
  plot_layout(heights = c(0.8, 1.3))
save_figure(figure, "Figure6", 210, 260)
