# Figure 1: clinical response and time-dependent changes across the six omics layers; Figure S1: QC principal components.
# Panel A (study design schematic) is not drawn here.
source("plotting/plot_utils.R")
library(patchwork)

# ---------------- Panel B: clinical trajectories and outcomes ----------------
clinical <- read_result("clinical", "clinical_outcomes.rds")
traj <- rbind_list(lapply(c("weight", "hba1c"), function(v) {
  d <- stats::aggregate(clinical$cohort[[v]], list(visit = clinical$cohort$visit),
                        function(x) c(mean = mean(x, na.rm = TRUE), se = stats::sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))))
  data.frame(variable = v, week = as.numeric(d$visit), mean = d$x[, "mean"], se = d$x[, "se"])
}))
traj <- traj[is.finite(traj$mean), ]   # HbA1c was not measured at week 1
p_b1 <- ggplot(traj, aes(week, mean, colour = variable, fill = variable)) +
  geom_ribbon(aes(ymin = mean - se, ymax = mean + se), alpha = 0.18, colour = NA) +
  geom_line(linewidth = 0.6) + geom_point(size = 1.1) +
  facet_wrap(~variable, ncol = 1, scales = "free_y", labeller = as_labeller(c(weight = "Weight (kg)", hba1c = "HbA1c (%)"))) +
  scale_colour_manual(values = c(weight = "#3C4650", hba1c = "#8B5A3C"), guide = "none") +
  scale_fill_manual(values = c(weight = "#3C4650", hba1c = "#8B5A3C"), guide = "none") +
  week_axis() + labs(x = "Week", y = "Mean ± SE")
outcomes <- data.frame(visit = rep(paste0("W", clinical$outcomes$visit), 2),
                       outcome = rep(c("T2D remission", "Weight loss ≥10%"), each = nrow(clinical$outcomes)),
                       pct = c(clinical$outcomes$remission_pct, clinical$outcomes$weight_loss10_pct))
p_b2 <- ggplot(outcomes, aes(visit, pct, fill = outcome)) +
  geom_col(position = position_dodge(0.75), width = 0.68) +
  geom_text(aes(label = sprintf("%.1f%%", pct)), position = position_dodge(0.75), vjust = -0.3, size = 2) +
  scale_fill_manual(values = c("#3C4650", "#BFC3C9"), name = NULL) +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0), labels = function(x) paste0(x, "%")) + labs(x = NULL, y = "Participants")

# ---------------- Panel C: overall time effect of all tested features ----------------
tm <- rbind_list(lapply(LAYERS, function(layer) {
  d <- load_time(layer)
  d <- d[order(d$feature_id), ]
  d$x <- match(layer, LAYERS) - 1 + (seq_len(nrow(d)) - 0.5) / nrow(d)
  d[, c("layer", "feature_id", "feature_name", "x", "fdr_global", "p_global", "beta_w24", "time_responsive")]
}))
tm$log_fdr <- pmin(-log10(tm$fdr_global), 20)   # axis truncated at 20
tm$direction <- ifelse(tm$beta_w24 >= 0, "Increased at W24", "Decreased at W24")
top3 <- rbind_list(lapply(split(tm, tm$layer), function(d) { d <- d[labelable(d$layer, d$feature_id), ]; head(d[order(d$p_global), ], 3) }))
p_c <- ggplot(tm[tm$time_responsive, ], aes(x, log_fdr)) +
  geom_point(aes(colour = layer, shape = direction), size = 0.9, alpha = 0.6) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, colour = "grey50", linewidth = 0.25) +
  ggrepel::geom_text_repel(data = top3, aes(label = feature_name, colour = layer), size = 1.8, min.segment.length = 0, seed = SEED) +
  scale_colour_manual(values = LAYER_COLORS, guide = "none") +
  scale_shape_manual(values = c(`Increased at W24` = 17, `Decreased at W24` = 6), name = NULL) +
  scale_x_continuous(breaks = seq_along(LAYERS) - 0.5, labels = LAYER_LABELS[LAYERS]) +
  labs(x = NULL, y = "−log10(FDR-adjusted P)")

# ---------------- Panel D: visits at which time-responsive features changed ----------------
paths <- rbind_list(lapply(LAYERS, function(layer) {
  d <- load_time(layer)
  s <- ifelse(as.matrix(d[, c("selected_w1", "selected_w12", "selected_w24")]), "Changed", "Other")
  counts <- as.data.frame(table(W1 = s[, 1], W12 = s[, 2], W24 = s[, 3]))
  counts <- counts[counts$Freq > 0, ]
  data.frame(layer = layer, path = seq_len(nrow(counts)), proportion = counts$Freq / nrow(d), counts[, 1:3])
}))
lodes <- ggalluvial::to_lodes_form(paths, axes = c("W1", "W12", "W24"), key = "visit", value = "status", id = "flow")
lodes$fill <- ifelse(lodes$status == "Changed", lodes$layer, "Other")
lodes$layer <- factor(lodes$layer, LAYERS)
p_d <- ggplot(lodes, aes(x = visit, stratum = status, alluvium = flow, y = proportion, fill = fill)) +
  ggalluvial::geom_flow(width = 0.22, alpha = 0.4) + ggalluvial::geom_stratum(width = 0.22, colour = "white") +
  facet_wrap(~layer, ncol = 3, labeller = as_labeller(LAYER_LABELS)) +
  scale_fill_manual(values = c(LAYER_COLORS, Other = "grey80"), guide = "none") +
  labs(x = NULL, y = NULL) + theme(axis.text.y = element_blank(), axis.ticks = element_blank(), axis.line = element_blank())

# ---------------- Panel E: combinations of visits among features changed at one or more visits ----------------
combos <- rbind_list(lapply(LAYERS, function(layer) {
  d <- load_time(layer)
  data.frame(layer = layer, as.data.frame(table(combination = d$visit_combination[d$visit_combination != "None"])))
}))
combos$combination <- factor(combos$combination, c("W1", "W1+W12", "W1+W12+W24", "W1+W24", "W12", "W12+W24", "W24"))
p_e <- ggplot(combos, aes(Freq, factor(layer, rev(LAYERS)), fill = combination)) +
  geom_col(position = position_fill(reverse = TRUE), width = 0.65) +
  scale_fill_manual(values = c("#1D4E89", "#3F7CBF", "#7FAEDB", "#BFD7EE", "#C2571A", "#EE9547", "#F8CB96"), name = "Visits") +
  scale_x_continuous(labels = scales::percent) + scale_y_discrete(labels = LAYER_LABELS) +
  labs(x = "Share of features changed at ≥1 visit", y = NULL)

figure <- ((p_b1 + labs(title = "(B) Clinical trajectories and response")) | p_b2) /
  (p_c + labs(title = "(C) Global time-dependent molecular changes")) /
  ((p_d + labs(title = "(D) Time-specific changes")) | (p_e + labs(title = "(E) Visit combinations of changes"))) +
  plot_layout(heights = c(1.1, 0.8, 1))
save_figure(figure, "Figure1", 210, 260)

# ---------------- Figure S1: principal components of study and pooled QC samples ----------------
qc <- read_result("preprocessed", "qc_pca.rds")
qc$participant <- ifelse(qc$sample_type == "QC", "QC", qc$participant_id)
ids <- sort(unique(qc$participant[qc$participant != "QC"]))
s1_colours <- c(stats::setNames(grDevices::hcl.colors(length(ids), "Dynamic"), ids), QC = "black")
S1_LABELS <- c(serum_proteome = "Serum proteins", serum_metabolome = "Serum metabolites", fecal_metabolome = "Fecal metabolites")
s1_panel <- function(layer) {
  d <- qc[qc$layer == layer, ]
  ggplot(d, aes(pc1, pc2, colour = participant)) +
    stat_ellipse(data = d[d$participant != "QC", ], aes(group = participant), level = 0.68, linewidth = 0.25) +
    geom_point(size = 0.8, alpha = 0.85) +
    scale_colour_manual(values = s1_colours, name = NULL) +
    labs(title = S1_LABELS[[layer]], x = sprintf("PC1 (%.1f%%)", 100 * d$pc1_var[1]), y = sprintf("PC2 (%.1f%%)", 100 * d$pc2_var[1]))
}
figure_s1 <- wrap_plots(lapply(names(S1_LABELS), s1_panel), nrow = 1, guides = "collect") & theme(legend.position = "right")
save_figure(figure_s1, "FigureS1", 183, 70)
