# Figure 7: prediction of week-24 outcomes from baseline profiles.
# Panel A (cross-validation workflow schematic) is not drawn here.
source("plotting/plot_utils.R")
library(patchwork)
pred <- lapply(stats::setNames(LAYERS, LAYERS), function(layer) read_result("prediction", paste0(layer, ".rds")))
comparison <- read_result("prediction", "comparison.rds")

# ---------------- Panels B-C: week-24 BMI and HbA1c ----------------
# Left: observed vs participant-averaged out-of-fold prediction of the omics model.
# Right: Pearson R (95% CI) for both models.
continuous_panel <- function(outcome, colour, label) {
  d <- rbind_list(lapply(LAYERS, function(l) data.frame(layer = factor(l, LAYERS), pred[[l]][[outcome]]$averaged)))
  stats_r <- rbind_list(lapply(LAYERS, function(l) { a <- d[d$layer == l, ]; t <- stats::cor.test(a$observed, a$Omics)
    data.frame(layer = factor(l, LAYERS), label = sprintf("R = %.2f\nP = %.2g", t$estimate, t$p.value)) }))
  scatter <- ggplot(d, aes(observed, Omics)) +
    geom_abline(linetype = 2, colour = "grey60", linewidth = 0.3) + geom_point(colour = colour, size = 0.8, alpha = 0.75) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = colour, linewidth = 0.45) +
    geom_text(data = stats_r, aes(-Inf, Inf, label = label), hjust = -0.08, vjust = 1.15, size = 1.65) +
    facet_wrap(~layer, ncol = 3, labeller = as_labeller(LAYER_LABELS)) +
    labs(x = paste("Observed W24", label), y = "Predicted (mean out-of-fold)")
  r <- rbind_list(lapply(LAYERS, function(l) { s <- pred[[l]][[outcome]]$performance; data.frame(layer = factor(l, LAYERS), s[s$metric == "pearson_r", ]) }))
  forest <- ggplot(r, aes(estimate, factor(layer, rev(LAYERS)), colour = model)) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey60", linewidth = 0.3) +
    geom_errorbar(aes(xmin = low95, xmax = high95), width = 0, orientation = "y", position = position_dodge(0.6)) +
    geom_point(size = 1.5, position = position_dodge(0.6)) +
    scale_colour_manual(values = MODEL_COLORS, labels = MODEL_LABELS, name = NULL) + scale_y_discrete(labels = LAYER_LABELS) +
    labs(x = "Pearson R (95% CI)", y = NULL) + theme(legend.position = "top")
  scatter + forest + plot_layout(widths = c(1.5, 1))
}
p_b <- continuous_panel("bmi", TRAIT_COLORS[["bmi"]], "BMI")
p_c <- continuous_panel("hba1c", TRAIT_COLORS[["hba1c"]], "HbA1c")

# ---------------- Panel D: week-24 remission ----------------
roc <- rbind_list(lapply(LAYERS, function(l) {
  a <- pred[[l]]$remission$averaged
  rbind_list(lapply(c("Base", "Omics"), function(m) {
    r <- pROC::roc(a$observed, a[[m]], levels = c(0, 1), direction = "<", quiet = TRUE)
    data.frame(layer = factor(l, LAYERS), model = m, fpr = 1 - r$specificities, tpr = r$sensitivities)
  }))
}))
roc <- roc[order(roc$layer, roc$model, roc$fpr, roc$tpr), ]
band <- rbind_list(lapply(comparison$remission, `[[`, "roc_band"))
band$layer <- factor(band$layer, LAYERS)
auc_labels <- rbind_list(lapply(LAYERS, function(l) { x <- comparison$remission[[l]]; a <- x$auc
  data.frame(layer = factor(l, LAYERS), label = sprintf("Base = %.2f\nOmics = %.2f", a$auc[1], a$auc[2]), p = sprintf("P = %.3f", x$p_delong)) }))
p_d1 <- ggplot(roc, aes(fpr, tpr, colour = model)) +
  geom_abline(linetype = 2, colour = "grey75", linewidth = 0.3) +
  geom_ribbon(data = band, aes(x = fpr, ymin = se_low95, ymax = se_high95, fill = model), inherit.aes = FALSE, alpha = 0.15) +
  geom_step(linewidth = 0.55, direction = "vh") +
  geom_text(data = auc_labels, aes(1, 0, label = label), inherit.aes = FALSE, hjust = 1, vjust = 0, size = 1.55) +
  facet_wrap(~layer, nrow = 1, labeller = as_labeller(LAYER_LABELS)) +
  scale_colour_manual(values = MODEL_COLORS, labels = MODEL_LABELS, name = NULL) + scale_fill_manual(values = MODEL_COLORS, guide = "none") +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + labs(x = "1 − specificity", y = "Sensitivity") + theme(legend.position = "top")
auc_repeats <- rbind_list(lapply(LAYERS, function(l) { m <- pred[[l]]$remission$metrics_by_repeat
  data.frame(layer = factor(l, LAYERS), model = rep(colnames(m), each = dim(m)[3]), auc = c(t(m["auc", , ]))) }))
p_d2 <- ggplot(auc_repeats, aes(model, auc, colour = model)) +
  geom_hline(yintercept = 0.5, linetype = 2, colour = "grey75", linewidth = 0.3) +
  geom_boxplot(width = 0.55, outlier.size = 0.25, fill = "white") +
  geom_text(data = auc_labels, aes(1.5, Inf, label = p), inherit.aes = FALSE, vjust = 1.2, size = 1.6) +
  facet_wrap(~layer, nrow = 1) + scale_x_discrete(labels = c(Base = "Base", Omics = "Omics")) +
  scale_colour_manual(values = MODEL_COLORS, guide = "none") + labs(x = NULL, y = "AUC") + theme(strip.text = element_blank())

# ---------------- Panel E: most frequently selected omics features ----------------
# Microbial names longer than 30 characters are abbreviated (Candidatus -> Ca., first taxon name to its initial).
shorten_microbe <- function(x, layer, max_chars = 30) {
  if (!layer %in% c("fecal_metagenome", "fecal_metatranscriptome")) return(x)
  long <- nchar(x) > max_chars; x[long] <- sub("^Candidatus ", "Ca. ", x[long])
  long <- nchar(x) > max_chars & !grepl("virus", x)
  x[long] <- sub("^(Ca[.] )?[[]?([A-Z])[a-z]+[]]? ", "\\1\\2. ", x[long])
  x
}
selection <- rbind_list(lapply(LAYERS, function(l) rbind_list(lapply(names(pred[[l]]), function(o) {
  f <- pred[[l]][[o]]$selection_frequency
  f <- head(f[labelable(l, names(f))], 10)
  data.frame(layer = factor(l, LAYERS), outcome = o, rank = seq_along(f), feature = shorten_microbe(feature_label(l, names(f)), l),
             frequency = as.numeric(f))
}))))
p_e <- ggplot(selection, aes(factor(layer, LAYERS), 11 - rank)) +
  geom_tile(aes(fill = frequency), width = 0.95, height = 0.92) + geom_text(aes(label = feature), size = 1.7) +
  facet_grid(factor(outcome, c("bmi", "hba1c", "remission"), c("BMI", "HbA1c", "Remission")) ~ ., switch = "y") +
  scale_fill_gradient(low = "#FBEDEE", high = MODEL_COLORS[["Omics"]], limits = c(0, 1), name = "Selection frequency") +
  scale_x_discrete(labels = LAYER_LABELS) + labs(x = NULL, y = NULL) +
  theme(axis.text.y = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(), strip.placement = "outside")

figure <- wrap_elements(p_b + plot_annotation(title = "(B) Week-24 BMI prediction")) / wrap_elements(p_c + plot_annotation(title = "(C) Week-24 HbA1c prediction")) /
  wrap_elements((p_d1 / p_d2) + plot_annotation(title = "(D) Week-24 remission prediction")) / (p_e + labs(title = "(E) Top selected omics features")) +
  plot_layout(heights = c(1.15, 1.15, 1.15, 1.8))
save_figure(figure, "Figure7", 210, 297)
