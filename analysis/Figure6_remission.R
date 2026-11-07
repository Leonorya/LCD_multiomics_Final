# Figure 6 and Tables S11: molecular features associated with T2D remission.
# Input : results/preprocessed/<layer>.rds, results/time_model/<layer>.rds, results/clinical_association/model2.rds,
#         data/microbiome/species_relative_abundance.csv
# Output: results/remission/*.rds
source("analysis/utils.R")

# ---------------- Part 1: settings ----------------
REMISSION_VISITS  <- c("0", "12", "24")   # no participant was in remission at baseline
MIN_MEAN_ABUNDANCE <- 0.001              # species compared at week 24 (mean relative abundance >= 0.1%)
TOP_SPECIES        <- 15                 # species shown in the composition plot

# ---------------- Part 2: remission GLMM (Figure 6A; Tables S11.1-S11.2) ----------------
# remission ~ omics_z + visit + age + sex + baseline BMI + baseline HbA1c + (1 | participant), binomial.
fit_remission_glmm <- function(x, d) {
  d$x <- x
  d <- d[stats::complete.cases(d), ]
  fit <- tryCatch(suppressWarnings(suppressMessages(lme4::glmer(remission ~ x + visit + age + sex + baseline_bmi + baseline_hba1c + (1 | id),
                    data = d, family = stats::binomial))), error = function(e) NULL)
  if (is.null(fit)) return(c(beta = NA, se = NA, p = NA))
  co <- summary(fit)$coefficients["x", ]
  c(beta = co[["Estimate"]], se = co[["Std. Error"]], p = co[["Pr(>|z|)"]])
}
remission <- rbind_list(lapply(LAYERS, function(layer) {
  prep <- load_layer(layer); meta <- prep$metadata; ids <- time_responsive_ids(layer)
  base <- meta[meta$visit == "0", ]
  keep <- meta$visit %in% REMISSION_VISITS
  d <- data.frame(remission = ifelse(meta$visit == "0", 0, meta$remission), id = factor(meta$participant_id),
                  visit = factor(meta$visit), age = meta$age, sex = factor(meta$sex),
                  baseline_bmi = base$bmi[match(meta$participant_id, base$participant_id)],
                  baseline_hba1c = base$hba1c[match(meta$participant_id, base$participant_id)])[keep, ]
  d$visit <- droplevels(d$visit)
  fits <- par_map(ids, function(id, z, d) fit_remission_glmm(z[, id], d), z = prep$z[keep, , drop = FALSE], d = d, packages = "lme4")
  res <- data.frame(layer = layer, feature_id = ids, rbind_list(fits))
  res$fdr <- stats::p.adjust(res$p, "BH")   # within layer
  # Standardized beta: signed rank-based normal score within layer.
  res$std_beta <- sign(res$beta) * stats::qnorm(0.5 + 0.5 * (rank(abs(res$beta), na.last = "keep") - 0.5) / sum(is.finite(res$beta)))
  res
}))

# ---------------- Part 3: KEGG analysis of remission-associated features (Table S11.3) ----------------
# Background: all tested features of the layer with a KEGG mapping.
remission_enrichment <- rbind_list(lapply(c("pbmc_transcriptome", "serum_proteome", "serum_metabolome", "fecal_metabolome"), function(layer) {
  r <- remission[remission$layer == layer, ]
  prep <- load_layer(layer)
  if (layer == "fecal_metabolome" && is.null(prep$annotation$kegg_id)) prep$annotation$kegg_id <- match_kegg_by_name(prep$annotation$feature_name)
  kegg_ora_layer(layer, list(Remission = r$feature_id[is_sig(r$fdr)]), colnames(prep$z), prep$annotation)
}))

# ---------------- Part 4: agreement with the longitudinal HbA1c association (Figure 6D; Table S11.5) ----------------
# Concordant: higher odds of remission with a negative HbA1c association, or lower odds with a positive association.
model2 <- read_result("clinical_association", "model2.rds")$model2
hba1c <- model2[model2$trait == "hba1c", ]
concordance <- remission[is_sig(remission$fdr), ]
concordance$beta_hba1c <- hba1c$beta[match(paste(concordance$layer, concordance$feature_id), paste(hba1c$layer, hba1c$feature_id))]
concordance$concordant <- sign(concordance$beta) == -sign(concordance$beta_hba1c)
concordance_by_layer <- stats::aggregate(concordant ~ layer, concordance, function(x) 100 * mean(x))

# ---------------- Part 5: gut species by remission status (Figure 6B-C; Table S11.4) ----------------
meta <- load_layer("fecal_metagenome")$metadata
species <- read.csv(file.path(DATA_DIR, "microbiome", "species_relative_abundance.csv"), check.names = FALSE)
rel <- as.matrix(species[match(meta$sample_id, species$sample_id), -1])
rel <- rel / rowSums(rel)
status_w24 <- meta$remission[meta$visit == "24"][match(meta$participant_id, meta$participant_id[meta$visit == "24"])]

# Week 24: Wilcoxon rank-sum test for species with a mean relative abundance of at least 0.1%.
w24 <- meta$visit == "24"
tested <- colnames(rel)[colMeans(rel) >= MIN_MEAN_ABUNDANCE]
species_w24 <- rbind_list(lapply(tested, function(s) {
  a <- rel[w24 & status_w24 %in% 1, s]; b <- rel[w24 & status_w24 %in% 0, s]
  data.frame(species = s, mean_remission = mean(a), mean_no_remission = mean(b),
             log2_fold_change = log2((mean(a) + 1e-5) / (mean(b) + 1e-5)), p = suppressWarnings(stats::wilcox.test(a, b)$p.value))
}))
species_w24$fdr <- stats::p.adjust(species_w24$p, "BH")

# Mean composition of the most abundant species by week-24 remission status and visit.
top <- names(sort(colMeans(rel), decreasing = TRUE))[seq_len(TOP_SPECIES)]
abundance <- cbind(rel[, top], Others = 1 - rowSums(rel[, top]))
composition <- stats::aggregate(abundance, list(remission_w24 = status_w24, visit = meta$visit), mean)

saveRDS(list(remission = remission, enrichment = remission_enrichment, concordance = concordance,
             concordance_by_layer = concordance_by_layer, species_w24 = species_w24, composition = composition),
        result_path("remission", "remission.rds"))
