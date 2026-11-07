# Figure 1 and Tables S1, S3 and S4: clinical outcomes, time effects of omics features, visit localization and
# effect modification by diet group, sex and pre-intervention metformin use.
# Input : data/clinical.csv, data/medication.csv, results/preprocessed/<layer>.rds
# Output: results/clinical/*.rds, results/time_model/<layer>.rds, results/effect_modification/features.rds
source("analysis/utils.R")

# ---------------- Part 1: settings ----------------
BASELINE_CONTINUOUS <- c("age", "diabetes_duration", "weight", "bmi", "waist", "body_fat_mass", "body_fat_pct", "muscle_mass",
                         "muscle_pct", "liver_fat", "pancreatic_fat", "hba1c", "fpg", "homa2_ir", "homa2_b", "homa2_s")
BASELINE_CATEGORICAL <- c("sex", "education", "n_glucose_lowering_drugs", "metformin", "family_history_diabetes",
                          "hypertension", "hyperlipidemia", "diet_group")
CATEGORICAL_S1.2 <- c("sex", "education", "n_glucose_lowering_drugs", "family_history_diabetes", "hypertension", "hyperlipidemia")
CATEGORICAL_S1.3 <- c("sex", "diet_group")
CHANGE_VARS <- c("weight", "bmi", "waist", "body_fat_pct", "muscle_mass", "liver_fat", "pancreatic_fat", "hba1c", "fpg", "homa2_ir")
MODIFIERS <- c("diet", "sex", "metformin")

clinical <- read.csv(file.path(DATA_DIR, "clinical.csv"))
clinical$visit <- as.character(clinical$visit)
medication <- read.csv(file.path(DATA_DIR, "medication.csv"))   # participant_id, metformin (0/1, before the intervention)

# ---------------- Part 2: clinical outcomes (Figure 1B; Tables S1.1-S1.3) ----------------
blood_ids <- unique(load_layer("serum_proteome")$metadata$participant_id)
stool_ids <- unique(load_layer("fecal_metagenome")$metadata$participant_id)
cohort <- clinical[clinical$participant_id %in% blood_ids, ]
baseline <- cohort[cohort$visit == "0", ]
baseline$metformin <- medication$metformin[match(baseline$participant_id, medication$participant_id)]
baseline$remission_w24 <- cohort$remission[cohort$visit == "24"][match(baseline$participant_id, cohort$participant_id[cohort$visit == "24"])]
baseline$fecal_included <- baseline$participant_id %in% stool_ids

# Week 24 versus baseline: paired t tests in all participants.
clinical_change <- rbind_list(lapply(CHANGE_VARS, function(v) {
  w0 <- cohort[cohort$visit == "0", c("participant_id", v)]; w24 <- cohort[cohort$visit == "24", c("participant_id", v)]
  m <- stats::na.omit(merge(w0, w24, by = "participant_id", suffixes = c("_w0", "_w24")))
  data.frame(variable = v, n = nrow(m), mean_w0 = mean(m[[2]]), mean_w24 = mean(m[[3]]),
             p = stats::t.test(m[[3]], m[[2]], paired = TRUE)$p.value)
}))

# Baseline characteristics (Tables S1.1-S1.3).
describe <- function(d) rbind(
  rbind_list(lapply(BASELINE_CONTINUOUS, function(v)
    data.frame(variable = v, level = "", value = sprintf("%.1f (%.1f)", mean(d[[v]], na.rm = TRUE), stats::sd(d[[v]], na.rm = TRUE))))),
  rbind_list(lapply(BASELINE_CATEGORICAL, function(v) { n <- table(d[[v]])
    data.frame(variable = v, level = names(n), value = sprintf("%d (%.1f)", as.integer(n), 100 * as.numeric(n) / sum(n))) })))
compare_baseline <- function(group, welch, categorical_vars, categorical_test) {
  continuous <- rbind_list(lapply(BASELINE_CONTINUOUS, function(v)
    data.frame(variable = v, p = stats::t.test(baseline[[v]] ~ group, var.equal = !welch)$p.value)))
  categorical <- rbind_list(lapply(categorical_vars, function(v)
    data.frame(variable = v, p = categorical_test(table(baseline[[v]], group))$p.value)))
  rbind(continuous, categorical)
}
baseline_all       <- describe(baseline)
baseline_by_fecal  <- lapply(split(baseline, baseline$fecal_included), describe)
baseline_by_remission <- lapply(split(baseline, baseline$remission_w24), describe)
baseline_fecal     <- compare_baseline(baseline$fecal_included, welch = FALSE, CATEGORICAL_S1.2, stats::chisq.test)
baseline_remission <- compare_baseline(baseline$remission_w24, welch = TRUE, CATEGORICAL_S1.3, stats::fisher.test)

# Remission and weight loss of at least 10% at weeks 12 and 24.
base_weight <- stats::setNames(baseline$weight, baseline$participant_id)
outcomes <- rbind_list(lapply(c("12", "24"), function(w) {
  d <- cohort[cohort$visit == w, ]
  data.frame(visit = w, remission_pct = 100 * mean(d$remission == 1, na.rm = TRUE),
             weight_loss10_pct = 100 * mean((base_weight[d$participant_id] - d$weight) / base_weight[d$participant_id] >= 0.10, na.rm = TRUE))
}))
saveRDS(list(change = clinical_change, baseline_all = baseline_all, baseline_by_fecal = baseline_by_fecal,
             baseline_by_remission = baseline_by_remission, baseline_fecal = baseline_fecal, baseline_remission = baseline_remission,
             outcomes = outcomes, cohort = cohort), result_path("clinical", "clinical_outcomes.rds"))

# ---------------- Part 3: time model and visit localization (Figure 1C-E; Tables S3) ----------------
for (layer in LAYERS) {
  prep <- load_layer(layer)
  fits <- par_map(colnames(prep$z), function(id, z, meta) fit_time_lmm(z[, id], meta), z = prep$z, meta = prep$metadata)
  res <- data.frame(feature_id = colnames(prep$z), layer = layer, rbind_list(fits))
  res$feature_name <- prep$annotation$feature_name[match(res$feature_id, prep$annotation$feature_id)]
  res$fdr_global <- stats::p.adjust(res$p_global, "BH")
  res$time_responsive <- is_sig(res$fdr_global)

  # Visits at which time-responsive features changed: contrast with baseline at P < 0.05.
  for (w in c(1, 12, 24)) res[[paste0("selected_w", w)]] <- res$time_responsive & is_sig(res[[paste0("p_w", w)]])
  selected <- as.matrix(res[, c("selected_w1", "selected_w12", "selected_w24")])
  res$visit_combination <- apply(selected, 1, function(s) if (any(s)) paste(c("W1", "W12", "W24")[s], collapse = "+") else "None")
  saveRDS(res, result_path("time_model", paste0(layer, ".rds")))
}

# ---------------- Part 4: effect modification of time-responsive features (Tables S4.1, S4.3, S4.4) ----------------
# Trajectory-cluster scores are tested in Figure3_fecal_diversity_trajectories.R after all layers are clustered.
modification <- rbind_list(lapply(LAYERS, function(layer) {
  prep <- load_layer(layer)
  meta <- prep$metadata
  meta$metformin <- medication$metformin[match(meta$participant_id, medication$participant_id)]
  ids <- time_responsive_ids(layer)
  rbind_list(lapply(MODIFIERS, function(modifier) {
    fits <- par_map(ids, function(id, z, meta, modifier) fit_time_modifier(z[, id], meta, modifier),
                    z = prep$z, meta = meta, modifier = modifier)
    res <- data.frame(layer = layer, modifier = modifier, feature_id = ids, rbind_list(fits))
    res$fdr_interaction <- stats::p.adjust(res$p_interaction, "BH")   # within layer and modifier
    res
  }))
}))
modification_summary <- rbind_list(lapply(split(modification, list(modification$layer, modification$modifier), drop = TRUE),
  function(d) data.frame(layer = d$layer[1], modifier = d$modifier[1], summarize_modification(d))))

# Pattern behind time x diet interactions: mean sign-aligned level by diet group and visit (Table S4.4).
interaction_pattern <- rbind_list(lapply(LAYERS, function(layer) {
  sig <- modification$feature_id[modification$layer == layer & modification$modifier == "diet" & is_sig(modification$fdr_interaction)]
  if (!length(sig)) return(NULL)
  prep <- load_layer(layer); tm <- load_time(layer)
  z <- sweep(prep$z[, sig, drop = FALSE], 2, sign(tm$beta_w24[match(sig, tm$feature_id)]), "*")
  m <- stats::aggregate(list(mean_z = rowMeans(z, na.rm = TRUE)),
                        list(diet_group = prep$metadata$diet_group, visit = prep$metadata$visit), mean, na.rm = TRUE)
  data.frame(layer = layer, n_features = length(sig), m)
}))
saveRDS(list(features = modification, summary = modification_summary, diet_pattern = interaction_pattern),
        result_path("effect_modification", "features.rds"))
