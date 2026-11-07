# Figure 5 and Tables S8-S10: longitudinal associations of time-responsive features with clinical traits, adiposity and
# glycemia feature sets (KEGG and correlation networks), adjustment for time-varying BMI, and mediation analyses.
# Input : results/preprocessed/<layer>.rds, results/time_model/<layer>.rds, results/navmix/<layer>.rds,
#         data/microbiome/species_kegg_pathways.csv
# Output: results/clinical_association/*.rds, results/mediation/*.rds
source("analysis/utils.R")

# ---------------- Part 1: settings ----------------
N_BOOT_MEDIATION <- 5000
N_BOOT_HIMA      <- 200
WINDOWS          <- c("12", "24")                  # change from baseline to week 12 and to week 24
EXPOSURES        <- c("delta_bmi", "delta_weight")

# ---------------- Part 2: longitudinal associations with clinical traits (Figure 5A-C; Tables S8.1-S8.2) ----------------
# clinical_z ~ omics_z + visit + age + sex + diet + (1 | participant); optionally + time-varying BMI (bmi_z).
fit_trait_lmm <- function(y, x, meta, bmi_z = NULL) {
  d <- model_frame(y, meta); d$x <- x
  if (!is.null(bmi_z)) d$bmi_z <- bmi_z
  d <- d[stats::complete.cases(d), ]; d$visit <- droplevels(d$visit)
  formula <- if (is.null(bmi_z)) y ~ x + visit + age + sex + diet + (1 | id) else y ~ x + bmi_z + visit + age + sex + diet + (1 | id)
  fit <- tryCatch(suppressMessages(lmerTest::lmer(formula, data = d, REML = FALSE, control = lme4::lmerControl(optimizer = "bobyqa"))),
                  error = function(e) NULL)
  if (is.null(fit)) return(c(beta = NA, se = NA, df = NA, p = NA))
  co <- summary(fit)$coefficients["x", ]
  c(beta = co[["Estimate"]], se = co[["Std. Error"]], df = co[["df"]], p = co[["Pr(>|t|)"]])
}
associate <- function(traits, adjust_bmi) rbind_list(lapply(LAYERS, function(layer) {
  prep <- load_layer(layer); ids <- time_responsive_ids(layer); meta <- prep$metadata
  bmi_z <- if (adjust_bmi) baseline_zscore_vector(meta$bmi, meta$visit) else NULL   # traits standardized within each layer's samples
  rbind_list(lapply(traits, function(trait) {
    y <- baseline_zscore_vector(meta[[trait]], meta$visit)
    fits <- par_map(ids, function(id, z, y, meta, bmi_z) fit_trait_lmm(y, z[, id], meta, bmi_z), z = prep$z, y = y, meta = meta, bmi_z = bmi_z)
    res <- data.frame(layer = layer, trait = trait, feature_id = ids, rbind_list(fits))
    res$fdr <- stats::p.adjust(res$p, "BH")   # within layer and trait
    res
  }))
}))
model2 <- associate(CLINICAL_TRAITS, adjust_bmi = FALSE)
significant <- model2[is_sig(model2$fdr), ]

# Adiposity- and glycemia-associated features (Table S9.1).
feature_sets <- unique(significant[, c("layer", "feature_id")])
feature_sets$node <- paste(feature_sets$layer, feature_sets$feature_id, sep = "::")
key <- paste(significant$layer, significant$feature_id, sep = "::")
feature_sets$adiposity <- feature_sets$node %in% key[significant$trait %in% ADIPOSITY_TRAITS]
feature_sets$glycemia  <- feature_sets$node %in% key[significant$trait %in% GLYCEMIC_TRAITS]
feature_sets$shared    <- feature_sets$adiposity & feature_sets$glycemia

# ---------------- Part 3: KEGG analysis of adiposity and glycemia feature sets (Figure 5D; Table S9.2) ----------------
# Background: all tested features of the layer with a KEGG mapping.
species_pathways <- read.csv(file.path(DATA_DIR, "microbiome", "species_kegg_pathways.csv"))
set_enrichment <- rbind_list(lapply(setdiff(LAYERS, "fecal_metatranscriptome"), function(layer) {
  f <- feature_sets[feature_sets$layer == layer, ]
  sets <- list(Adiposity = f$feature_id[f$adiposity], Glycemia = f$feature_id[f$glycemia])
  sets <- sets[lengths(sets) >= 2]
  if (!length(sets)) return(NULL)
  prep <- load_layer(layer)
  if (layer == "fecal_metabolome" && is.null(prep$annotation$kegg_id)) prep$annotation$kegg_id <- match_kegg_by_name(prep$annotation$feature_name)
  kegg_ora_layer(layer, sets, colnames(prep$z), prep$annotation, species_pathways)
}))

# ---------------- Part 4: correlation networks of the feature sets (Figure 5E-F; Table S9.3) ----------------
# One pooled BH adjustment over all pairs of associated features; each set uses the pairs among its own members.
tests <- correlate_pairs(wide_feature_table(feature_sets), feature_sets$node)
tests$cross_layer <- sub("::.*$", "", tests$node1) != sub("::.*$", "", tests$node2)
set_networks <- lapply(c(Adiposity = "adiposity", Glycemia = "glycemia"), function(flag) {
  nodes <- feature_sets$node[feature_sets[[flag]]]
  t <- tests[tests$node1 %in% nodes & tests$node2 %in% nodes, ]
  edges <- t[is_sig(t$fdr), ]
  layer_n <- table(sub("::.*$", "", nodes))
  cross_pairs <- (length(nodes)^2 - sum(layer_n^2)) / 2
  list(edges = edges, summary = data.frame(nodes = length(nodes), edges = nrow(edges), density = nrow(edges) / nrow(t),
       cross_edges = sum(edges$cross_layer), cross_density = sum(edges$cross_layer) / cross_pairs))
})
saveRDS(list(model2 = model2, feature_sets = feature_sets, set_enrichment = set_enrichment, set_networks = set_networks),
        result_path("clinical_association", "model2.rds"))

# ---------------- Part 5: adjustment for time-varying BMI (Tables S9.4-S9.5) ----------------
bmi_adjusted <- associate(GLYCEMIC_TRAITS, adjust_bmi = TRUE)
ref <- model2[match(paste(bmi_adjusted$layer, bmi_adjusted$trait, bmi_adjusted$feature_id),
                    paste(model2$layer, model2$trait, model2$feature_id)), ]
bmi_adjusted$beta_model2 <- ref$beta
bmi_adjusted$fdr_model2 <- ref$fdr
bmi_adjusted$change_pct <- 100 * (1 - bmi_adjusted$beta / bmi_adjusted$beta_model2)
bmi_adjusted$category <- ifelse(!is_sig(bmi_adjusted$fdr_model2), "Not significant in Model 2",
  ifelse(is_sig(bmi_adjusted$fdr) & sign(bmi_adjusted$beta) == sign(bmi_adjusted$beta_model2), "BMI-independent", "BMI-dependent"))
saveRDS(bmi_adjusted, result_path("clinical_association", "bmi_adjusted.rds"))

# ---------------- Part 6: mediation by trajectory-cluster scores (Tables S10.1-S10.2) ----------------
# Change in BMI or weight -> change in mediator -> change in HbA1c; linear models adjusted for baseline HbA1c, the baseline
# mediator, age, sex and diet group; indirect effect from participant bootstrap resamples.
window_data <- function(meta, mediator, end) {
  m <- data.frame(meta[, c("participant_id", "visit", "hba1c", "bmi", "weight", "age", "sex", "diet_group")], mediator = mediator)
  d <- merge(m[m$visit == "0", ], m[m$visit == end, c("participant_id", "hba1c", "bmi", "weight", "mediator")],
             by = "participant_id", suffixes = c("", "_end"))
  data.frame(participant_id = d$participant_id, delta_hba1c = d$hba1c_end - d$hba1c, delta_bmi = d$bmi_end - d$bmi,
             delta_weight = d$weight_end - d$weight, delta_m = d$mediator_end - d$mediator,
             baseline_hba1c = d$hba1c, baseline_m = d$mediator, age = d$age, sex = factor(d$sex), diet = factor(d$diet_group))
}
covariates <- function(serial = FALSE)
  c("baseline_hba1c", if (serial) c("baseline_m1", "baseline_m2") else "baseline_m", "age", "sex", "diet")
bootstrap_p <- function(x) { x <- x[is.finite(x)]; min(1, 2 * min((sum(x <= 0) + 1) / (length(x) + 1), (sum(x >= 0) + 1) / (length(x) + 1))) }

# paths: list of (response, predictors); effect(coefs) returns the indirect effect from the fitted path coefficients.
bootstrap_indirect <- function(d, paths, effect, n_boot, seed) {
  d <- d[stats::complete.cases(d[, unique(unlist(paths))]), ]
  designs <- lapply(paths, function(p) list(x = stats::model.matrix(stats::reformulate(p[-1]), d), y = d[[p[1]]]))
  estimate <- function(i) effect(lapply(designs, function(m) stats::lm.fit(m$x[i, , drop = FALSE], m$y[i])$coefficients))
  set.seed(seed)
  draws <- replicate(n_boot, estimate(sample.int(nrow(d), replace = TRUE)))
  ci <- stats::quantile(draws, c(0.025, 0.975), na.rm = TRUE)
  c(n = nrow(d), indirect = estimate(seq_len(nrow(d))), low95 = ci[[1]], high95 = ci[[2]], p = bootstrap_p(draws))
}
simple_mediation <- function(d, exposure, n_boot, seed) {
  cov <- covariates()
  bootstrap_indirect(d, list(c("delta_m", exposure, cov), c("delta_hba1c", exposure, "delta_m", cov)),
                     function(b) b[[1]][[exposure]] * b[[2]][["delta_m"]], n_boot, seed)
}

# Mediators: all trajectory-cluster scores (noise group excluded).
candidates <- list()
for (layer in LAYERS) {
  meta <- load_layer(layer)$metadata
  scores <- read_result("navmix", paste0(layer, ".rds"))$cluster_scores
  scores <- scores[match(meta$sample_id, scores$sample_id), ]
  for (end in WINDOWS) for (cl in setdiff(names(scores), c("sample_id", "participant_id", "visit")))
    candidates[[paste(layer, end, cl)]] <- list(layer = layer, window = end, cluster = cl, data = window_data(meta, scores[[cl]], end))
}
jobs <- expand.grid(candidate = names(candidates), exposure = EXPOSURES, stringsAsFactors = FALSE)
fits <- par_map(seq_len(nrow(jobs)), function(i, jobs, candidates, n_boot) {
  simple_mediation(candidates[[jobs$candidate[i]]]$data, jobs$exposure[i], n_boot, SEED + i)
}, jobs = jobs, candidates = candidates, n_boot = N_BOOT_MEDIATION, packages = character())
info <- rbind_list(lapply(candidates[jobs$candidate], function(c) data.frame(layer = c$layer, window = c$window, mediator = c$cluster)))
mediation <- data.frame(info, exposure = jobs$exposure, rbind_list(fits))
mediation$fdr <- stats::ave(mediation$p, mediation$layer, mediation$window, mediation$exposure, FUN = function(p) stats::p.adjust(p, "BH"))

# ---------------- Part 7: serial mediation, metagenomic cluster -> fecal metabolite cluster (Table S10.2) ----------------
serial_mediation <- function(d, exposure, n_boot, seed) {
  cov <- covariates(serial = TRUE)
  bootstrap_indirect(d, list(c("delta_m1", exposure, cov), c("delta_m2", exposure, "delta_m1", cov),
                             c("delta_hba1c", exposure, "delta_m1", "delta_m2", cov)),
                     function(b) b[[1]][[exposure]] * b[[2]][["delta_m1"]] * b[[3]][["delta_m2"]], n_boot, seed)
}
pick <- function(layer, end) Filter(function(c) c$layer == layer && c$window == end, candidates)
serial_jobs <- rbind_list(lapply(WINDOWS, function(end) expand.grid(window = end,
  m1 = names(pick("fecal_metagenome", end)), m2 = names(pick("fecal_metabolome", end)), exposure = EXPOSURES, stringsAsFactors = FALSE)))
serial_fits <- par_map(seq_len(nrow(serial_jobs)), function(i, jobs, candidates, n_boot) {
  a <- candidates[[jobs$m1[i]]]$data; b <- candidates[[jobs$m2[i]]]$data
  d <- merge(a, b[, c("participant_id", "delta_m", "baseline_m")], by = "participant_id", suffixes = c("1", "2"))
  serial_mediation(d, jobs$exposure[i], n_boot, SEED + 100000 + i)
}, jobs = serial_jobs, candidates = candidates, n_boot = N_BOOT_MEDIATION, packages = character())
serial <- data.frame(serial_jobs, rbind_list(serial_fits))
serial[c("m1", "m2")] <- lapply(serial[c("m1", "m2")], function(k) sub("^.* ", "", k))   # cluster names
serial$fdr <- stats::ave(serial$p, serial$window, serial$exposure, FUN = function(p) stats::p.adjust(p, "BH"))
saveRDS(list(mediation = mediation, serial = serial), result_path("mediation", "cluster_mediation.rds"))

# ---------------- Part 8: high-dimensional mediation with HIMA (Tables S10.3-S10.4) ----------------
# De-biased lasso HIMA: sure independence screening of the top ceiling(2n / log n) mediators, de-biased lasso for the
# outcome model, least squares for the mediator models and joint significance by HDMT FDR. Mediators: change in every
# time-responsive feature of a layer, residualized on its baseline level. Covariates: baseline HbA1c, age, sex, diet group.
hima_dblasso <- function(X, M, Y, COV) {
  n <- nrow(M); p <- ncol(M)
  X <- HIMA:::process_var(X, TRUE); M <- HIMA:::process_var(M, TRUE); COV <- HIMA:::process_var(COV, TRUE)
  MZX <- cbind(M, COV, X); XZ <- cbind(X, COV); q <- ncol(COV)
  d0 <- min(p, ceiling(2 * n / log(n)))
  beta_sis  <- vapply(seq_len(p), function(i) stats::lsfit(MZX[, c(i, (p + 1):(p + q + 1))], Y)$coefficients[2], numeric(1))
  alpha_sis <- vapply(seq_len(p), function(i) stats::lsfit(XZ, M[, i])$coefficients[2], numeric(1))
  id <- order(-abs(alpha_sis * beta_sis))[seq_len(d0)]
  dl <- suppressMessages(hdi::lasso.proj(x = MZX[, c(id, (p + 1):(p + q + 1))], y = Y, family = "gaussian"))
  a <- t(vapply(id, function(i) { f <- stats::lsfit(XZ, M[, i]); e <- stats::coef(f)[2]; s <- stats::ls.diag(f)$std.err[2]
    c(e, 2 * (1 - stats::pnorm(abs(e) / s))) }, numeric(2)))
  PA <- cbind(a[, 2], dl$pval[seq_along(id)]) + matrix(stats::runif(2 * length(id), 0, 1e-10), length(id), 2)
  np <- HIMA:::null_estimation(PA)
  fdr <- HDMT::fdr_est(np$alpha00, np$alpha01, np$alpha10, np$alpha1, np$alpha2, PA, exact = 0)
  data.frame(feature_id = colnames(M)[id], alpha = a[, 1], beta = dl$bhat[seq_along(id)], hdmt_fdr = as.numeric(fdr),
             significant = as.numeric(fdr) <= 0.05)
}
hima_data <- function(prep, ids, end) {
  m <- prep$metadata; b <- which(m$visit == "0"); e <- which(m$visit == end)
  pid <- intersect(m$participant_id[b], m$participant_id[e])
  b <- b[match(pid, m$participant_id[b])]; e <- e[match(pid, m$participant_id[e])]
  d <- data.frame(delta_hba1c = m$hba1c[e] - m$hba1c[b], delta_bmi = m$bmi[e] - m$bmi[b], delta_weight = m$weight[e] - m$weight[b],
                  baseline_hba1c = m$hba1c[b], age = m$age[b], sex = as.numeric(factor(m$sex[b])) - 1, diet = as.numeric(factor(m$diet_group[b])) - 1)
  dm <- prep$z[e, ids, drop = FALSE] - prep$z[b, ids, drop = FALSE]; bm <- prep$z[b, ids, drop = FALSE]
  keep <- stats::complete.cases(d) & stats::complete.cases(dm)
  res <- vapply(seq_along(ids), function(j) stats::residuals(stats::lm(dm[keep, j] ~ bm[keep, j])), numeric(sum(keep)))
  colnames(res) <- ids
  list(d = d[keep, ], M = res[, apply(res, 2, stats::sd) > 0, drop = FALSE])
}
hima <- list()
for (layer in LAYERS) for (end in WINDOWS) for (exposure in EXPOSURES) {
  w <- hima_data(load_layer(layer), time_responsive_ids(layer), end)
  x_values <- w$d[[exposure]]; y_values <- w$d$delta_hba1c
  cov_values <- as.matrix(w$d[, c("baseline_hba1c", "age", "sex", "diet")])
  set.seed(SEED)
  main <- hima_dblasso(x_values, w$M, y_values, cov_values)
  # Stability across participant bootstrap resamples.
  boots <- par_map(seq_len(N_BOOT_HIMA), function(r, exposure, mediators, outcome, cov) {
    set.seed(SEED + r); i <- sample.int(length(outcome), replace = TRUE)
    f <- tryCatch(hima_dblasso(exposure[i], mediators[i, , drop = FALSE], outcome[i], cov[i, , drop = FALSE]), error = function(e) NULL)
    if (is.null(f)) NULL else f$feature_id[f$significant]
  }, exposure = x_values, mediators = w$M, outcome = y_values, cov = cov_values, packages = "HIMA")
  valid <- !vapply(boots, is.null, logical(1))
  freq <- table(unlist(boots[valid])) / sum(valid)
  main$bootstrap_significant_freq <- ifelse(is.na(freq[main$feature_id]), 0, freq[main$feature_id])
  hima[[paste(layer, end, exposure)]] <- data.frame(layer = layer, window = end, exposure = exposure, n = length(y_values),
                                                    n_boot_valid = sum(valid), main)
}
saveRDS(rbind_list(hima), result_path("mediation", "hima.rds"))
