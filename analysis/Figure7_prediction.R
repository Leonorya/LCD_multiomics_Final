# Figure 7 and Tables S12: prediction of week-24 BMI, HbA1c and remission from baseline profiles.
# Base model: age + sex + corresponding baseline clinical measure (HbA1c for HbA1c and remission, BMI for BMI);
# linear or logistic regression.
# Omics model: base model variables (unpenalized) + all baseline features of one omics layer (LASSO).
# Repeated nested cross-validation; preprocessing and penalty tuning within outer training folds.
# Input : results/preprocessed/<layer>.rds
# Output: results/prediction/<layer>.rds, results/prediction/comparison.rds
source("analysis/utils.R")

# ---------------- Part 1: settings ----------------
N_REPEATS <- 100    # repeats of outer stratified 5-fold cross-validation
OUTER_K   <- 5
INNER_K   <- 4      # inner folds for LASSO penalty selection
N_BOOT    <- 2000   # participant bootstrap resamples
OUTCOMES <- list(bmi       = list(baseline = "bmi",   family = "gaussian"),
                 hba1c     = list(baseline = "hba1c", family = "gaussian"),
                 remission = list(baseline = "hba1c", family = "binomial"))

# ---------------- Part 2: cross-validation functions ----------------
# Folds balanced across outcome classes (binary) or across the median split (continuous).
make_folds <- function(y, k, seed, binary) {
  set.seed(seed)
  strata <- if (binary) y else as.integer(y > stats::median(y))
  fold <- integer(length(y)); offset <- 0
  for (s in sort(unique(strata))) {
    idx <- which(strata == s); idx <- idx[sample.int(length(idx))]
    fold[idx] <- ((seq_along(idx) - 1 + offset) %% k) + 1
    offset <- (offset + length(idx)) %% k
  }
  fold
}
# Median imputation, removal of constant predictors and standardization using the training rows only.
prepare_fold <- function(train, test) {
  med <- apply(train, 2, stats::median, na.rm = TRUE)
  for (j in seq_len(ncol(train))) { train[is.na(train[, j]), j] <- med[j]; test[is.na(test[, j]), j] <- med[j] }
  center <- colMeans(train); sds <- apply(train, 2, stats::sd); keep <- sds > 0
  list(train = scale(train[, keep, drop = FALSE], center[keep], sds[keep]),
       test  = scale(test[, keep, drop = FALSE], center[keep], sds[keep]))
}
fit_base <- function(train, y, test, family) {
  fit <- stats::glm(y ~ ., data = data.frame(y = y, train), family = family)
  as.numeric(stats::predict(fit, newdata = data.frame(test), type = "response"))
}
run_repeat <- function(r, y, predictors, family) {
  X <- predictors
  binary <- family == "binomial"
  outer_k <- if (binary) min(OUTER_K, min(table(y))) else OUTER_K
  folds <- make_folds(y, outer_k, SEED + r, binary)
  pred <- data.frame(repeat_id = r, row = seq_along(y), observed = y, Base = NA_real_, Omics = NA_real_)
  selected <- character()
  for (k in seq_len(outer_k)) {
    train <- folds != k; test <- folds == k
    b <- prepare_fold(X$base[train, , drop = FALSE], X$base[test, , drop = FALSE])
    pred$Base[test] <- fit_base(b$train, y[train], b$test, family)
    f <- prepare_fold(X$full[train, , drop = FALSE], X$full[test, , drop = FALSE])
    inner_k <- if (binary) min(INNER_K, min(table(y[train]))) else INNER_K
    clinical <- grepl("^clinical_", colnames(f$train))
    set.seed(SEED + 1000 * r + k)
    fit <- glmnet::cv.glmnet(f$train, y[train], family = family, alpha = 1, standardize = FALSE,
                             foldid = make_folds(y[train], inner_k, SEED + 1000 * r + k, binary),
                             type.measure = if (binary) "deviance" else "mse", penalty.factor = as.numeric(!clinical))
    pred$Omics[test] <- as.numeric(stats::predict(fit, newx = f$test, s = "lambda.min", type = "response"))
    co <- as.matrix(stats::coef(fit, s = "lambda.min"))
    selected <- c(selected, setdiff(rownames(co)[co[, 1] != 0], c("(Intercept)", colnames(f$train)[clinical])))
  }
  list(predictions = pred, selected = selected, n_fits = outer_k)
}
metric_values <- function(observed, predicted, binary) {
  if (binary) {
    cls <- predicted >= 0.5
    c(auc = as.numeric(pROC::auc(observed, predicted, levels = c(0, 1), direction = "<", quiet = TRUE)),
      sensitivity = mean(cls[observed == 1]), specificity = mean(!cls[observed == 0]))
  } else {
    c(pearson_r = stats::cor(observed, predicted), rmse = sqrt(mean((observed - predicted)^2)),
      mae = mean(abs(observed - predicted)), calibration_slope = unname(stats::coef(stats::lm(observed ~ predicted))[2]))
  }
}
repeat_metrics <- function(p, binary) sapply(c("Base", "Omics"), function(m) metric_values(p$observed, p[[m]], binary))

# Performance of participant-averaged out-of-fold predictions with 95% CIs (Table S12.1).
averaged_performance <- function(a, binary) {
  rbind_list(lapply(c("Base", "Omics"), function(m) {
    est <- metric_values(a$observed, a[[m]], binary)
    boot <- replicate(N_BOOT, {
      i <- if (binary) unlist(lapply(split(seq_along(a$observed), a$observed), function(g) g[sample.int(length(g), replace = TRUE)]))
           else sample.int(nrow(a), replace = TRUE)
      metric_values(a$observed[i], a[[m]][i], binary)
    })
    out <- data.frame(model = m, metric = names(est), estimate = est,
                      low95 = apply(boot, 1, stats::quantile, 0.025, na.rm = TRUE), high95 = apply(boot, 1, stats::quantile, 0.975, na.rm = TRUE))
    if (binary) out[out$metric == "auc", c("low95", "high95")] <-
      as.numeric(pROC::ci.auc(pROC::roc(a$observed, a[[m]], levels = c(0, 1), direction = "<", quiet = TRUE), method = "delong"))[c(1, 3)]
    out
  }))
}

# ---------------- Part 3: repeated nested cross-validation per layer and outcome (Figure 7B-E; Tables S12.1, S12.3-S12.5) ----------------
for (layer in LAYERS) {
  prep <- load_layer(layer)
  base <- prep$metadata$visit == "0"
  meta0 <- prep$metadata[base, ]
  omics0 <- prep$transformed[base, , drop = FALSE]
  colnames(omics0) <- paste0("omics_", colnames(omics0))
  w24 <- prep$metadata[prep$metadata$visit == "24", ]
  w24 <- w24[match(meta0$participant_id, w24$participant_id), ]

  results <- lapply(names(OUTCOMES), function(outcome) {
    o <- OUTCOMES[[outcome]]; binary <- o$family == "binomial"
    y <- w24[[outcome]]; ok <- !is.na(y)
    clinical <- cbind(clinical_age = meta0$age, clinical_sex = as.numeric(factor(meta0$sex)) - 1,
                      clinical_baseline = meta0[[o$baseline]])
    predictors <- list(base = clinical[ok, ], full = cbind(clinical, omics0)[ok, ])
    reps <- par_map(seq_len(N_REPEATS), function(r, y, predictors, family) run_repeat(r, y, predictors, family),
                    y = y[ok], predictors = predictors, family = o$family, packages = c("glmnet", "pROC"))
    predictions <- rbind_list(lapply(reps, `[[`, "predictions"))
    metrics <- lapply(split(predictions, predictions$repeat_id), repeat_metrics, binary = binary)
    metrics <- array(unlist(metrics), c(dim(metrics[[1]]), N_REPEATS), c(dimnames(metrics[[1]]), list(NULL)))

    # Out-of-fold predictions averaged over repeats for each participant.
    averaged <- stats::aggregate(cbind(observed, Base, Omics) ~ row, predictions, mean)
    averaged$participant_id <- meta0$participant_id[ok][averaged$row]
    set.seed(SEED)
    performance <- averaged_performance(averaged, binary)
    # Calibration of averaged remission probabilities: quartiles of predicted probability.
    calibration <- if (binary) rbind_list(lapply(c("Base", "Omics"), function(m) {
      bin <- ceiling(4 * rank(averaged[[m]], ties.method = "first") / nrow(averaged))
      data.frame(model = m, stats::aggregate(cbind(predicted = averaged[[m]], observed = averaged$observed), list(bin = bin), mean))
    })) else NULL
    n_fits <- sum(vapply(reps, `[[`, numeric(1), "n_fits"))
    selection <- sort(table(sub("^omics_", "", unlist(lapply(reps, `[[`, "selected")))) / n_fits, decreasing = TRUE)
    list(outcome = outcome, metrics_by_repeat = metrics, performance = performance, averaged = averaged,
         calibration = calibration, selection_frequency = selection)
  })
  saveRDS(stats::setNames(results, names(OUTCOMES)), result_path("prediction", paste0(layer, ".rds")))
}

# ---------------- Part 4: model comparisons on participant-averaged out-of-fold predictions ----------------
# Continuous outcomes: Williams' test for dependent overlapping correlations and bootstrap differences in R and RMSE (Table S12.2).
williams_test <- function(r12, r13, r23, n) {   # r12 = r(observed, model), r13 = r(observed, reference), r23 = r(model, reference)
  det_r <- 1 - r12^2 - r13^2 - r23^2 + 2 * r12 * r13 * r23; r_bar <- (r12 + r13) / 2
  t <- (r12 - r13) * sqrt((n - 1) * (1 + r23) / (2 * (n - 1) / (n - 3) * det_r + r_bar^2 * (1 - r23)^3))
  2 * stats::pt(-abs(t), n - 3)
}
rmse <- function(y, p) sqrt(mean((y - p)^2))
# Remission: paired DeLong test and DeLong 95% CI of the AUC (Table S12.3), bootstrap ROC bands (Figure 7D) and decision
# curves with non-remission as the event (Table S12.6).
net_benefit <- function(event, p, thresholds) vapply(thresholds, function(t)
  mean(p >= t & event == 1) - mean(p >= t & event == 0) * t / (1 - t), numeric(1))

continuous <- list(); remission <- list()
set.seed(SEED)
for (layer in LAYERS) {
  res <- read_result("prediction", paste0(layer, ".rds"))
  for (outcome in c("bmi", "hba1c")) {
    a <- res[[outcome]]$averaged; y <- a$observed
    boot <- replicate(N_BOOT, { i <- sample.int(length(y), replace = TRUE)
      c(stats::cor(y[i], a$Omics[i]) - stats::cor(y[i], a$Base[i]), rmse(y[i], a$Omics[i]) - rmse(y[i], a$Base[i])) })
    continuous[[paste(layer, outcome)]] <- data.frame(layer = layer, outcome = outcome, n = length(y),
      r_base = stats::cor(y, a$Base), r_omics = stats::cor(y, a$Omics),
      p_williams = williams_test(stats::cor(y, a$Omics), stats::cor(y, a$Base), stats::cor(a$Omics, a$Base), length(y)),
      delta_r = stats::cor(y, a$Omics) - stats::cor(y, a$Base),
      delta_r_low95 = stats::quantile(boot[1, ], 0.025), delta_r_high95 = stats::quantile(boot[1, ], 0.975),
      delta_rmse = rmse(y, a$Omics) - rmse(y, a$Base),
      delta_rmse_low95 = stats::quantile(boot[2, ], 0.025), delta_rmse_high95 = stats::quantile(boot[2, ], 0.975))
  }
  a <- res$remission$averaged
  roc <- lapply(c(Base = "Base", Omics = "Omics"), function(m) pROC::roc(a$observed, a[[m]], levels = c(0, 1), direction = "<", quiet = TRUE))
  auc_ci <- lapply(roc, function(r) as.numeric(pROC::ci.auc(r, method = "delong")))
  bands <- lapply(names(roc), function(m) {
    b <- pROC::ci.se(roc[[m]], specificities = seq(0, 1, 0.02), boot.n = N_BOOT, stratified = TRUE, progress = "none")
    data.frame(layer = layer, model = m, fpr = 1 - seq(0, 1, 0.02), se_low95 = b[, 1], se_high95 = b[, 3])
  })
  thresholds <- seq(0.01, 0.99, 0.01); event <- 1 - a$observed   # event = non-remission
  dca <- data.frame(layer = layer, threshold = thresholds,
                    Base = net_benefit(event, 1 - a$Base, thresholds), Omics = net_benefit(event, 1 - a$Omics, thresholds),
                    treat_all = mean(event) - (1 - mean(event)) * thresholds / (1 - thresholds), treat_none = 0)
  remission[[layer]] <- list(
    auc = data.frame(layer = layer, model = names(auc_ci), auc = sapply(auc_ci, `[`, 2), low95 = sapply(auc_ci, `[`, 1), high95 = sapply(auc_ci, `[`, 3)),
    p_delong = pROC::roc.test(roc$Base, roc$Omics, method = "delong", paired = TRUE)$p.value,
    roc_band = rbind_list(bands), dca = dca,
    useful_thresholds = sapply(c("Base", "Omics"), function(m) sum(dca[[m]] > pmax(dca$treat_all, dca$treat_none))))
}
saveRDS(list(continuous = rbind_list(continuous), remission = remission), result_path("prediction", "comparison.rds"))
