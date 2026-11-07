# Shared settings and helper functions sourced by every analysis script.
# Run all scripts from the repository root, e.g. `Rscript analysis/Figure1_time_effects.R`.

# ---------------- Part 1: settings ----------------
DATA_DIR    <- "data"
RESULTS_DIR <- "results"
N_CORES     <- 10                       # maximum number of parallel workers
SEED        <- 2026
VISITS      <- c("0", "1", "12", "24")  # study weeks

LAYERS <- c("pbmc_transcriptome", "serum_proteome", "serum_metabolome",
            "fecal_metagenome", "fecal_metatranscriptome", "fecal_metabolome")
BLOOD_LAYERS <- LAYERS[1:3]
STOOL_LAYERS <- LAYERS[4:6]
LAYER_TRANSFORM <- c(pbmc_transcriptome = "log2p1", serum_proteome = "log2", serum_metabolome = "log2",
                     fecal_metagenome = "clr", fecal_metatranscriptome = "clr", fecal_metabolome = "none")
NAVMIX_K <- c(pbmc_transcriptome = 5, serum_proteome = 3, serum_metabolome = 5,
              fecal_metagenome = 7, fecal_metatranscriptome = 3, fecal_metabolome = 5)

CLINICAL_TRAITS  <- c("bmi", "liver_fat", "pancreatic_fat", "hba1c", "fpg", "homa2_ir")
ADIPOSITY_TRAITS <- CLINICAL_TRAITS[1:3]
GLYCEMIC_TRAITS  <- CLINICAL_TRAITS[4:6]

# ---------------- Part 2: input/output ----------------
result_path <- function(...) {
  path <- file.path(RESULTS_DIR, ...)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  path
}
read_result <- function(...) readRDS(file.path(RESULTS_DIR, ...))
load_layer  <- function(layer) read_result("preprocessed", paste0(layer, ".rds"))
load_time   <- function(layer) read_result("time_model", paste0(layer, ".rds"))
time_responsive_ids <- function(layer) with(load_time(layer), feature_id[time_responsive])

# Run fun over x on up to N_CORES workers. All functions and settings defined in the global environment are copied to
# the workers; data are passed through `...`.
par_map <- function(x, fun, ..., packages = c("lme4", "lmerTest")) {
  cl <- parallel::makeCluster(min(N_CORES, parallel::detectCores()))
  on.exit(parallel::stopCluster(cl))
  parallel::clusterCall(cl, function(p) for (pkg in p) suppressMessages(library(pkg, character.only = TRUE)), packages)
  shared <- Filter(function(n) is.function(get(n, globalenv())) || grepl("^[A-Z_0-9]+$", n), ls(globalenv()))
  parallel::clusterExport(cl, shared, envir = globalenv())
  parallel::parLapply(cl, x, fun, ...)
}
rbind_list <- function(x) { out <- do.call(rbind, x); rownames(out) <- NULL; out }
is_sig <- function(p, alpha = 0.05) !is.na(p) & p < alpha

# ---------------- Part 3: preprocessing ----------------
transform_layer <- function(x, method) {
  switch(method,
    log2p1 = log2(x + 1),
    log2   = log2(x),
    clr    = { lx <- log(x); lx - rowMeans(lx, na.rm = TRUE) },
    none   = x)
}

# Standardize each feature by the mean and SD of the baseline samples.
baseline_zscore <- function(x, visit) {
  base <- as.character(visit) == "0"
  center <- colMeans(x[base, , drop = FALSE], na.rm = TRUE)
  scale <- apply(x[base, , drop = FALSE], 2, stats::sd, na.rm = TRUE)
  zero_variance <- !is.finite(scale) | scale <= 0
  z <- sweep(sweep(x, 2, center, "-"), 2, scale, "/")
  list(z = z[, !zero_variance, drop = FALSE], zero_variance = zero_variance)
}
baseline_zscore_vector <- function(x, visit) {
  base <- as.character(visit) == "0" & is.finite(x)
  (x - mean(x[base])) / stats::sd(x[base])
}

# ---------------- Part 4: longitudinal mixed models ----------------
model_frame <- function(y, meta) {
  data.frame(y = as.numeric(y), id = factor(meta$participant_id),
             visit = factor(as.character(meta$visit), levels = VISITS),
             age = meta$age, sex = factor(meta$sex), diet = factor(meta$diet_group))
}

# Time model: omics_z ~ visit + age + sex + diet + (1 | participant); type III F test of visit (Satterthwaite).
fit_time_lmm <- function(y, meta) {
  d <- model_frame(y, meta)
  d <- d[stats::complete.cases(d), ]
  d$visit <- droplevels(d$visit)
  out <- stats::setNames(rep(NA_real_, 10), c(outer(c("beta_w", "se_w", "p_w"), c(1, 12, 24), paste0), "p_global"))
  fit <- tryCatch(suppressMessages(lmerTest::lmer(y ~ visit + age + sex + diet + (1 | id), data = d, REML = FALSE,
                    control = lme4::lmerControl(optimizer = "bobyqa"))), error = function(e) NULL)
  if (is.null(fit)) return(out)
  coefs <- summary(fit)$coefficients
  for (w in c(1, 12, 24)) {
    term <- paste0("visit", w)
    if (term %in% rownames(coefs)) out[paste0(c("beta_w", "se_w", "p_w"), w)] <- coefs[term, c("Estimate", "Std. Error", "Pr(>|t|)")]
  }
  out["p_global"] <- stats::anova(fit, type = 3, ddf = "Satterthwaite")["visit", "Pr(>F)"]
  out
}

# Effect modification: omics_z ~ visit * modifier + age + other covariates + (1 | participant).
# modifier is "diet", "sex" or "metformin"; returns the interaction P value and the week-24 change within each level.
fit_time_modifier <- function(y, meta, modifier) {
  d <- model_frame(y, meta)
  d$metformin <- factor(if (is.null(meta$metformin)) NA else meta$metformin, levels = c(0, 1))
  d <- d[stats::complete.cases(d[, c("y", "id", "visit", "age", "sex", "diet", modifier)]), ]
  d$visit <- droplevels(d$visit)
  covariates <- setdiff(c("sex", "diet"), modifier)
  formula <- stats::reformulate(c(paste("visit *", modifier), "age", covariates, "(1 | id)"), response = "y")
  out <- c(p_interaction = NA, beta_w24_level1 = NA, se_w24_level1 = NA, beta_w24_level2 = NA, se_w24_level2 = NA,
           beta_w24_difference = NA, se_w24_difference = NA)
  fit <- tryCatch(suppressWarnings(suppressMessages(lmerTest::lmer(formula, data = d, REML = FALSE,
                    control = lme4::lmerControl(optimizer = "bobyqa")))), error = function(e) NULL)
  if (is.null(fit)) return(out)
  out["p_interaction"] <- stats::anova(fit, type = 3)[paste0("visit:", modifier), "Pr(>F)"]
  b <- lme4::fixef(fit); V <- as.matrix(stats::vcov(fit))
  interaction <- paste0("visit24:", modifier, levels(d[[modifier]])[2])
  contrast <- function(w) { v <- stats::setNames(rep(0, length(b)), names(b)); v[names(w)] <- w
    c(sum(v * b), sqrt(drop(t(v) %*% V %*% v))) }
  out[c("beta_w24_level1", "se_w24_level1")] <- contrast(c(visit24 = 1))
  out[c("beta_w24_level2", "se_w24_level2")] <- contrast(stats::setNames(c(1, 1), c("visit24", interaction)))
  out[c("beta_w24_difference", "se_w24_difference")] <- contrast(stats::setNames(1, interaction))
  out
}

# Summary of effect modification across features: tested, significant interactions and week-24 direction agreement.
summarize_modification <- function(d) {
  data.frame(tested = nrow(d), interaction_fdr05 = sum(d$fdr_interaction < 0.05, na.rm = TRUE),
             correlation_w24 = stats::cor(d$beta_w24_level1, d$beta_w24_level2, use = "complete.obs"),
             same_direction_pct = 100 * mean(sign(d$beta_w24_level1) == sign(d$beta_w24_level2), na.rm = TRUE))
}

# ---------------- Part 5: trajectory clustering (NAvMix) ----------------
# Clusters the mean standardized trajectories of time-responsive features; label K + 1 is the noise group.
run_navmix <- function(prep, feature_ids, k, seed) {
  trajectories <- sapply(VISITS, function(v) colMeans(prep$z[prep$metadata$visit == v, feature_ids, drop = FALSE], na.rm = TRUE))
  set.seed(seed)
  fixed <- navmix::navmix(trajectories, K = k, select_K = FALSE, pj_ini = 0.05, plot = FALSE, plot_radial = FALSE, reorder_traits = FALSE)
  set.seed(seed + 100)
  profile <- navmix::navmix(trajectories, K = min(20, nrow(trajectories) - 1), select_K = TRUE, pj_ini = 0.05,
                            plot = FALSE, plot_radial = FALSE, reorder_traits = FALSE)
  raw <- fixed$fit$z
  assignments <- data.frame(feature_id = feature_ids, cluster = ifelse(raw > k, "Noise", paste0("C", raw)))
  clusters <- sort(setdiff(unique(assignments$cluster), "Noise"))
  scores <- prep$metadata[, c("sample_id", "participant_id", "visit")]
  for (cl in clusters) scores[[cl]] <- rowMeans(prep$z[, assignments$feature_id[assignments$cluster == cl], drop = FALSE], na.rm = TRUE)
  list(assignments = assignments, trajectories = trajectories, cluster_scores = scores,
       bic_profile = data.frame(k = seq_along(profile$BIC), bic = profile$BIC))
}

# Features changed at all three post-baseline visits: direction agreement and median absolute beta per visit.
all_visit_summary <- function(layer) {
  d <- load_time(layer)
  d <- d[d$selected_w1 & d$selected_w12 & d$selected_w24, ]
  b <- as.matrix(d[, c("beta_w1", "beta_w12", "beta_w24")])
  data.frame(layer = layer, n_features = nrow(d), same_direction = sum(abs(rowSums(sign(b))) == 3),
             median_abs_beta_w1 = stats::median(abs(b[, 1])), median_abs_beta_w12 = stats::median(abs(b[, 2])),
             median_abs_beta_w24 = stats::median(abs(b[, 3])))
}

# ---------------- Part 6: KEGG over-representation analysis ----------------
# KEGG pathway membership for human genes (Entrez IDs) or compounds (KEGG compound IDs).
kegg_terms <- function(type) {
  links <- if (type == "gene") KEGGREST::keggLink("pathway", "hsa") else KEGGREST::keggLink("pathway", "compound")
  names_vec <- if (type == "gene") KEGGREST::keggList("pathway", "hsa") else KEGGREST::keggList("pathway")
  list(term2gene = unique(data.frame(term = sub("^path:", "", unname(links)), gene = sub("^(hsa|cpd):", "", names(links)))),
       term2name = data.frame(term = sub("^path:", "", names(names_vec)), name = sub(" - Homo sapiens \\(human\\)$", "", unname(names_vec))))
}

# Map features to KEGG identifiers: transcripts (Ensembl) and proteins (gene symbols) to Entrez IDs;
# metabolites use the KEGG compound IDs in the annotation table.
map_to_kegg <- function(layer, annotation) {
  if (layer == "pbmc_transcriptome" || layer == "serum_proteome") {
    keytype <- if (layer == "pbmc_transcriptome") "ENSEMBL" else "SYMBOL"
    key <- if (layer == "pbmc_transcriptome") sub("\\..*$", "", annotation$feature_id) else annotation$gene_symbol
    m <- AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, keys = unique(stats::na.omit(key)), keytype = keytype, columns = "ENTREZID")
    out <- data.frame(feature_id = annotation$feature_id, kegg_id = m$ENTREZID[match(key, m[[keytype]])])
  } else {
    out <- data.frame(feature_id = annotation$feature_id, kegg_id = annotation$kegg_id)
  }
  out[!is.na(out$kegg_id) & out$kegg_id != "", ]
}

# Fecal metabolites: KEGG compound IDs by exact matching of normalized names to KEGG compound names and synonyms.
match_kegg_by_name <- function(feature_names) {
  compounds <- KEGGREST::keggList("compound")
  normalize <- function(x) gsub("[^a-z0-9]", "", tolower(trimws(x)))
  synonyms <- strsplit(unname(compounds), ";")
  catalogue <- unique(data.frame(name = normalize(unlist(synonyms)),
                                 kegg_id = rep(sub("^cpd:", "", names(compounds)), lengths(synonyms))))
  ambiguous <- unique(catalogue$name[duplicated(catalogue$name)])   # names shared by several compounds are not used
  catalogue <- catalogue[!catalogue$name %in% ambiguous, ]
  catalogue$kegg_id[match(normalize(feature_names), catalogue$name)]
}

# Hypergeometric test (clusterProfiler) for each feature set against all tested features of the layer with a KEGG mapping.
run_kegg_ora <- function(feature_sets, mapping, universe, terms) {
  universe_ids <- unique(mapping$kegg_id[mapping$feature_id %in% universe])
  rbind_list(lapply(names(feature_sets), function(set) {
    ids <- unique(mapping$kegg_id[mapping$feature_id %in% feature_sets[[set]]])
    if (length(ids) < 2) return(NULL)
    e <- as.data.frame(clusterProfiler::enricher(ids, universe = universe_ids, TERM2GENE = terms$term2gene, TERM2NAME = terms$term2name,
                                                  pvalueCutoff = 1, qvalueCutoff = 1, pAdjustMethod = "BH", minGSSize = 2, maxGSSize = 1000))
    if (!nrow(e)) return(NULL)
    data.frame(set = set, pathway_id = e$ID, pathway = e$Description, gene_ratio = e$GeneRatio,
               count = e$Count, p = e$pvalue, fdr = e$p.adjust, kegg_ids = e$geneID)
  }))
}

# Metagenomic species: hypergeometric test against species-level KEGG pathway annotations of the gene catalogue.
run_species_ora <- function(feature_sets, universe, species_pathways) {
  sp <- unique(species_pathways[species_pathways$feature_id %in% universe, c("feature_id", "pathway_id", "pathway")])
  size <- table(sp$pathway_id)
  sp <- sp[sp$pathway_id %in% names(size)[size >= 2 & size <= 1000], ]
  universe <- unique(sp$feature_id)
  rbind_list(lapply(names(feature_sets), function(set) {
    query <- intersect(feature_sets[[set]], universe)
    if (length(query) < 2) return(NULL)
    e <- rbind_list(lapply(split(sp, sp$pathway_id), function(t) {
      hits <- intersect(query, t$feature_id)
      data.frame(set = set, pathway_id = t$pathway_id[1], pathway = t$pathway[1],
                 gene_ratio = paste0(length(hits), "/", length(query)), count = length(hits),
                 p = stats::phyper(length(hits) - 1, nrow(t), length(universe) - nrow(t), length(query), lower.tail = FALSE),
                 kegg_ids = paste(hits, collapse = "/"))
    }))
    e$fdr <- stats::p.adjust(e$p, "BH")
    e
  }))
}

# KEGG pathway class and subclass.
kegg_hierarchy <- function() {
  lines <- readLines("https://rest.kegg.jp/get/br:br08901")
  level1 <- level2 <- NA_character_; out <- list()
  for (l in lines) {
    if (grepl("^A", l)) level1 <- gsub("<[^>]+>", "", sub("^A", "", l))
    else if (grepl("^B ", l)) level2 <- trimws(sub("^B", "", l))
    else if (grepl("^C ", l)) out[[length(out) + 1]] <- data.frame(map = sub("^C\\s+(\\d{5}).*$", "\\1", l),
                                                                  pathway_class = level1, pathway_subclass = level2)
  }
  rbind_list(out)
}

# One entry point per layer. Enriched: FDR-adjusted P < 0.05; nominal: P < 0.05 with at least 2 features.
kegg_ora_layer <- function(layer, feature_sets, universe, annotation, species_pathways = NULL) {
  e <- if (layer == "fecal_metagenome") {
    run_species_ora(feature_sets, universe, species_pathways)
  } else {
    type <- if (layer %in% c("pbmc_transcriptome", "serum_proteome")) "gene" else "compound"
    run_kegg_ora(feature_sets, map_to_kegg(layer, annotation), universe, kegg_terms(type))
  }
  if (is.null(e)) return(NULL)
  e$layer <- layer
  h <- kegg_hierarchy()
  e[c("pathway_class", "pathway_subclass")] <- h[match(sub("^[a-z]+", "", e$pathway_id), h$map), c("pathway_class", "pathway_subclass")]
  e$tier <- ifelse(e$fdr < 0.05, "Enriched", ifelse(e$p < 0.05 & e$count >= 2, "Nominal", "None"))
  e
}

# ---------------- Part 7: repeated-measures correlation networks ----------------
# Partial repeated-measures correlation: participant intercepts and visit means are removed from both features.
partial_rmcorr <- function(x, y, id, visit) {
  ok <- stats::complete.cases(x, y, id, visit)
  design <- stats::model.matrix(~ factor(id[ok]) + factor(visit[ok]))
  q <- qr(design)
  res <- qr.resid(q, cbind(x[ok], y[ok]))
  df <- sum(ok) - q$rank - 1
  r <- sum(res[, 1] * res[, 2]) / sqrt(sum(res[, 1]^2) * sum(res[, 2]^2))
  c(r = r, df = df, p = 2 * stats::pt(-abs(r * sqrt(df / (1 - r^2))), df))
}

# All pairwise partial correlations among the columns `nodes` of a wide table (pairwise complete observations).
correlate_pairs <- function(wide, nodes) {
  pairs <- utils::combn(nodes, 2, simplify = FALSE)
  chunks <- split(pairs, ceiling(seq_along(pairs) / 2000))
  res <- par_map(chunks, function(chunk, data) {
    t(vapply(chunk, function(p) partial_rmcorr(data[[p[1]]], data[[p[2]]], data$participant_id, data$visit), numeric(3)))
  }, data = wide, packages = character())
  out <- data.frame(node1 = unlist(lapply(pairs, `[`, 1)), node2 = unlist(lapply(pairs, `[`, 2)), do.call(rbind, res))
  out$fdr <- stats::p.adjust(out$p, "BH")
  out
}

# Join the baseline-standardized levels of selected features from several layers by participant and visit.
wide_feature_table <- function(features) {   # features: data.frame(layer, feature_id, node)
  tables <- lapply(split(features, features$layer), function(f) {
    prep <- load_layer(f$layer[1])
    d <- prep$metadata[, c("participant_id", "visit")]
    d[f$node] <- as.data.frame(prep$z[, f$feature_id, drop = FALSE])
    d
  })
  Reduce(function(a, b) merge(a, b, by = c("participant_id", "visit"), all = TRUE), tables)
}
