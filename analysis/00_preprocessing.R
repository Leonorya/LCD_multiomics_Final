# Preprocessing of the six omics layers (Methods: omics data preprocessing; Figure S1).
# Input : data/clinical.csv, data/omics/<layer>.csv, data/omics/<layer>_annotation.csv, data/sample_lists/*.txt,
#         data/qc/<layer>_with_pooled_qc.csv (serum proteome, serum metabolome and fecal metabolome)
# Output: results/preprocessed/<layer>.rds, results/preprocessed/qc_pca.rds
source("analysis/utils.R")

# ---------------- Part 1: settings ----------------
MIN_PREVALENCE <- 0.10   # microbial species present in fewer samples are removed before CLR
ZERO_REPLACE   <- 1e-8   # replacement value for zero abundances before CLR
SAMPLE_SET <- c(pbmc_transcriptome = "blood", serum_proteome = "blood", serum_metabolome = "blood",
                fecal_metagenome = "stool", fecal_metatranscriptome = "stool", fecal_metabolome = "stool")
QC_LAYERS <- c("serum_proteome", "serum_metabolome", "fecal_metabolome")

clinical <- read.csv(file.path(DATA_DIR, "clinical.csv"))
clinical$visit <- as.character(clinical$visit)

# ---------------- Part 2: transformation and baseline standardization ----------------
for (layer in LAYERS) {
  omics <- read.csv(file.path(DATA_DIR, "omics", paste0(layer, ".csv")), check.names = FALSE)
  annotation <- read.csv(file.path(DATA_DIR, "omics", paste0(layer, "_annotation.csv")))
  included <- readLines(file.path(DATA_DIR, "sample_lists", paste0(SAMPLE_SET[[layer]], "_samples.txt")))
  omics <- omics[omics$sample_id %in% included & omics$sample_id %in% clinical$sample_id, ]
  meta <- clinical[match(omics$sample_id, clinical$sample_id), ]

  x <- as.matrix(omics[, -1])
  rownames(x) <- omics$sample_id
  if (LAYER_TRANSFORM[[layer]] == "clr") {
    x <- x[, colMeans(x > 0, na.rm = TRUE) >= MIN_PREVALENCE, drop = FALSE]
    x[x == 0] <- ZERO_REPLACE
  }
  transformed <- transform_layer(x, LAYER_TRANSFORM[[layer]])
  scaled <- baseline_zscore(transformed, meta$visit)   # features with zero variance at baseline are excluded
  kept <- colnames(scaled$z)

  saveRDS(list(layer = layer, metadata = meta, annotation = annotation[match(kept, annotation$feature_id), ],
               transformed = transformed[, kept, drop = FALSE], z = scaled$z),
          result_path("preprocessed", paste0(layer, ".rds")))
}

# ---------------- Part 3: principal component analysis of study and pooled QC samples (Figure S1) ----------------
qc_pca <- rbind_list(lapply(QC_LAYERS, function(layer) {
  d <- read.csv(file.path(DATA_DIR, "qc", paste0(layer, "_with_pooled_qc.csv")), check.names = FALSE)
  x <- transform_layer(as.matrix(d[, -(1:2)]), if (layer == "fecal_metabolome") "none" else "log2")
  x <- x[, colSums(!is.finite(x)) == 0 & apply(x, 2, stats::sd) > 0, drop = FALSE]
  pca <- stats::prcomp(x, center = TRUE, scale. = TRUE)
  meta <- load_layer(layer)$metadata
  data.frame(layer = layer, sample_id = d$sample_id, sample_type = d$sample_type,
             participant_id = meta$participant_id[match(d$sample_id, meta$sample_id)],
             pc1 = pca$x[, 1], pc2 = pca$x[, 2],
             pc1_var = summary(pca)$importance[2, 1], pc2_var = summary(pca)$importance[2, 2])
}))
saveRDS(qc_pca, result_path("preprocessed", "qc_pca.rds"))
