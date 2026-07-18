suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

root <- Sys.getenv("SCARVCELL_ROOT", unset = normalizePath(getwd(), winslash = "/", mustWork = TRUE))
stage <- file.path(Sys.getenv("SCARVCELL_STAGE", unset = file.path(root, ".stage")), "external", "GSE191067")
tables <- file.path(root, "results", "tables")
set.seed(20260712)

reference <- readRDS(file.path(root, "data", "processed", "integrated",
                               "scar_scrna_compact_label_transfer_reference.rds"))
model <- read.csv(file.path(tables, "fssi_frozen_model.csv"), stringsAsFactors = FALSE)
rownames(model) <- model$gene
markers <- read.csv(file.path(stage, "locked_marker_sets.csv"), stringsAsFactors = FALSE)
matrix_dir <- file.path(stage, "locked_feature_matrix")
counts <- Matrix::readMM(file.path(matrix_dir, "matrix.mtx"))
feature_table <- read.delim(file.path(matrix_dir, "features.tsv"), header = FALSE, stringsAsFactors = FALSE)
barcode_table <- read.delim(file.path(matrix_dir, "barcodes.tsv"), header = FALSE, stringsAsFactors = FALSE)
rownames(counts) <- make.unique(feature_table[[2]])
colnames(counts) <- barcode_table[[1]]
counts <- as(counts, "dgCMatrix")

sample_id <- sub("_.*$", "", colnames(counts))
condition_map <- c(HK1 = "keloid", HK2 = "keloid", HK3 = "keloid",
                   `HK-NS1` = "perilesional_skin", `HK-NS2` = "perilesional_skin", `HK-NS3` = "perilesional_skin",
                   HNS1 = "normal_skin", HNS2 = "normal_skin", HNS3 = "normal_skin",
                   HNSR1 = "normal_scar", HNSR2 = "normal_scar", HNSR3 = "normal_scar")

score_cells <- list()
all_meta <- list()
all_expr <- list()
audit <- list()
for (sid in unique(sample_id)) {
  cells <- colnames(counts)[sample_id == sid]
  query <- CreateSeuratObject(counts = counts[, cells, drop = FALSE], project = "GSE191067")
  query <- NormalizeData(query, verbose = FALSE)
  transfer_features <- intersect(VariableFeatures(reference), rownames(query))
  anchors <- FindTransferAnchors(reference = reference, query = query, reference.reduction = "pca",
                                 dims = 1:30, features = transfer_features, verbose = FALSE)
  prediction <- TransferData(anchorset = anchors, refdata = reference$broad_cell_type,
                             dims = 1:30, verbose = FALSE)
  query <- AddMetaData(query, prediction)
  norm <- LayerData(query, assay = "RNA", layer = "data")
  meta <- data.frame(cell_id = colnames(query), sample_id = sid,
                     condition = unname(condition_map[sid]), predicted_id = query$predicted.id,
                     prediction_score = query$prediction.score.max, stringsAsFactors = FALSE)

  expressed_count <- function(set_name) {
    genes <- intersect(markers$gene[markers$set == set_name], rownames(norm))
    if (!length(genes)) return(rep(0L, ncol(norm)))
    Matrix::colSums(norm[genes, , drop = FALSE] > 0)
  }
  meta$fibroblast_marker_n <- expressed_count("fibroblast")
  meta$immune_marker_n <- expressed_count("immune")
  meta$endothelial_marker_n <- expressed_count("endothelial")
  meta$epithelial_marker_n <- expressed_count("epithelial")
  meta$neural_marker_n <- expressed_count("neural")
  meta$marker_pure_fibroblast <- meta$fibroblast_marker_n >= 3 & meta$immune_marker_n <= 1 &
    meta$endothelial_marker_n <= 1 & meta$epithelial_marker_n <= 1 & meta$neural_marker_n <= 1

  fibro <- meta$cell_id[meta$predicted_id == "fibroblast" & meta$prediction_score >= 0.4]
  if (length(fibro)) {
    genes <- intersect(model$gene, rownames(norm))
    m <- model[genes, , drop = FALSE]
    z <- sweep(as.matrix(norm[genes, fibro, drop = FALSE]), 1, m$training_mean, "-")
    z <- sweep(z, 1, m$training_sd, "/")
    w <- m$weight / sum(abs(m$weight))
    ew <- m$equal_weight / sum(abs(m$equal_weight))
    fmeta <- meta[match(fibro, meta$cell_id), , drop = FALSE]
    fmeta$fssi_frozen_weighted <- as.numeric(crossprod(w, z))
    fmeta$fssi_frozen_equal_weight <- as.numeric(crossprod(ew, z))
    fmeta$pathological_component <- colMeans(z[m$program == "pathological", , drop = FALSE])
    fmeta$repair_component_reversed <- -colMeans(z[m$program == "repair", , drop = FALSE])
    fib5 <- intersect(c("ADAM12", "COMP", "POSTN"), rownames(z))
    fmeta$Fib5_three_marker_baseline <- colMeans(z[fib5, , drop = FALSE])
    ecm <- intersect(markers$gene[markers$set == "collagen_ecm"], rownames(norm))
    fmeta$collagen_ECM_raw <- Matrix::colMeans(norm[ecm, fibro, drop = FALSE])
    score_cells[[sid]] <- fmeta
  }
  all_meta[[sid]] <- meta
  all_expr[[sid]] <- norm
  audit[[sid]] <- data.frame(
    dataset_id = "GSE191067", sample_id = sid, condition = unname(condition_map[sid]),
    cells_after_source_qc = ncol(query), transfer_features = length(transfer_features),
    fibroblasts_threshold_0.4 = sum(meta$predicted_id == "fibroblast" & meta$prediction_score >= 0.4),
    fibroblasts_threshold_0.5 = sum(meta$predicted_id == "fibroblast" & meta$prediction_score >= 0.5),
    fibroblasts_threshold_0.6 = sum(meta$predicted_id == "fibroblast" & meta$prediction_score >= 0.6),
    myeloid_threshold_0.5 = sum(meta$predicted_id == "myeloid" & meta$prediction_score >= 0.5),
    stringsAsFactors = FALSE
  )
  rm(query, prediction, anchors, z)
  gc(verbose = FALSE)
}

scores <- do.call(rbind, score_cells)
rownames(scores) <- NULL
write.csv(scores, gzfile(file.path(tables, "GSE191067_frozen_fibroblast_cell_scores.csv.gz")), row.names = FALSE)
write.csv(do.call(rbind, audit), file.path(tables, "GSE191067_projection_audit.csv"), row.names = FALSE)

unit_rows <- list()
numeric_scores <- c("fssi_frozen_weighted", "fssi_frozen_equal_weight", "pathological_component",
                    "repair_component_reversed", "Fib5_three_marker_baseline", "collagen_ECM_raw")
for (threshold in c(0.4, 0.5, 0.6)) {
  for (purity in c("reference_transfer", "marker_purity_subset")) {
    keep <- scores$prediction_score >= threshold
    if (purity == "marker_purity_subset") keep <- keep & scores$marker_pure_fibroblast
    sub <- scores[keep, , drop = FALSE]
    if (!nrow(sub)) next
    means <- aggregate(sub[, numeric_scores, drop = FALSE],
                       by = sub[, c("sample_id", "condition"), drop = FALSE], FUN = mean)
    means$threshold <- threshold
    means$annotation_sensitivity <- purity
    means$fibroblast_cell_n <- as.integer(table(sub$sample_id)[means$sample_id])
    unit_rows[[length(unit_rows) + 1]] <- means
  }
}
units <- do.call(rbind, unit_rows)
units$dataset_id <- "GSE191067"
write.csv(units, file.path(tables, "GSE191067_frozen_sample_scores_sensitivity.csv"), row.names = FALSE)

# Target-blind focused communication in the new cohort.
prior <- readRDS(file.path(root, "data", "external", "nichenet", "ligand_target_matrix_nsga2r_final.rds"))
axes <- data.frame(
  axis_id = c("POSTN_ITGAV_ITGB1", "TGFB1_TGFBR", "CXCL12_CXCR4", "CSF1_CSF1R", "PDGFB_PDGFRA", "MIF_CD74"),
  ligand = c("POSTN", "TGFB1", "CXCL12", "CSF1", "PDGFB", "MIF"),
  receptor = c("ITGAV;ITGB1", "TGFBR1;TGFBR2", "CXCR4", "CSF1R", "PDGFRA", "CD74"),
  direction = c("fibro_to_myeloid", "fibro_to_myeloid", "fibro_to_myeloid", "fibro_to_myeloid", "myeloid_to_fibro", "fibro_to_myeloid"),
  expected_delta = c("positive", "positive", "negative", "negative", "negative", "positive"),
  stringsAsFactors = FALSE
)
gene_stats <- function(expr, genes, cells) {
  genes <- intersect(strsplit(genes, ";", fixed = TRUE)[[1]], rownames(expr))
  cells <- intersect(cells, colnames(expr))
  if (!length(genes) || !length(cells)) return(c(avg = NA_real_, pct = NA_real_))
  av <- Matrix::rowMeans(expr[genes, cells, drop = FALSE])
  pct <- Matrix::rowMeans(expr[genes, cells, drop = FALSE] > 0)
  c(avg = min(av), pct = min(pct))
}
lr_score <- function(expr, ligand, receptor, source, target) {
  lig <- gene_stats(expr, ligand, source); rec <- gene_stats(expr, receptor, target)
  lig[["avg"]] * rec[["avg"]] * sqrt(lig[["pct"]] * rec[["pct"]])
}
comm_rows <- list()
for (sid in names(all_meta)) {
  meta <- all_meta[[sid]]; expr <- all_expr[[sid]]
  fibro <- scores[scores$sample_id == sid & scores$prediction_score >= 0.5, , drop = FALSE]
  myeloid <- meta$cell_id[meta$predicted_id == "myeloid" & meta$prediction_score >= 0.5]
  for (i in seq_len(nrow(axes))) {
    ax <- axes[i, ]
    top_targets <- character()
    if (ax$ligand %in% colnames(prior)) {
      values <- prior[, ax$ligand]
      top_targets <- names(sort(values, decreasing = TRUE))[seq_len(min(50L, length(values)))]
    }
    removed <- unique(c(ax$ligand, strsplit(ax$receptor, ";", fixed = TRUE)[[1]], top_targets))
    path_genes <- setdiff(model$gene[model$program == "pathological"], removed)
    repair_genes <- setdiff(model$gene[model$program == "repair"], removed)
    eligible <- intersect(fibro$cell_id, colnames(expr))
    if (!length(eligible) || !length(path_genes) || !length(repair_genes)) next
    state <- Matrix::colMeans(expr[path_genes, eligible, drop = FALSE]) -
      Matrix::colMeans(expr[repair_genes, eligible, drop = FALSE])
    q <- quantile(state, c(0.25, 0.75), na.rm = TRUE)
    low <- eligible[state <= q[[1]]]; high <- eligible[state >= q[[2]]]
    if (min(length(low), length(high), length(myeloid)) < 10) next
    if (ax$direction == "fibro_to_myeloid") {
      high_score <- lr_score(expr, ax$ligand, ax$receptor, high, myeloid)
      low_score <- lr_score(expr, ax$ligand, ax$receptor, low, myeloid)
    } else {
      high_score <- lr_score(expr, ax$ligand, ax$receptor, myeloid, high)
      low_score <- lr_score(expr, ax$ligand, ax$receptor, myeloid, low)
    }
    delta <- high_score - low_score
    comm_rows[[length(comm_rows) + 1]] <- data.frame(
      dataset_id = "GSE191067", sample_id = sid, condition = unname(condition_map[sid]),
      axis_id = ax$axis_id, ligand = ax$ligand, receptor = ax$receptor,
      expected_delta = ax$expected_delta, high_fibro_cells = length(high), low_fibro_cells = length(low),
      myeloid_cells = length(myeloid), downstream_target_n_removed = length(top_targets),
      high_score = high_score, low_score = low_score, delta_high_minus_low = delta,
      expected_direction_support = ifelse(ax$expected_delta == "positive", delta > 0, delta < 0),
      stringsAsFactors = FALSE)
  }
}
communication <- do.call(rbind, comm_rows)
write.csv(communication, file.path(tables, "GSE191067_target_blind_communication_by_sample.csv"), row.names = FALSE)
cat("Projected", nrow(scores), "fibroblasts at threshold >=0.4 across", length(unique(scores$sample_id)), "samples\n")
cat("Communication rows:", nrow(communication), "\n")
