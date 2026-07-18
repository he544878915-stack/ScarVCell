suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(SingleCellExperiment)
  library(scDblFinder)
})

root <- Sys.getenv("SCARVCELL_ROOT", unset = normalizePath(getwd(), winslash = "/", mustWork = TRUE))
stage <- file.path(Sys.getenv("SCARVCELL_STAGE", unset = file.path(root, ".stage")), "external", "GSE191067")
tables <- file.path(root, "results", "tables")
set.seed(20260712)

reference <- readRDS(file.path(root, "data", "processed", "integrated",
                               "scar_scrna_compact_label_transfer_reference.rds"))
matrix_dir <- file.path(stage, "locked_feature_matrix")
counts <- Matrix::readMM(gzfile(file.path(matrix_dir, "matrix.mtx.gz")))
feature_table <- read.delim(gzfile(file.path(matrix_dir, "features.tsv.gz")), header = FALSE, stringsAsFactors = FALSE)
barcode_table <- read.delim(gzfile(file.path(matrix_dir, "barcodes.tsv.gz")), header = FALSE, stringsAsFactors = FALSE)
rownames(counts) <- make.unique(feature_table[[2]])
colnames(counts) <- barcode_table[[1]]
counts <- as(counts, "dgCMatrix")
markers <- read.csv(file.path(stage, "locked_marker_sets.csv"), stringsAsFactors = FALSE)
sample_id <- sub("_.*$", "", colnames(counts))

marker_sets <- list(
  fibroblast = markers$gene[markers$set == "fibroblast"],
  endothelial = markers$gene[markers$set == "endothelial"],
  epithelial = markers$gene[markers$set == "epithelial"],
  neural = markers$gene[markers$set == "neural"],
  myeloid = markers$gene[markers$set == "myeloid"],
  lymphoid = c("CD3D", "CD3E", "CD79A", "MS4A1", "NKG7")
)

collapse_transfer <- function(x) {
  out <- rep("other", length(x))
  out[x == "fibroblast"] <- "fibroblast"
  out[x %in% c("endothelial", "lymphatic_endothelial")] <- "endothelial"
  out[x == "keratinocyte"] <- "epithelial"
  out[x == "neural_glial"] <- "neural"
  out[x == "myeloid"] <- "myeloid"
  out[x %in% c("T_NK", "mast")] <- "lymphoid"
  out
}

all_meta <- list()
for (sid in unique(sample_id)) {
  message("QC sample ", sid)
  cells <- colnames(counts)[sample_id == sid]
  sample_counts <- counts[, cells, drop = FALSE]
  query <- CreateSeuratObject(counts = sample_counts, project = "GSE191067")
  query <- NormalizeData(query, verbose = FALSE)
  transfer_features <- intersect(VariableFeatures(reference), rownames(query))
  anchors <- FindTransferAnchors(reference = reference, query = query, reference.reduction = "pca",
                                 dims = 1:30, features = transfer_features, verbose = FALSE)
  prediction <- TransferData(anchorset = anchors, refdata = reference$broad_cell_type,
                             dims = 1:30, verbose = FALSE)
  norm <- LayerData(query, assay = "RNA", layer = "data")
  score_matrix <- sapply(marker_sets, function(genes) {
    genes <- intersect(genes, rownames(norm))
    if (!length(genes)) return(rep(0, ncol(norm)))
    Matrix::colMeans(norm[genes, , drop = FALSE])
  })
  marker_label <- colnames(score_matrix)[max.col(score_matrix, ties.method = "first")]
  marker_margin <- apply(score_matrix, 1, function(x) sort(x, decreasing = TRUE)[1] - sort(x, decreasing = TRUE)[2])

  sce <- SingleCellExperiment(list(counts = sample_counts))
  sce <- scDblFinder(sce, samples = rep(sid, ncol(sce)), verbose = FALSE)
  meta <- data.frame(
    cell_id = cells, sample_id = sid,
    transfer_label = prediction$predicted.id,
    transfer_score = prediction$prediction.score.max,
    transfer_collapsed = collapse_transfer(prediction$predicted.id),
    marker_label = marker_label,
    marker_margin = marker_margin,
    scDblFinder_score = colData(sce)$scDblFinder.score,
    scDblFinder_class = colData(sce)$scDblFinder.class,
    stringsAsFactors = FALSE
  )
  all_meta[[sid]] <- meta
  rm(query, anchors, prediction, norm, score_matrix, sce, sample_counts)
  gc(verbose = FALSE)
}

meta <- do.call(rbind, all_meta)
rownames(meta) <- NULL
write.csv(meta, gzfile(file.path(tables, "GSE191067_annotation_doublet_cell_audit.csv.gz")), row.names = FALSE)

conf <- as.data.frame.matrix(table(meta$transfer_collapsed, meta$marker_label))
conf$transfer_collapsed <- rownames(conf)
rownames(conf) <- NULL
write.csv(conf, file.path(tables, "GSE191067_transfer_marker_confusion.csv"), row.names = FALSE)

summary <- aggregate(cbind(doublet = meta$scDblFinder_class == "doublet",
                           marker_concordant = meta$transfer_collapsed == meta$marker_label),
                     by = list(sample_id = meta$sample_id, transfer_collapsed = meta$transfer_collapsed), mean)
counts_summary <- aggregate(meta$cell_id, by = list(sample_id = meta$sample_id,
                                                    transfer_collapsed = meta$transfer_collapsed), length)
names(counts_summary)[3] <- "cell_n"
summary <- merge(summary, counts_summary, by = c("sample_id", "transfer_collapsed"))
write.csv(summary, file.path(tables, "GSE191067_annotation_doublet_summary.csv"), row.names = FALSE)

fibro_scores <- read.csv(gzfile(file.path(tables, "GSE191067_frozen_fibroblast_cell_scores.csv.gz")),
                         stringsAsFactors = FALSE)
fibro_scores <- merge(fibro_scores, meta[, c("cell_id", "scDblFinder_class")], by = "cell_id", all.x = TRUE)
fibro_scores <- fibro_scores[fibro_scores$prediction_score >= 0.5 & fibro_scores$scDblFinder_class == "singlet", ]
score_columns <- c("fssi_frozen_weighted", "fssi_frozen_equal_weight", "pathological_component",
                   "repair_component_reversed", "Fib5_three_marker_baseline")
singlet_units <- aggregate(fibro_scores[, score_columns],
                           by = fibro_scores[, c("sample_id", "condition")], mean)
singlet_units$fibroblast_singlet_n <- as.integer(table(fibro_scores$sample_id)[singlet_units$sample_id])
write.csv(singlet_units, file.path(tables, "GSE191067_singlet_only_sample_scores.csv"), row.names = FALSE)

audit <- data.frame(
  item = c("full_transcriptome_de_novo_annotation", "independent_marker_concordance",
           "locked_feature_scDblFinder_sensitivity", "ambient_RNA_from_empty_droplets",
           "author_label_confusion_matrix"),
  status = c("not_possible_from_released_locked_workflow_without_materially_expanding_matrix",
             "completed", "completed", "not_evaluable", "not_evaluable"),
  reason = c("public combined matrix can be streamed, but the frozen workflow retained 2435 predeclared features; no full Seurat object was released",
             "six independent broad marker modules compared with frozen reference transfer",
             "scDblFinder applied per matrix sample to the 2435-feature count matrix; sensitivity only",
             "raw unfiltered droplet matrix and empty droplets were not released",
             "cell-level author annotations or labelled Seurat object were not released"),
  stringsAsFactors = FALSE
)
write.csv(audit, file.path(tables, "GSE191067_QC_capability_audit.csv"), row.names = FALSE)
cat("Completed GSE191067 annotation and doublet sensitivity for", nrow(meta), "cells\n")
