suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

root <- Sys.getenv("SCARVCELL_ROOT", unset = normalizePath(getwd(), winslash = "/", mustWork = TRUE))
stage <- file.path(Sys.getenv("SCARVCELL_STAGE", unset = file.path(root, ".stage")), "postlock_locked")
tables <- file.path(root, "results", "tables")
args <- commandArgs(trailingOnly = TRUE)
if (!length(args) || !args[1] %in% c("GSE282885", "GSE335482")) {
  stop("Provide GSE282885 or GSE335482")
}
dataset <- args[1]
set.seed(20260712)

annotation <- read.csv(file.path(tables, "scar_scrna_initial_cluster_annotation_marker_supported.csv"),
                       stringsAsFactors = FALSE)
cluster_to_type <- setNames(annotation$broad_cell_type, annotation$cluster)
compact_reference_file <- file.path(root, "data", "processed", "integrated",
                                    "scar_scrna_compact_label_transfer_reference.rds")
if (file.exists(compact_reference_file)) {
  reference <- readRDS(compact_reference_file)
} else {
  reference <- readRDS(file.path(root, "data", "processed", "integrated", "scar_scrna_initial_harmony.rds"))
  reference$broad_cell_type <- unname(cluster_to_type[as.character(reference$seurat_clusters)])
  reference <- subset(reference, cells = colnames(reference)[!is.na(reference$broad_cell_type)])
  reference_cells <- unlist(lapply(split(colnames(reference), reference$broad_cell_type), function(x) {
    sample(x, min(length(x), 800L))
  }), use.names = FALSE)
  reference <- subset(reference, cells = reference_cells)
  reference <- JoinLayers(reference, assay = "RNA")
  reference <- NormalizeData(reference, verbose = FALSE)
  reference <- FindVariableFeatures(reference, nfeatures = 2500, verbose = FALSE)
  reference <- ScaleData(reference, features = VariableFeatures(reference), verbose = FALSE)
  reference <- RunPCA(reference, features = VariableFeatures(reference), npcs = 30, verbose = FALSE)
  saveRDS(reference, compact_reference_file)
}

model <- read.csv(file.path(tables, "fssi_frozen_model.csv"), stringsAsFactors = FALSE)
rownames(model) <- model$gene

design <- if (dataset == "GSE282885") {
  data.frame(
    sample_id = paste0("GSM86521", 65:70),
    participant_id = paste0("unverified_", paste0("GSM86521", 65:70)),
    condition = c(rep("normal_skin", 3), rep("keloid", 3)),
    inference_unit = "library",
    stringsAsFactors = FALSE
  )
} else {
  data.frame(
    sample_id = paste0("GSM98144", 25:32),
    participant_id = rep(paste0("participant_", 1:4), each = 2),
    condition = rep(c("hypercellular_zone", "infiltrating_zone"), 4),
    inference_unit = "participant",
    stringsAsFactors = FALSE
  )
}

score_sample <- function(row) {
  sample_id <- row$sample_id
  count_dir <- file.path(stage, dataset, "matrices", sample_id)
  counts <- Read10X(data.dir = count_dir, gene.column = 2, unique.features = TRUE)
  if (is.list(counts)) counts <- counts[[1]]
  query <- CreateSeuratObject(counts = counts, project = dataset, min.cells = 3, min.features = 200)
  query[["percent.mt"]] <- PercentageFeatureSet(query, pattern = "^MT-")
  before <- ncol(query)
  query <- subset(query, subset = nFeature_RNA <= 7000 & percent.mt < 20)
  after <- ncol(query)
  if (after < 20) {
    audit <- data.frame(dataset_id = dataset, sample_id = sample_id,
      participant_id = row$participant_id, condition = row$condition,
      inference_unit = row$inference_unit, cells_after_min_features = before,
      cells_after_qc = after, predicted_fibroblasts = 0, genes_present = NA_integer_,
      status = "insufficient_cells_after_fixed_qc")
    return(list(scores = NULL, audit = audit))
  }
  query <- NormalizeData(query, verbose = FALSE)
  anchors <- FindTransferAnchors(reference = reference, query = query, reference.reduction = "pca",
                                 dims = 1:30, features = VariableFeatures(reference), verbose = FALSE)
  prediction <- TransferData(anchorset = anchors, refdata = reference$broad_cell_type,
                             dims = 1:30, verbose = FALSE)
  query <- AddMetaData(query, prediction)
  fibro <- colnames(query)[query$predicted.id == "fibroblast" & query$prediction.score.max >= 0.5]
  if (!length(fibro)) {
    audit <- data.frame(dataset_id = dataset, sample_id = sample_id,
      participant_id = row$participant_id, condition = row$condition,
      inference_unit = row$inference_unit, cells_after_min_features = before,
      cells_after_qc = after, predicted_fibroblasts = 0, genes_present = NA_integer_,
      status = "no_fibroblasts_at_locked_threshold")
    return(list(scores = NULL, audit = audit))
  }
  norm <- LayerData(query, assay = "RNA", layer = "data")[, fibro, drop = FALSE]
  genes <- intersect(model$gene, rownames(norm))
  m <- model[genes, , drop = FALSE]
  z <- sweep(as.matrix(norm[genes, , drop = FALSE]), 1, m$training_mean, "-")
  z <- sweep(z, 1, m$training_sd, "/")
  w <- m$weight / sum(abs(m$weight))
  ew <- m$equal_weight / sum(abs(m$equal_weight))
  scores <- data.frame(
    dataset_id = dataset, sample_id = sample_id, participant_id = row$participant_id,
    condition = row$condition, inference_unit = row$inference_unit, cell_id = fibro,
    prediction_score = query$prediction.score.max[fibro],
    fssi_frozen_weighted = as.numeric(crossprod(w, z)),
    fssi_frozen_equal_weight = as.numeric(crossprod(ew, z)),
    pathological_component = colMeans(z[m$program == "pathological", , drop = FALSE]),
    repair_component_reversed = -colMeans(z[m$program == "repair", , drop = FALSE]),
    genes_present = length(genes), stringsAsFactors = FALSE
  )
  fib5 <- intersect(c("ADAM12", "COMP", "POSTN"), rownames(z))
  scores$Fib5_three_marker_baseline <- colMeans(z[fib5, , drop = FALSE])
  gene_z <- matrix(NA_real_, nrow = ncol(z), ncol = nrow(model),
                   dimnames = list(colnames(z), paste0("frozen_z__", model$gene)))
  gene_z[, paste0("frozen_z__", rownames(z))] <- t(z)
  scores <- cbind(scores, as.data.frame(gene_z, stringsAsFactors = FALSE))
  audit <- data.frame(dataset_id = dataset, sample_id = sample_id,
    participant_id = row$participant_id, condition = row$condition,
    inference_unit = row$inference_unit, cells_after_min_features = before,
    cells_after_qc = after, predicted_fibroblasts = length(fibro),
    genes_present = length(genes), status = "projected_at_locked_threshold")
  rm(counts, query, norm, z)
  gc(verbose = FALSE)
  list(scores = scores, audit = audit)
}

outputs <- lapply(seq_len(nrow(design)), function(i) score_sample(design[i, ]))
audits <- do.call(rbind, lapply(outputs, `[[`, "audit"))
score_list <- lapply(outputs, `[[`, "scores")
score_list <- score_list[!vapply(score_list, is.null, logical(1))]
write.csv(audits, file.path(tables, paste0("postlock_", dataset, "_projection_audit.csv")), row.names = FALSE)
if (!length(score_list)) stop("No samples yielded fibroblasts at the locked threshold")
scores <- do.call(rbind, score_list)
write.csv(scores, gzfile(file.path(tables, paste0("postlock_", dataset, "_frozen_cell_scores.csv.gz"))),
          row.names = FALSE)
numeric <- c("fssi_frozen_weighted", "fssi_frozen_equal_weight", "pathological_component",
             "repair_component_reversed", "Fib5_three_marker_baseline",
             paste0("frozen_z__", model$gene))
unit <- aggregate(scores[, numeric, drop = FALSE],
                  by = scores[, c("dataset_id", "sample_id", "participant_id", "condition", "inference_unit"), drop = FALSE],
                  FUN = mean, na.rm = TRUE)
cell_n <- aggregate(cell_id ~ dataset_id + sample_id + participant_id + condition + inference_unit,
                    data = scores, FUN = length)
names(cell_n)[6] <- "fibroblast_cell_n"
unit <- merge(unit, cell_n,
              by = c("dataset_id", "sample_id", "participant_id", "condition", "inference_unit"), all.x = TRUE)
write.csv(unit, file.path(tables, paste0("postlock_", dataset, "_frozen_unit_scores.csv")), row.names = FALSE)
print(audits)
print(unit[, c("sample_id", "participant_id", "condition", "fssi_frozen_weighted", "fibroblast_cell_n")])
