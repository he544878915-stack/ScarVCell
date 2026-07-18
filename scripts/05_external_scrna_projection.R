suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

root <- Sys.getenv("SCARVCELL_ROOT", unset = normalizePath(getwd(), winslash = "/", mustWork = TRUE))
stage <- file.path(Sys.getenv("SCARVCELL_STAGE", unset = file.path(root, ".stage")), "external_scrna")
tables <- file.path(root, "results", "tables")
args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("Provide one dataset accession")
dataset <- args[1]
set.seed(20260711)

annotation <- read.csv(file.path(tables, "scar_scrna_initial_cluster_annotation_marker_supported.csv"),
                       stringsAsFactors = FALSE)
cluster_to_type <- setNames(annotation$broad_cell_type, annotation$cluster)
compact_reference_file <- file.path(root, "data", "processed", "integrated", "scar_scrna_compact_label_transfer_reference.rds")
if (file.exists(compact_reference_file)) {
  reference <- readRDS(compact_reference_file)
} else {
  reference <- readRDS(file.path(root, "data", "processed", "integrated", "scar_scrna_initial_harmony.rds"))
  reference$broad_cell_type <- unname(cluster_to_type[as.character(reference$seurat_clusters)])
  reference <- subset(reference, cells = colnames(reference)[!is.na(reference$broad_cell_type)])
  # The compact reference is sampled without looking at any external expression.
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

read_triplet <- function(matrix_file, feature_file, barcode_file) {
  ReadMtx(mtx = matrix_file, features = feature_file, cells = barcode_file,
          feature.column = 2, cell.column = 1, unique.features = TRUE)
}

score_query <- function(counts, sample_id, donor_id, condition) {
  query <- CreateSeuratObject(counts = counts, project = dataset, min.cells = 3, min.features = 200)
  query[["percent.mt"]] <- PercentageFeatureSet(query, pattern = "^MT-")
  before <- ncol(query)
  query <- subset(query, subset = nFeature_RNA <= 7000 & percent.mt < 20)
  after <- ncol(query)
  if (after < 20) {
    return(list(scores = NULL, audit = data.frame(dataset_id = dataset, sample_id = sample_id,
      donor_id = donor_id, condition = condition, cells_after_min_features = before,
      cells_after_qc = after, predicted_fibroblasts = 0, status = "insufficient_cells_after_fixed_qc")))
  }
  query <- NormalizeData(query, verbose = FALSE)
  anchors <- FindTransferAnchors(reference = reference, query = query, reference.reduction = "pca",
                                 dims = 1:30, features = VariableFeatures(reference), verbose = FALSE)
  prediction <- TransferData(anchorset = anchors, refdata = reference$broad_cell_type,
                             dims = 1:30, verbose = FALSE)
  query <- AddMetaData(query, prediction)
  fibro <- colnames(query)[query$predicted.id == "fibroblast" & query$prediction.score.max >= 0.5]
  if (!length(fibro)) {
    return(list(scores = NULL, audit = data.frame(dataset_id = dataset, sample_id = sample_id,
      donor_id = donor_id, condition = condition, cells_after_min_features = before,
      cells_after_qc = after, predicted_fibroblasts = 0, status = "no_fibroblasts_at_locked_threshold")))
  }
  data <- LayerData(query, assay = "RNA", layer = "data")[, fibro, drop = FALSE]
  genes <- intersect(model$gene, rownames(data))
  m <- model[genes, , drop = FALSE]
  z <- sweep(as.matrix(data[genes, , drop = FALSE]), 1, m$training_mean, "-")
  z <- sweep(z, 1, m$training_sd, "/")
  w <- m$weight / sum(abs(m$weight))
  ew <- m$equal_weight / sum(abs(m$equal_weight))
  score <- data.frame(dataset_id = dataset, sample_id = sample_id, donor_id = donor_id,
                      condition = condition, cell_id = fibro,
                      prediction_score = query$prediction.score.max[fibro],
                      fssi_frozen_weighted = as.numeric(crossprod(w, z)),
                      fssi_frozen_equal_weight = as.numeric(crossprod(ew, z)),
                      pathological_component = colMeans(z[m$program == "pathological", , drop = FALSE]),
                      repair_component_reversed = -colMeans(z[m$program == "repair", , drop = FALSE]),
                      genes_present = length(genes), stringsAsFactors = FALSE)
  fib5 <- intersect(c("ADAM12", "COMP", "POSTN"), rownames(z))
  score$Fib5_three_marker_baseline <- colMeans(z[fib5, , drop = FALSE])
  gene_z <- matrix(NA_real_, nrow = ncol(z), ncol = nrow(model),
                   dimnames = list(colnames(z), paste0("frozen_z__", model$gene)))
  gene_z[, paste0("frozen_z__", rownames(z))] <- t(z)
  gene_z <- as.data.frame(gene_z, stringsAsFactors = FALSE)
  score <- cbind(score, gene_z)
  audit <- data.frame(dataset_id = dataset, sample_id = sample_id, donor_id = donor_id,
                      condition = condition, cells_after_min_features = before,
                      cells_after_qc = after, predicted_fibroblasts = length(fibro),
                      status = "projected_at_locked_threshold", stringsAsFactors = FALSE)
  rm(query, data, z)
  gc(verbose = FALSE)
  list(scores = score, audit = audit)
}

file_triplets <- function(directory) {
  matrices <- list.files(directory, pattern = "matrix\\.mtx\\.gz$", full.names = TRUE)
  lapply(matrices, function(mtx) {
    prefix <- sub("_matrix\\.mtx\\.gz$", "", mtx)
    feature <- c(paste0(prefix, "_features.tsv.gz"), paste0(prefix, "_genes.tsv.gz"))
    feature <- feature[file.exists(feature)][1]
    list(matrix = mtx, feature = feature, barcode = paste0(prefix, "_barcodes.tsv.gz"),
         prefix = basename(prefix))
  })
}

sample_design <- function(accession, prefix) {
  if (accession == "GSE156326") {
    gsm <- sub("_.*", "", prefix)
    condition <- if (grepl("human_scar", prefix)) "hypertrophic_scar" else "normal_skin"
    return(c(sample_id = gsm, donor_id = gsm, condition = condition))
  }
  if (accession == "GSE181316") {
    gsm <- sub("_.*", "", prefix)
    condition <- if (grepl("keloid", prefix)) "keloid" else if (grepl("scar", prefix)) "normal_scar" else "healthy_skin"
    donor <- if (grepl("keloid_3[LR]", prefix)) "keloid_donor_3_collapsed" else gsm
    return(c(sample_id = gsm, donor_id = donor, condition = condition))
  }
  if (accession == "GSE151177") {
    gsm <- sub("_.*", "", prefix)
    condition <- if (grepl("Control", prefix, ignore.case = TRUE)) "healthy_control" else "psoriasis"
    return(c(sample_id = gsm, donor_id = gsm, condition = condition))
  }
  stop("Unsupported triplet dataset")
}

results <- list()
if (dataset %in% c("GSE156326", "GSE181316", "GSE151177")) {
  dataset_dir <- if (dataset == "GSE181316") paste0(dataset, "_filtered") else dataset
  triplets <- file_triplets(file.path(stage, dataset_dir))
  if (dataset == "GSE156326") {
    triplets <- triplets[vapply(triplets, function(x) grepl("_human_", x$prefix), logical(1))]
  }
  if (dataset == "GSE151177") {
    triplets <- triplets[vapply(triplets, function(x) grepl("_Control0[1-5]$", x$prefix), logical(1))]
  }
  for (i in seq_along(triplets)) {
    triplet <- triplets[[i]]
    design <- sample_design(dataset, triplet$prefix)
    counts <- read_triplet(triplet$matrix, triplet$feature, triplet$barcode)
    results[[i]] <- score_query(counts, design[["sample_id"]], design[["donor_id"]], design[["condition"]])
    rm(counts)
    gc(verbose = FALSE)
  }
} else if (dataset == "GSE130973") {
  counts <- read_triplet(file.path(stage, "GSE130973_matrix_filtered_full.mtx.gz"),
                         file.path(stage, "GSE130973_genes_filtered.tsv.gz"),
                         file.path(stage, "GSE130973_barcodes_filtered.tsv.gz"))
  suffix <- sub(".*-", "", colnames(counts))
  for (i in sort(unique(suffix))) {
    cells <- colnames(counts)[suffix == i]
    results[[length(results) + 1L]] <- score_query(counts[, cells, drop = FALSE],
      paste0("GSE130973_donor_", i), paste0("GSE130973_donor_", i), "healthy_skin_reference")
  }
} else stop("Unsupported dataset")

score_list <- lapply(results, `[[`, "scores")
score_list <- score_list[!vapply(score_list, is.null, logical(1))]
audits <- do.call(rbind, lapply(results, `[[`, "audit"))
write.csv(audits, file.path(tables, paste0("external_scrna_", dataset, "_projection_audit.csv")), row.names = FALSE)
if (length(score_list)) {
  scores <- do.call(rbind, score_list)
  write.csv(scores, gzfile(file.path(tables, paste0("external_scrna_", dataset, "_frozen_cell_scores.csv.gz"))), row.names = FALSE)
  donor <- aggregate(cbind(fssi_frozen_weighted, fssi_frozen_equal_weight,
                           pathological_component, repair_component_reversed,
                           Fib5_three_marker_baseline) ~ dataset_id + donor_id + condition,
                     data = scores, FUN = mean)
  cell_n <- aggregate(cell_id ~ dataset_id + donor_id + condition, data = scores, FUN = length)
  names(cell_n)[4] <- "fibroblast_cell_n"
  donor <- merge(donor, cell_n, by = c("dataset_id", "donor_id", "condition"))
  write.csv(donor, file.path(tables, paste0("external_scrna_", dataset, "_frozen_donor_scores.csv")), row.names = FALSE)
} else {
  empty <- data.frame(dataset_id = character(), sample_id = character(), donor_id = character(),
                      condition = character(), cell_id = character(), prediction_score = numeric(),
                      fssi_frozen_weighted = numeric(), fssi_frozen_equal_weight = numeric())
  write.csv(empty, gzfile(file.path(tables, paste0("external_scrna_", dataset, "_frozen_cell_scores.csv.gz"))), row.names = FALSE)
}
print(audits)
