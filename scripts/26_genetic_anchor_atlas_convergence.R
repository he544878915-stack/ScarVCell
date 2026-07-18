#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
})

root <- Sys.getenv("SCARVCELL_ROOT", unset = normalizePath(file.path(dirname(commandArgs()[1]), "..")))
tables <- file.path(root, "results", "tables")
object_file <- file.path(root, "data", "processed", "integrated", "scar_fibroblast_reference.rds")
anchor_file <- file.path(tables, "Greene2025_locked_susceptibility_anchor.csv")
model_file <- file.path(tables, "fssi_frozen_model.csv")
set.seed(20260713)

anchor <- read.csv(anchor_file, check.names = FALSE)
anchor$included_in_signed_score <- tolower(as.character(anchor$included_in_signed_score)) == "true"
anchor <- anchor[anchor$included_in_signed_score, , drop = FALSE]
model <- read.csv(model_file, check.names = FALSE)
object <- readRDS(object_file)
object <- JoinLayers(object, assay = "RNA")
metadata <- object@meta.data
metadata$unit_id <- paste(metadata$sample_id, metadata$seurat_clusters, sep = "__C")
units <- unique(metadata$unit_id)

data <- LayerData(object[["RNA"]], layer = "data")
design <- sparse.model.matrix(~ 0 + factor(metadata$unit_id, levels = units))
colnames(design) <- units
average <- data %*% design
average <- sweep(average, 2, Matrix::colSums(design), "/")
average <- as.matrix(average)

unit_meta <- unique(metadata[, c("unit_id", "dataset_id", "sample_id", "seurat_clusters")])
unit_meta <- unit_meta[match(colnames(average), unit_meta$unit_id), ]
unit_meta$cluster <- as.character(unit_meta$seurat_clusters)

within_sample_z <- function(matrix, meta) {
  output <- matrix
  for (sample in unique(meta$sample_id)) {
    columns <- which(meta$sample_id == sample)
    block <- matrix[, columns, drop = FALSE]
    gene_mean <- rowMeans(block)
    gene_sd <- apply(block, 1, sd)
    gene_sd[!is.finite(gene_sd) | gene_sd == 0] <- 1
    output[, columns] <- sweep(sweep(block, 1, gene_mean, "-"), 1, gene_sd, "/")
  }
  output
}

z <- within_sample_z(average, unit_meta)
available_anchor <- intersect(anchor$gene, rownames(z))
anchor_use <- anchor[match(available_anchor, anchor$gene), ]
anchor_use$signed_weight <- anchor_use$signed_weight / sum(abs(anchor_use$signed_weight))
genetic_score <- as.numeric(crossprod(anchor_use$signed_weight, z[available_anchor, , drop = FALSE]))

available_fssi <- intersect(model$gene, rownames(average))
model_use <- model[match(available_fssi, model$gene), ]
fssi_standardised <- sweep(average[available_fssi, , drop = FALSE], 1, model_use$training_mean, "-")
fssi_standardised <- sweep(fssi_standardised, 1, model_use$training_sd, "/")
fssi_weight <- model_use$weight / sum(abs(model_use$weight))
fssi_score <- as.numeric(crossprod(fssi_weight, fssi_standardised))

scores <- cbind(unit_meta[, c("dataset_id", "sample_id", "cluster", "unit_id")],
                genetic_susceptibility_expression_score = genetic_score,
                frozen_fssi_pseudobulk_score = fssi_score,
                anchor_genes_present = length(available_anchor),
                anchor_genes_total = nrow(anchor))
write.csv(scores, file.path(tables, "atlas_donor_cluster_genetic_anchor_scores.csv"), row.names = FALSE)

pathological <- c("1", "5", "11")
repair <- "2"
donor_rows <- list()
for (sample in unique(scores$sample_id)) {
  current <- scores[scores$sample_id == sample, ]
  path <- current[current$cluster %in% pathological, ]
  ref <- current[current$cluster == repair, ]
  if (nrow(path) && nrow(ref)) {
    donor_rows[[length(donor_rows) + 1]] <- data.frame(
      dataset_id = current$dataset_id[1], sample_id = sample,
      pathological_clusters_present = paste(sort(unique(path$cluster)), collapse = ";"),
      genetic_anchor_difference = mean(path$genetic_susceptibility_expression_score) -
        mean(ref$genetic_susceptibility_expression_score),
      fssi_difference = mean(path$frozen_fssi_pseudobulk_score) - mean(ref$frozen_fssi_pseudobulk_score)
    )
  }
}
donors <- do.call(rbind, donor_rows)
write.csv(donors, file.path(tables, "atlas_donor_genetic_anchor_contrasts.csv"), row.names = FALSE)

t_interval <- function(values) {
  estimate <- mean(values)
  error <- qt(0.975, df = length(values) - 1) * sd(values) / sqrt(length(values))
  c(estimate = estimate, low = estimate - error, high = estimate + error)
}
exact_signflip <- function(values) {
  signs <- as.matrix(expand.grid(rep(list(c(-1, 1)), length(values))))
  mean(abs(as.numeric(signs %*% values) / length(values)) >= abs(mean(values)) - 1e-12)
}

interval <- t_interval(donors$genetic_anchor_difference)
rho <- suppressWarnings(cor.test(donors$genetic_anchor_difference, donors$fssi_difference,
                                 method = "spearman", exact = FALSE))

# The null preserves each anchor's signed weight and samples genes from the same
# atlas-abundance ventile after excluding FSSI and anchor genes.
gene_mean <- rowMeans(average)
finite <- is.finite(gene_mean)
breaks <- unique(quantile(gene_mean[finite], probs = seq(0, 1, 0.05), na.rm = TRUE))
expression_bin <- cut(gene_mean, breaks = breaks, include.lowest = TRUE, labels = FALSE)
names(expression_bin) <- rownames(average)
exclude <- union(anchor$gene, model$gene)
pool <- setdiff(rownames(average)[is.finite(expression_bin)], exclude)

donor_delta_by_gene <- matrix(NA_real_, nrow = nrow(z), ncol = nrow(donors),
                              dimnames = list(rownames(z), donors$sample_id))
for (index in seq_len(nrow(donors))) {
  current <- unit_meta$sample_id == donors$sample_id[index]
  path_columns <- which(current & unit_meta$cluster %in% pathological)
  repair_columns <- which(current & unit_meta$cluster == repair)
  donor_delta_by_gene[, index] <- rowMeans(z[, path_columns, drop = FALSE]) -
    rowMeans(z[, repair_columns, drop = FALSE])
}

iterations <- 10000
null_mean <- numeric(iterations)
for (iteration in seq_len(iterations)) {
  sampled <- character(nrow(anchor_use))
  for (index in seq_len(nrow(anchor_use))) {
    bin <- expression_bin[anchor_use$gene[index]]
    candidates <- setdiff(pool[expression_bin[pool] == bin], sampled)
    if (!length(candidates)) stop("No expression-matched candidate gene available")
    sampled[index] <- sample(candidates, 1)
  }
  null_donor <- as.numeric(crossprod(anchor_use$signed_weight,
                                    donor_delta_by_gene[sampled, , drop = FALSE]))
  null_mean[iteration] <- mean(null_donor)
}
observed <- mean(donors$genetic_anchor_difference)
null_summary <- data.frame(
  observed_mean_difference = observed,
  matched_null_mean = mean(null_mean),
  matched_null_sd = sd(null_mean),
  matched_null_q025 = quantile(null_mean, 0.025),
  matched_null_q975 = quantile(null_mean, 0.975),
  empirical_p_greater = (1 + sum(null_mean >= observed)) / (iterations + 1),
  empirical_p_two_sided = (1 + sum(abs(null_mean) >= abs(observed))) / (iterations + 1),
  iterations = iterations
)
write.csv(data.frame(iteration = seq_len(iterations), null_mean_difference = null_mean),
          gzfile(file.path(tables, "atlas_genetic_anchor_matched_null.csv.gz")), row.names = FALSE)
write.csv(null_summary, file.path(tables, "atlas_genetic_anchor_matched_null_summary.csv"), row.names = FALSE)

by_dataset <- aggregate(genetic_anchor_difference ~ dataset_id, donors,
                        function(x) c(n = length(x), mean = mean(x), positive = sum(x > 0)))
by_dataset <- do.call(data.frame, by_dataset)
names(by_dataset) <- c("dataset_id", "donors", "mean_difference", "positive_donors")
write.csv(by_dataset, file.path(tables, "atlas_genetic_anchor_by_dataset.csv"), row.names = FALSE)

summary <- data.frame(
  donors = nrow(donors), datasets = length(unique(donors$dataset_id)),
  anchor_genes_present = length(available_anchor), anchor_genes_total = nrow(anchor),
  mean_pathological_minus_repair = interval["estimate"],
  ci95_low = interval["low"], ci95_high = interval["high"],
  exact_signflip_p_two_sided = exact_signflip(donors$genetic_anchor_difference),
  positive_donors = sum(donors$genetic_anchor_difference > 0),
  spearman_with_fssi_difference = unname(rho$estimate),
  spearman_p = rho$p.value,
  matched_null_p_greater = null_summary$empirical_p_greater,
  matched_null_p_two_sided = null_summary$empirical_p_two_sided
)
write.csv(summary, file.path(tables, "atlas_genetic_anchor_convergence_summary.csv"), row.names = FALSE)
print(summary)
print(by_dataset)
