suppressPackageStartupMessages({
  library(Seurat)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

root <- Sys.getenv("SCARVCELL_ROOT", unset = normalizePath(getwd(), winslash = "/", mustWork = TRUE))
tables <- file.path(root, "results", "tables")
stage <- file.path(Sys.getenv("SCARVCELL_STAGE", unset = file.path(root, ".stage")), "external", "GSE191067")
dir.create(stage, recursive = TRUE, showWarnings = FALSE)

reference <- readRDS(file.path(root, "data", "processed", "integrated",
                               "scar_scrna_compact_label_transfer_reference.rds"))
hvg <- VariableFeatures(reference)
if (!length(hvg)) stop("The compact reference has no stored variable features")
model <- read.csv(file.path(tables, "fssi_frozen_model.csv"), stringsAsFactors = FALSE)

marker_sets <- list(
  fibroblast = c("COL1A1", "COL1A2", "COL3A1", "COL6A1", "DCN", "LUM", "PDGFRA", "PDGFRB"),
  immune = c("PTPRC", "CD3D", "CD3E", "CD79A", "MS4A1", "NKG7", "LST1", "TYROBP"),
  endothelial = c("PECAM1", "VWF", "KDR", "EMCN", "RAMP2", "PLVAP"),
  epithelial = c("EPCAM", "KRT5", "KRT14", "KRT1", "KRT10", "KRT19"),
  neural = c("S100B", "SOX10", "PLP1", "MPZ", "SCD", "PMP22"),
  myeloid = c("LST1", "TYROBP", "AIF1", "FCER1G", "CTSS", "CD74", "CSF1R"),
  communication = c("POSTN", "ITGAV", "ITGB1", "TGFB1", "TGFBR1", "TGFBR2",
                    "CXCL12", "CXCR4", "CSF1", "CSF1R", "PDGFB", "PDGFRA", "MIF", "CD74"),
  collagen_ecm = c("COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL12A1", "FN1", "POSTN")
)
prior_path <- file.path(root, "data", "external", "nichenet", "ligand_target_matrix_nsga2r_final.rds")
ligand_targets <- character()
if (file.exists(prior_path)) {
  prior <- readRDS(prior_path)
  for (ligand in c("POSTN", "TGFB1", "CXCL12", "CSF1", "PDGFB", "MIF")) {
    if (ligand %in% colnames(prior)) {
      values <- prior[, ligand]
      ligand_targets <- c(ligand_targets, names(sort(values, decreasing = TRUE))[seq_len(min(50L, length(values)))])
    }
  }
  marker_sets$nichenet_top50_targets <- unique(ligand_targets)
}
desired <- unique(toupper(c(hvg, model$gene, unlist(marker_sets, use.names = FALSE))))

mapping <- AnnotationDbi::select(org.Hs.eg.db, keys = desired, keytype = "SYMBOL",
                                columns = c("ENSEMBL", "SYMBOL"))
mapping <- mapping[!is.na(mapping$ENSEMBL) & !is.na(mapping$SYMBOL), c("ENSEMBL", "SYMBOL")]
mapping$ENSEMBL <- sub("\\..*$", "", mapping$ENSEMBL)
mapping$SYMBOL <- toupper(mapping$SYMBOL)
mapping <- unique(mapping)
mapping$priority <- match(mapping$SYMBOL, desired)
mapping <- mapping[order(mapping$priority, mapping$ENSEMBL), ]
mapping$priority <- NULL
write.csv(mapping, file.path(stage, "selected_ensembl_symbol_map.csv"), row.names = FALSE)

all_map <- AnnotationDbi::select(org.Hs.eg.db, keys = keys(org.Hs.eg.db, keytype = "ENSEMBL"),
                                keytype = "ENSEMBL", columns = "SYMBOL")
all_map <- all_map[!is.na(all_map$SYMBOL), c("ENSEMBL", "SYMBOL")]
all_map$ENSEMBL <- sub("\\..*$", "", all_map$ENSEMBL)
all_map$SYMBOL <- toupper(all_map$SYMBOL)
all_map <- unique(all_map)
write.csv(all_map, file.path(stage, "all_ensembl_symbol_map.csv"), row.names = FALSE)

membership <- do.call(rbind, lapply(names(marker_sets), function(set_name) {
  data.frame(set = set_name, gene = marker_sets[[set_name]], stringsAsFactors = FALSE)
}))
write.csv(membership, file.path(stage, "locked_marker_sets.csv"), row.names = FALSE)
writeLines(desired, file.path(stage, "desired_symbols.txt"))
cat("Reference HVGs:", length(hvg), "\nDesired symbols:", length(desired),
    "\nMapped symbols:", length(unique(mapping$SYMBOL)), "\n")
