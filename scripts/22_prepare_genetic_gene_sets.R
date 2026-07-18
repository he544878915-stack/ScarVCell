suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

root_env <- Sys.getenv("SCARVCELL_ROOT", unset = "")
if (nzchar(root_env)) {
  root <- gsub("\\\\", "/", root_env)
} else {
  script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  script_path <- sub("^--file=", "", script_arg[[1]])
  root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)
}
model <- read.csv(file.path(root, "results", "tables", "fssi_frozen_model.csv"), stringsAsFactors = FALSE)

sets <- list(
  FSSI_ALL = model$gene,
  FSSI_PATHOLOGICAL = model$gene[model$program == "pathological"],
  FSSI_REPAIR = model$gene[model$program == "repair"],
  FSSI_BOOTSTRAP_CORE = c("ADAM12", "COMP", "POSTN", "ASPN", "APOD", "PI16", "DCN", "ABCA8"),
  EXTENDED_CANDIDATES = c("RUNX2", "POSTN", "TGFB1", "CSF1", "CXCL12", "MIF", "APOD", "PI16", "PTGDS", "PDGFB")
)
symbols <- unique(unlist(sets))
entrez <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = symbols, column = "ENTREZID",
                                keytype = "SYMBOL", multiVals = "first")
mapping <- data.frame(symbol = symbols, entrez_id = unname(entrez[symbols]), stringsAsFactors = FALSE)
if (anyNA(mapping$entrez_id)) stop("Unmapped symbols: ", paste(mapping$symbol[is.na(mapping$entrez_id)], collapse = ", "))

out_dir <- file.path(root, ".stage", "genetics", "MAGMA", "analysis")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(mapping, file.path(out_dir, "locked_symbol_entrez_mapping.csv"), row.names = FALSE)
lines <- vapply(names(sets), function(name) {
  ids <- mapping$entrez_id[match(sets[[name]], mapping$symbol)]
  paste(c(name, ids), collapse = " ")
}, character(1))
writeLines(lines, file.path(out_dir, "locked_gene_sets_magma.txt"), useBytes = TRUE)

membership <- do.call(rbind, lapply(names(sets), function(name) {
  data.frame(set = name, symbol = sets[[name]],
             entrez_id = mapping$entrez_id[match(sets[[name]], mapping$symbol)], stringsAsFactors = FALSE)
}))
write.csv(membership, file.path(root, "results", "tables", "locked_genetic_gene_sets.csv"), row.names = FALSE)
cat("Prepared", length(sets), "locked gene sets and", nrow(mapping), "unique genes\n")
