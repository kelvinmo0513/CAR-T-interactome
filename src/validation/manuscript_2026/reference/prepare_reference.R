#!/usr/bin/env Rscript
# Recreate the selected SCT reference from an external annotated Seurat object.
# The external cell map fixes the historical selection and order; no resampling.
suppressPackageStartupMessages({library(Seurat); library(jsonlite)})
args <- commandArgs(trailingOnly=TRUE)
usage <- paste("Usage: Rscript prepare_reference.R --rds FILE",
               "--selected-cells CSV --output NEW_DIRECTORY")
if (length(args) != 6 || anyDuplicated(args[seq(1, 6, 2)])) stop(usage)
opt <- setNames(as.list(args[seq(2, 6, 2)]), args[seq(1, 6, 2)])
if (!setequal(names(opt), c("--rds", "--selected-cells", "--output"))) stop(usage)
out <- opt[["--output"]]
script_arg <- grep("^--file=", commandArgs(trailingOnly=FALSE), value=TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[1]), mustWork=TRUE)
repository <- dirname(dirname(dirname(dirname(dirname(script_path)))))
out_absolute <- path.expand(out)
if (!startsWith(out_absolute, "/")) out_absolute <- file.path(getwd(), out_absolute)
out_absolute <- normalizePath(out_absolute, mustWork=FALSE)
if (identical(out_absolute, repository) || startsWith(out_absolute, paste0(repository, "/"))) {
  stop("Reference outputs include cell identifiers; choose a directory outside the repository")
}
if (file.exists(out) || dir.exists(out)) stop("Output must be a new directory")
obj <- readRDS(opt[["--rds"]])
map <- read.csv(opt[["--selected-cells"]], stringsAsFactors=FALSE, check.names=FALSE)
required <- c("source_cell", "phenotype", "source_group")
if (!all(required %in% names(map))) stop("Selection map lacks required columns")
stopifnot(nrow(map) == 3759, !anyNA(map[required]),
          !anyDuplicated(map$source_cell))
rna <- LayerData(obj, assay="RNA", layer="counts")
sct <- LayerData(obj, assay="SCT", layer="data")
stopifnot(identical(colnames(rna), rownames(obj@meta.data)),
          identical(colnames(rna), colnames(sct)),
          all(c("predicted_celltype_MI", "orig.ident") %in% names(obj@meta.data)),
          all(c("CCL8", "CCL13") %in% rownames(rna)))
labels <- as.character(obj$predicted_celltype_MI)
labels[labels == "NK cells"] <- "T cells"
labels[labels %in% c("Monocytes/cDC", "pDC")] <- "Myeloid DC"
gate <- as.vector(rna["CCL8", ] > 0 & rna["CCL13", ] > 0)
stopifnot(sum(gate) == 2, all(labels[gate] == "Myeloid DC"),
          length(unique(as.character(obj$orig.ident)[gate])) == 1)
labels[gate] <- "CCL8+CCL13+ Myeloid DC"
selected <- match(map$source_cell, colnames(sct))
stopifnot(!anyNA(selected), identical(labels[selected], map$phenotype),
          identical(as.character(obj$orig.ident)[selected], map$source_group),
          all(which(gate) %in% selected), nrow(sct) == 20674,
          !("CCL13" %in% rownames(sct)), !anyDuplicated(rownames(sct)))
x <- as.matrix(sct[sort(rownames(sct)), selected, drop=FALSE])
# Historical code clamped negatives, but none exist in this selected SCT layer.
# Fail rather than silently transform a different input object.
stopifnot(all(is.finite(x)), all(x >= 0), ncol(x) == 3759)
colnames(x) <- map$phenotype
names(dimnames(x)) <- c("GeneSymbol", "")
dir.create(out, recursive=TRUE)
write.table(x, file.path(out, "Kelvin_DP_SCT_reference.txt"), sep="\t",
            quote=FALSE, row.names=TRUE, col.names=NA)
write.csv(map[required], file.path(out, "selected_cells_LOCAL_ONLY.csv"), row.names=FALSE)
write.csv(as.data.frame(table(map$phenotype)), file.path(out, "reference_counts.csv"), row.names=FALSE)
audit <- list(analysis_label="oversampled two-cell AND reference—exploratory",
  passed=TRUE, RNA_features=nrow(rna), SCT_features=nrow(sct),
  selected_cells=ncol(x), background_cells=3757L, unique_target_cells=2L,
  target_source_groups=1L, CCL13_in_SCT=FALSE, expression_layer="SCT/data",
  gate="RNA counts CCL8 > 0 AND CCL13 > 0 in merged Myeloid DC",
  selection="external fixed cell map; order preserved; no new sampling",
  historical_selection_seed=20260930L, original_Kelvin_sampling_seed="unknown",
  negative_values_clamped=0L, finite_nonnegative=TRUE, clinical_outcomes_used=FALSE,
  dependencies=list(R=R.version.string, Seurat=as.character(packageVersion("Seurat")),
                    SeuratObject=as.character(packageVersion("SeuratObject"))))
write_json(audit, file.path(out,"reference_audit.json"), pretty=TRUE, auto_unbox=TRUE)
capture.output(sessionInfo(), file=file.path(out, "R_sessionInfo.txt"))
cat("PASS: 20,674 features; 3,757 background cells and 2 target cells.\n")
