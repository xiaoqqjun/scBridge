
#' Export Seurat Object to Bridge Intermediate Files
#'
#' @param object A Seurat object
#' @param bridge_dir Directory to store intermediate files. If NULL, uses tempdir.
#' @param slots Character vector of assay slots to export. Default c("counts", "data").
#' @param assay Name of the assay to export. Default is the default assay.
#' @param verbose Logical, print progress messages. Default TRUE.
#' @return The path to the bridge directory (invisibly).
#' @export
seurat_to_bridge <- function(object, bridge_dir = NULL,
                              slots = c("counts", "data"),
                              assay = NULL, verbose = TRUE) {

  if (!inherits(object, "Seurat")) {
    stop("Input must be a Seurat object.")
  }

  if (is.null(assay)) assay <- DefaultAssay(object)

  if (is.null(bridge_dir)) {
    bridge_dir <- file.path(tempdir(), paste0("scBridge_", format(Sys.time(), "%Y%m%d_%H%M%S")))
  }
  dir.create(bridge_dir, showWarnings = FALSE, recursive = TRUE)

  if (verbose) message("=== scBridge: Seurat -> Bridge ===")
  if (verbose) message("Output: ", bridge_dir)

  # 1. 表达矩阵 - 使用兼容函数
  for (s in slots) {
    mat <- tryCatch(
      .get_assay_data_compat(object, layer = s, assay = assay),
      error = function(e) NULL
    )
    if (is.null(mat) || prod(dim(mat)) == 0) {
      if (verbose) message("  Slot \'", s, "\': skipped (empty)")
      next
    }
    if (!inherits(mat, "dgCMatrix")) mat <- as(mat, "dgCMatrix")
    Matrix::writeMM(mat, file.path(bridge_dir, paste0("matrix_", s, ".mtx")))
    if (verbose) message("  Slot \'", s, "\': ", nrow(mat), " x ", ncol(mat))
  }

  # 2. 基因名和细胞名
  writeLines(rownames(object), file.path(bridge_dir, "genes.txt"))
  writeLines(colnames(object), file.path(bridge_dir, "barcodes.txt"))
  if (verbose) message("  Genes: ", nrow(object), " | Cells: ", ncol(object))

  # 3. Metadata（全部列）
  write.csv(object@meta.data, file.path(bridge_dir, "metadata.csv"))
  if (verbose) message("  Metadata: ", ncol(object@meta.data), " columns")

  # 4. 降维结果
  reductions <- names(object@reductions)
  if (length(reductions) > 0) {
    for (rd in reductions) {
      emb <- Embeddings(object, reduction = rd)
      rd_key <- Key(object[[rd]])
      write.csv(emb, file.path(bridge_dir, paste0("reduction_", rd, ".csv")))
      cat(rd, "\t", rd_key, "\t", ncol(emb), "\n",
          file = file.path(bridge_dir, "reduction_keys.txt"),
          append = (rd != reductions[1]))
    }
    if (verbose) message("  Reductions: ", paste(reductions, collapse = ", "))
  }

  # 5. Variable features
  var_features <- VariableFeatures(object, assay = assay)
  if (length(var_features) > 0) {
    writeLines(var_features, file.path(bridge_dir, "variable_features.txt"))
    if (verbose) message("  Variable features: ", length(var_features))
  }

  if (verbose) {
    message("=== Export complete ===")
    message("Files: ", paste(list.files(bridge_dir), collapse = ", "))
  }

  invisible(bridge_dir)
}

