
#' scBridge: Seamless Conversion Between Seurat and AnnData
#'
#' @description
#' Provides reliable bidirectional conversion between Seurat (R) and
#' AnnData (Python) objects. Bypasses known bugs in SeuratDisk and MuDataSeurat
#' by using a safe MTX + CSV intermediate format with automatic temporary
#' directory management.
#'
#' @docType package
#' @name scBridge-package
"_PACKAGE"


#' Check Python dependencies
#' @noRd
.check_python_deps <- function() {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("Package reticulate is required. Install with: install.packages(\"reticulate\")")
  }

  required_modules <- c("scanpy", "pandas", "numpy", "scipy")
  missing <- character(0)

  for (mod in required_modules) {
    if (!reticulate::py_module_available(mod)) {
      missing <- c(missing, mod)
    }
  }

  if (length(missing) > 0) {
    stop("Missing Python modules: ", paste(missing, collapse = ", "),
         "\nInstall with: pip install ", paste(missing, collapse = " "),
         call. = FALSE)
  }
}


#' Get assay data with Seurat v4/v5 compatibility
#' @param object A Seurat object
#' @param layer Layer/slot name (e.g., "counts", "data")
#' @param assay Assay name
#' @return The assay data matrix
#' @noRd
.get_assay_data_compat <- function(object, layer, assay = NULL) {
  if (is.null(assay)) assay <- DefaultAssay(object)

  # 检测 Seurat 版本
  seurat_v5 <- packageVersion("Seurat") >= "5.0.0"

  if (seurat_v5) {
    # V5 API: use layer parameter
    return(GetAssayData(object, assay = assay, layer = layer))
  } else {
    # V4 API: use slot parameter
    return(GetAssayData(object, assay = assay, slot = layer))
  }
}


#' Set assay data with Seurat v4/v5 compatibility
#' @param object A Seurat object
#' @param layer Layer/slot name
#' @param value New data matrix
#' @param assay Assay name
#' @return Updated Seurat object
#' @noRd
.set_assay_data_compat <- function(object, layer, value, assay = NULL) {
  if (is.null(assay)) assay <- DefaultAssay(object)

  # 检测 Seurat 版本
  seurat_v5 <- packageVersion("Seurat") >= "5.0.0"

  if (seurat_v5) {
    # V5 API: use layer parameter
    return(SetAssayData(object, assay = assay, layer = layer, new.data = value))
  } else {
    # V4 API: use slot parameter
    return(SetAssayData(object, assay = assay, slot = layer, new.data = value))
  }
}


#' Ensure matrix is dgCMatrix format
#' readMM may return dgTMatrix which is not compatible with Seurat slots
#' @param mat A matrix object
#' @return A dgCMatrix
#' @noRd
.ensure_dgCMatrix <- function(mat) {
  if (inherits(mat, "dgCMatrix")) {
    return(mat)
  } else if (inherits(mat, "dgTMatrix")) {
    return(as(mat, "CsparseMatrix"))
  } else if (inherits(mat, "Matrix")) {
    return(as(mat, "CsparseMatrix"))
  } else {
    return(Matrix::Matrix(mat, sparse = TRUE))
  }
}

