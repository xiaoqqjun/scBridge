
#' Convert Seurat Object to h5ad File
#'
#' One-step conversion from Seurat to h5ad format.
#'
#' @param object A Seurat object.
#' @param file Output h5ad file path.
#' @param slots Character vector of assay slots to export. Default c("counts", "data").
#' @param assay Name of the assay to export. Default is the default assay.
#' @param verbose Logical, print progress messages. Default TRUE.
#' @return The output file path (invisibly).
#' @export
#' @examples
#' \dontrun{
#' seurat2py(test_obj, "output.h5ad")
#' }
seurat2py <- function(object, file, slots = c("counts", "data"),
                     assay = NULL, verbose = TRUE) {
  bridge_dir <- seurat_to_bridge(object, bridge_dir = NULL, slots = slots,
                                  assay = assay, verbose = verbose)
  bridge_to_h5ad(bridge_dir, output_file = file, verbose = verbose)
  unlink(bridge_dir, recursive = TRUE)
  if (verbose) message("Temp files cleaned.")
  invisible(file)
}


#' Convert h5ad File to Seurat Object
#'
#' One-step conversion from h5ad to Seurat format.
#'
#' @param file Input h5ad file path.
#' @param verbose Logical, print progress messages. Default TRUE.
#' @return A Seurat object.
#' @export
#' @examples
#' \dontrun{
#' sobj <- py2seurat("input.h5ad")
#' DimPlot(sobj, reduction = "umap")
#' }
py2seurat <- function(file, verbose = TRUE) {
  bridge_dir <- h5ad_to_bridge(file, bridge_dir = NULL, verbose = verbose)
  sobj <- bridge_to_seurat(bridge_dir, verbose = verbose)
  unlink(bridge_dir, recursive = TRUE)
  if (verbose) message("Temp files cleaned.")
  return(sobj)
}

