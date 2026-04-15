
#' Convert h5ad File to Bridge Intermediate Files
#'
#' @param h5ad_file Path to the input h5ad file.
#' @param bridge_dir Directory to store intermediate files. If NULL, uses tempdir.
#' @param verbose Logical, print progress messages. Default TRUE.
#' @return The path to the bridge directory (invisibly).
#' @export
h5ad_to_bridge <- function(h5ad_file, bridge_dir = NULL, verbose = TRUE) {

  if (!file.exists(h5ad_file)) stop("File not found: ", h5ad_file)
  .check_python_deps()

  if (is.null(bridge_dir)) {
    bridge_dir <- file.path(tempdir(), paste0("scBridge_", format(Sys.time(), "%Y%m%d_%H%M%S")))
  }
  dir.create(bridge_dir, showWarnings = FALSE, recursive = TRUE)

  if (verbose) message("=== scBridge: h5ad -> Bridge ===")

  # 用 py_run_string 处理所有 Python 操作，避免 S4 类问题
  h5ad_norm <- normalizePath(h5ad_file, winslash = "/")
  bridge_norm <- normalizePath(bridge_dir, winslash = "/", mustWork = FALSE)

  # 读取 h5ad 并导出
  reticulate::py_run_string(sprintf("
import scanpy as sc
import scipy.io
import scipy.sparse
import pandas as pd
import numpy as np
import os

_scb_adata = sc.read_h5ad(\'%s\')
_scb_bridge = \'%s\'

# 1. 表达矩阵 (genes x cells)
_scb_X = _scb_adata.X
if not scipy.sparse.issparse(_scb_X):
    _scb_X = scipy.sparse.csc_matrix(_scb_X)
scipy.io.mmwrite(os.path.join(_scb_bridge, \'matrix_X.mtx\'), _scb_X.T)

# 2. Layers
for lk in list(_scb_adata.layers.keys()):
    _scb_layer = _scb_adata.layers[lk]
    if not scipy.sparse.issparse(_scb_layer):
        _scb_layer = scipy.sparse.csc_matrix(_scb_layer)
    scipy.io.mmwrite(os.path.join(_scb_bridge, f\'layer_{lk}.mtx\'), _scb_layer.T)

# 3. 基因名和细胞名
with open(os.path.join(_scb_bridge, \'genes.txt\'), \'w\') as f:
    f.write(\'\\n\'.join(_scb_adata.var_names.tolist()) + \'\\n\')
with open(os.path.join(_scb_bridge, \'barcodes.txt\'), \'w\') as f:
    f.write(\'\\n\'.join(_scb_adata.obs_names.tolist()) + \'\\n\')

# 4. Metadata
_scb_adata.obs.to_csv(os.path.join(_scb_bridge, \'metadata.csv\'))

# 5. 降维
for ok in list(_scb_adata.obsm.keys()):
    _scb_emb = pd.DataFrame(_scb_adata.obsm[ok], index=_scb_adata.obs_names)
    _scb_emb.to_csv(os.path.join(_scb_bridge, f\'reduction_{ok}.csv\'))

# 6. Variable features
if \'highly_variable\' in _scb_adata.var.columns:
    _scb_hvg = _scb_adata.var_names[_scb_adata.var[\'highly_variable\']].tolist()
    if len(_scb_hvg) > 0:
        with open(os.path.join(_scb_bridge, \'variable_features.txt\'), \'w\') as f:
            f.write(\'\\n\'.join(_scb_hvg) + \'\\n\')

# 收集信息供 R 端使用
_scb_info = {
    \'n_obs\': _scb_adata.n_obs,
    \'n_vars\': _scb_adata.n_vars,
    \'obsm_keys\': list(_scb_adata.obsm.keys()),
    \'layer_keys\': list(_scb_adata.layers.keys()),
    \'meta_cols\': list(_scb_adata.obs.columns)
}

del _scb_adata, _scb_X
", h5ad_norm, bridge_norm))

  info <- reticulate::py$`_scb_info`
  reticulate::py_run_string("del _scb_info")

  if (verbose) {
    message("  Read: ", h5ad_file)
    message("  Shape: ", info$n_obs, " x ", info$n_vars)
    message("  Metadata: ", length(info$meta_cols), " columns")
    if (length(info$obsm_keys) > 0)
      message("  Reductions: ", paste(info$obsm_keys, collapse = ", "))
    if (length(info$layer_keys) > 0)
      message("  Layers: ", paste(info$layer_keys, collapse = ", "))
    message("=== Export complete ===")
    message("Files: ", paste(list.files(bridge_dir), collapse = ", "))
  }

  invisible(bridge_dir)
}

