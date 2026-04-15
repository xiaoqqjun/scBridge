
#' Convert Bridge Intermediate Files to h5ad
#'
#' @param bridge_dir Path to the bridge intermediate directory.
#' @param output_file Path for the output h5ad file.
#' @param verbose Logical, print progress messages. Default TRUE.
#' @return The output file path (invisibly).
#' @export
bridge_to_h5ad <- function(bridge_dir, output_file, verbose = TRUE) {

  if (!dir.exists(bridge_dir)) stop("Bridge directory not found: ", bridge_dir)
  if (verbose) message("=== scBridge: Bridge -> h5ad ===")

  .check_python_deps()

  bridge_norm <- normalizePath(bridge_dir, winslash = "/")
  out_norm <- normalizePath(dirname(output_file), winslash = "/", mustWork = FALSE)
  out_path <- paste0(out_norm, "/", basename(output_file))

  # 全部在 Python 端完成，避免 reticulate S4 类问题
  reticulate::py_run_string(sprintf("
import scanpy as sc
import scipy.io
import scipy.sparse
import pandas as pd
import numpy as np
import os
import glob

_scb_bridge = \'%s\'
_scb_output = \'%s\'

# 1. 找主矩阵 (counts 优先)
_scb_mtx_files = glob.glob(os.path.join(_scb_bridge, \'matrix_*.mtx\'))
_scb_mtx_names = {os.path.basename(f): f for f in _scb_mtx_files}

if \'matrix_counts.mtx\' in _scb_mtx_names:
    _scb_primary = \'counts\'
else:
    _scb_primary = os.path.basename(_scb_mtx_files[0]).replace(\'matrix_\', \'\').replace(\'.mtx\', \'\')

_scb_X = scipy.sparse.csc_matrix(scipy.io.mmread(os.path.join(_scb_bridge, f\'matrix_{_scb_primary}.mtx\')).T)

# 2. 基因名和细胞名
with open(os.path.join(_scb_bridge, \'genes.txt\')) as f:
    _scb_genes = [line.strip() for line in f if line.strip()]
with open(os.path.join(_scb_bridge, \'barcodes.txt\')) as f:
    _scb_barcodes = [line.strip() for line in f if line.strip()]

# 3. Metadata
_scb_meta = pd.read_csv(os.path.join(_scb_bridge, \'metadata.csv\'), index_col=0)

# 4. 组装 AnnData
_scb_adata = sc.AnnData(X=_scb_X)
_scb_adata.obs_names = pd.Index(_scb_barcodes)
_scb_adata.var_names = pd.Index(_scb_genes)
_scb_adata.obs = _scb_meta

# 5. 其他 slots 作为 layers
for _scb_fname, _scb_fpath in _scb_mtx_names.items():
    _scb_slot = _scb_fname.replace(\'matrix_\', \'\').replace(\'.mtx\', \'\')
    if _scb_slot == _scb_primary:
        continue
    _scb_adata.layers[_scb_slot] = scipy.sparse.csc_matrix(scipy.io.mmread(_scb_fpath).T)

# 6. 降维结果
_scb_rd_files = glob.glob(os.path.join(_scb_bridge, \'reduction_*.csv\'))
for _scb_rf in _scb_rd_files:
    _scb_rd_name = os.path.basename(_scb_rf).replace(\'reduction_\', \'\').replace(\'.csv\', \'\')
    _scb_emb = pd.read_csv(_scb_rf, index_col=0).values
    _scb_adata.obsm[f\'X_{_scb_rd_name}\'] = _scb_emb

# 7. Variable features
_scb_vf_path = os.path.join(_scb_bridge, \'variable_features.txt\')
if os.path.exists(_scb_vf_path):
    with open(_scb_vf_path) as f:
        _scb_hvg = [line.strip() for line in f if line.strip()]
    _scb_adata.var[\'highly_variable\'] = _scb_adata.var_names.isin(_scb_hvg)

# 8. 保存
_scb_adata.write_h5ad(_scb_output)

# 收集信息
_scb_info = {
    \'n_obs\': int(_scb_adata.n_obs),
    \'n_vars\': int(_scb_adata.n_vars),
    \'primary_slot\': _scb_primary,
    \'layers\': list(_scb_adata.layers.keys()),
    \'obsm\': list(_scb_adata.obsm.keys()),
    \'hvg\': int(_scb_adata.var[\'highly_variable\'].sum()) if \'highly_variable\' in _scb_adata.var.columns else 0
}

del _scb_adata, _scb_X
", bridge_norm, out_path))

  info <- reticulate::py$`_scb_info`
  reticulate::py_run_string("del _scb_info")

  if (verbose) {
    message("  X (", info$primary_slot, "): ", info$n_obs, " x ", info$n_vars)
    if (length(info$layers) > 0)
      message("  Layers: ", paste(info$layers, collapse = ", "))
    if (length(info$obsm) > 0)
      message("  obsm: ", paste(info$obsm, collapse = ", "))
    if (info$hvg > 0)
      message("  HVG: ", info$hvg)
    message("=== Saved: ", out_path, " ===")
  }

  invisible(output_file)
}

