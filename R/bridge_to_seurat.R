
#' Convert Bridge Intermediate Files to Seurat Object
#'
#' @param bridge_dir Path to the bridge intermediate directory.
#' @param verbose Logical, print progress messages. Default TRUE.
#' @return A Seurat object.
#' @export
bridge_to_seurat <- function(bridge_dir, verbose = TRUE) {

  if (!dir.exists(bridge_dir)) stop("Bridge directory not found: ", bridge_dir)
  if (verbose) message("=== scBridge: Bridge -> Seurat ===")

  # 1. 基因名和细胞名
  genes <- readLines(file.path(bridge_dir, "genes.txt"))
  barcodes <- readLines(file.path(bridge_dir, "barcodes.txt"))
  # 去掉可能的空行
  genes <- genes[nchar(genes) > 0]
  barcodes <- barcodes[nchar(barcodes) > 0]
  if (verbose) message("  Genes: ", length(genes), " | Cells: ", length(barcodes))

  # 2. 表达矩阵
  mtx_files <- list.files(bridge_dir, pattern = "^(matrix_|layer_).*\\.mtx$")

  primary_file <- NULL
  for (candidate in c("matrix_counts.mtx", "matrix_X.mtx", "layer_counts.mtx")) {
    if (candidate %in% mtx_files) {
      primary_file <- candidate
      break
    }
  }
  if (is.null(primary_file)) primary_file <- mtx_files[1]
  primary_slot <- gsub("matrix_|layer_|\\.mtx", "", primary_file)

  X <- Matrix::readMM(file.path(bridge_dir, primary_file))
  X <- .ensure_dgCMatrix(X)
  rownames(X) <- genes
  colnames(X) <- barcodes
  if (verbose) message("  Matrix (", primary_slot, "): ", nrow(X), " x ", ncol(X))

  # 3. 创建 Seurat 对象
  sobj <- CreateSeuratObject(counts = X)

  # 4. 其他矩阵 - 使用兼容函数
  for (mf in mtx_files) {
    if (mf == primary_file) next
    slot_name <- gsub("matrix_|layer_|\\.mtx", "", mf)
    extra_mat <- Matrix::readMM(file.path(bridge_dir, mf))
    extra_mat <- .ensure_dgCMatrix(extra_mat)
    rownames(extra_mat) <- genes
    colnames(extra_mat) <- barcodes
    if (slot_name == "data") {
      sobj <- .set_assay_data_compat(sobj, layer = "data", value = extra_mat)
      if (verbose) message("  Slot \'data\': loaded")
    } else if (slot_name == "scale.data") {
      sobj <- .set_assay_data_compat(sobj, layer = "scale.data", value = as.matrix(extra_mat))
      if (verbose) message("  Slot \'scale.data\': loaded")
    }
  }

  # 5. Metadata
  meta_file <- file.path(bridge_dir, "metadata.csv")
  if (file.exists(meta_file)) {
    meta <- read.csv(meta_file, row.names = 1)
    auto_cols <- c("orig.ident", "nCount_RNA", "nFeature_RNA")
    meta <- meta[, !colnames(meta) %in% auto_cols, drop = FALSE]
    if (ncol(meta) > 0) sobj <- AddMetaData(sobj, metadata = meta)
    if (verbose) message("  Metadata: ", ncol(meta), " columns added")
  }

  # 6. 降维结果
  rd_files <- list.files(bridge_dir, pattern = "^reduction_.*\\.csv$")
  if (length(rd_files) > 0) {
    # 读取 key 映射
    key_map <- list()
    key_file <- file.path(bridge_dir, "reduction_keys.txt")
    if (file.exists(key_file)) {
      key_lines <- readLines(key_file)
      for (line in key_lines) {
        parts <- strsplit(trimws(line), "\\t")[[1]]
        if (length(parts) >= 2) key_map[[parts[1]]] <- parts[2]
      }
    }

    for (f in rd_files) {
      rd_name_raw <- gsub("reduction_|\\.csv", "", f)
      emb <- as.matrix(read.csv(file.path(bridge_dir, f), row.names = 1))
      rd_name <- gsub("^X_", "", rd_name_raw)

      if (rd_name %in% names(key_map)) {
        rd_key <- key_map[[rd_name]]
      } else {
        rd_key <- .infer_reduction_key(rd_name)
      }

      colnames(emb) <- paste0(rd_key, seq_len(ncol(emb)))

      sobj[[rd_name]] <- CreateDimReducObject(
        embeddings = emb, key = rd_key, assay = DefaultAssay(sobj)
      )
      if (verbose) message("  Reduction \'", rd_name, "\' (key=", rd_key, "): ",
                           nrow(emb), " x ", ncol(emb))
    }
  }

  # 7. Variable features
  vf_file <- file.path(bridge_dir, "variable_features.txt")
  if (file.exists(vf_file)) {
    var_features <- readLines(vf_file)
    var_features <- var_features[nchar(var_features) > 0]
    VariableFeatures(sobj) <- var_features
    if (verbose) message("  Variable features: ", length(var_features))
  }

  # 8. 自动设置 Idents
  # 按优先级选择：seurat_clusters > orig.ident > 第一个非自动生成的列
  idents_cols <- c("seurat_clusters", "orig.ident")
  meta_cols <- colnames(sobj@meta.data)

  # 找到合适的列
  idents_col <- NULL
  for (col in idents_cols) {
    if (col %in% meta_cols) {
      idents_col <- col
      break
    }
  }

  # 如果没有找到，选择第一个非自动生成的列
  if (is.null(idents_col)) {
    auto_cols <- c("orig.ident", "nCount_RNA", "nFeature_RNA", "nFeature_Spatial",
                    "percent.mt", "percent.ribo", "percent hb")
    non_auto_cols <- meta_cols[!meta_cols %in% auto_cols]
    if (length(non_auto_cols) > 0) {
      idents_col <- non_auto_cols[1]
    }
  }

  # 设置 Idents
  if (!is.null(idents_col)) {
    Idents(sobj) <- idents_col
    if (verbose) message("  Idents set to: ", idents_col)
  }

  if (verbose) message("=== Import complete ===")
  return(sobj)
}


#' Infer Seurat reduction key from reduction name
#' @noRd
.infer_reduction_key <- function(rd_name) {
  known_keys <- list(
    pca     = "PC_",
    umap    = "UMAP_",
    tsne    = "tSNE_",
    harmony = "harmony_",
    scvi    = "scVI_",
    lsi     = "LSI_"
  )
  rd_lower <- tolower(rd_name)
  if (rd_lower %in% names(known_keys)) return(known_keys[[rd_lower]])
  return(paste0(rd_name, "_"))
}

