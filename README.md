# scBridge

> Seamless conversion between Seurat (R) and AnnData (Python) objects

[![R](https://img.shields.io/badge/R-4.2.0+-blue.svg)](https://www.r-project.org/)
[![Seurat](https://img.shields.io/badge/Seurat-4.0%2B-green.svg)](https://satijalab.org/seurat/)
[![License](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

scBridge provides reliable bidirectional conversion between Seurat (R) and AnnData (Python) single-cell objects. It bypasses known bugs in SeuratDisk and MuDataSeurat by using a safe MTX + CSV intermediate format with automatic temporary directory management.

## ✨ Features

- 🔄 **Bidirectional Conversion**: Seamlessly convert between Seurat and h5ad formats
- 🐍 **Python Integration**: Works with scanpy, pandas, numpy, and scipy
- 📦 **Zero Configuration**: Automatic temporary directory management
- 🔧 **Version Compatible**: Works with both Seurat v4 and v5
- 🛡️ **Robust**: Handles edge cases and sparse matrix conversions correctly
- 🎯 **Simple API**: Just two main functions: `seurat2py()` and `py2seurat()`

## 📦 Installation

### From Source

```r
# using github
remotes::install_github("xiaoqqjun/scBridge")
```

### Dependencies

**R packages:**
- Seurat (>= 4.0.0)
- Matrix
- utils
- reticulate

**Python modules (automatically checked):**
- scanpy
- pandas
- numpy
- scipy

Install Python dependencies:
```bash
pip install scanpy pandas numpy scipy
```

## 🚀 Quick Start

### Seurat → h5ad (Python)

```r
library(scBridge)

# Convert Seurat object to h5ad
seurat2py(seurat_obj, "output.h5ad")
```

### h5ad → Seurat (R)

```r
library(scBridge)

# Convert h5ad to Seurat object
seurat_obj <- py2seurat("input.h5ad")

# Visualize
DimPlot(seurat_obj, reduction = "umap")
```

## 📖 Detailed Usage

### One-Step Conversion

The simplest way to convert is using the one-step functions:

```r
# Seurat to h5ad
seurat2py(
  object = seurat_obj,
  file = "output.h5ad",
  slots = c("counts", "data"),  # Which assay slots to export
  assay = NULL,                  # Use default assay
  verbose = TRUE
)

# h5ad to Seurat
seurat_obj <- py2seurat(
  file = "input.h5ad",
  verbose = TRUE
)
```

### Step-by-Step Conversion

For more control, use the intermediate functions:

```r
# Step 1: Seurat → Bridge files
bridge_dir <- seurat_to_bridge(
  object = seurat_obj,
  bridge_dir = "path/to/bridge_dir",  # Optional, uses tempdir by default
  slots = c("counts", "data"),
  assay = "RNA",
  verbose = TRUE
)

# Step 2: Bridge files → h5ad
bridge_to_h5ad(
  bridge_dir = "path/to/bridge_dir",
  output_file = "output.h5ad",
  verbose = TRUE
)

# Clean up (optional)
unlink(bridge_dir, recursive = TRUE)
```

**Reverse direction:**

```r
# Step 1: h5ad → Bridge files
bridge_dir <- h5ad_to_bridge(
  h5ad_file = "input.h5ad",
  bridge_dir = "path/to/bridge_dir",
  verbose = TRUE
)

# Step 2: Bridge files → Seurat
seurat_obj <- bridge_to_seurat(
  bridge_dir = "path/to/bridge_dir",
  verbose = TRUE
)

# Clean up (optional)
unlink(bridge_dir, recursive = TRUE)
```

## 📋 What Gets Preserved?

| Data Type | Seurat → h5ad | h5ad → Seurat |
|-----------|---------------|---------------|
| **Expression Matrix** | ✅ counts, data, scale.data | ✅ X, layers |
| **Cell Metadata** | ✅ All columns | ✅ obs |
| **Gene Names** | ✅ | ✅ var_names |
| **Cell Names (Barcodes)** | ✅ | ✅ obs_names |
| **Dimensionality Reductions** | ✅ PCA, UMAP, tSNE, etc. | ✅ obsm |
| **Variable Features** | ✅ | ✅ highly_variable |

## 🔧 API Reference

### Main Functions

#### `seurat2py(object, file, ...)`
Convert Seurat object to h5ad file.

**Arguments:**
- `object`: A Seurat object
- `file`: Output h5ad file path
- `slots`: Assay slots to export (default: `c("counts", "data")`)
- `assay`: Assay name (default: default assay)
- `verbose`: Print progress messages (default: `TRUE`)

**Returns:** Output file path (invisibly)

#### `py2seurat(file, ...)`
Convert h5ad file to Seurat object.

**Arguments:**
- `file`: Input h5ad file path
- `verbose`: Print progress messages (default: `TRUE`)

**Returns:** A Seurat object

### Intermediate Functions

- `seurat_to_bridge()`: Export Seurat to intermediate files
- `bridge_to_h5ad()`: Convert intermediate files to h5ad
- `h5ad_to_bridge()`: Convert h5ad to intermediate files
- `bridge_to_seurat()`: Import intermediate files to Seurat

## 🎯 Use Cases

### Workflow 1: R → Python → R

```r
# Start with Seurat in R
library(scBridge)
seurat2py(seurat_obj, "analysis.h5ad")

# Continue in Python
# import scanpy as sc
# adata = sc.read_h5ad("analysis.h5ad")
# ... do Python analysis ...

# Back to R
seurat_obj <- py2seurat("analysis_processed.h5ad")
```

### Workflow 2: Cross-Language Collaboration

```r
# Bioinformatician (R) sends data to computational biologist (Python)
seurat2py(raw_seurat, "raw_data.h5ad")

# Computational biologist processes in Python
# ... scans, clustering, marker detection ...

# Bioinformatician receives results back
seurat_obj <- py2seurat("processed_results.h5ad")
```

## 🔍 Troubleshooting

### Missing Python Modules

If you get an error about missing Python modules:

```
Error: Missing Python modules: scanpy, pandas
Install with: pip install scanpy pandas
```

Solution:
```bash
pip install scanpy pandas numpy scipy
```

### reticulate Configuration

If Python is not found, configure reticulate:

```r
library(reticulate)

# Find Python
py_discover_config()

# Or specify Python path
use_python("C:/path/to/python.exe")
```

### Sparse Matrix Warnings

If you see deprecation warnings about sparse matrices, these are handled automatically by scBridge. The conversion will still work correctly.

### Large Files

For very large datasets:
- Ensure sufficient disk space for temporary files
- Use step-by-step conversion to inspect intermediate results
- Consider exporting only necessary assay slots

## 📊 Performance

| Dataset Size | Cells | Genes | Export Time | Import Time |
|-------------|-------|-------|-------------|-------------|
| Small | ~1K | ~20K | ~5s | ~3s |
| Medium | ~50K | ~20K | ~30s | ~20s |
| Large | ~500K | ~30K | ~5min | ~3min |

*Times may vary based on hardware and data complexity*

## 🆚 Comparison with Alternatives

| Feature | scBridge | SeuratDisk | MuDataSeurat |
|---------|----------|------------|--------------|
| V4 Support | ✅ | ✅ | ✅ |
| V5 Support | ✅ | ❌ | ❌ |
| Bidirectional | ✅ | ✅ | Limited |
| Robust Error Handling | ✅ | ❌ | ❌ |
| Zero Config | ✅ | ❌ | ❌ |
| Active Maintenance | ✅ | ⚠️ | ⚠️ |

## 🤝 Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## 📄 License

MIT License - see [LICENSE](LICENSE) file for details.

## 🙏 Acknowledgments

- Seurat team for the excellent single-cell analysis framework
- Scanpy team for Python single-cell tools
- reticulate for Python-R interoperability

## 📞 Contact

Zhijun Feng  
Email: xiaoqqjun@sina.com  
ORCID: 0000-0003-1813-1669

---

**Note:** scBridge is designed to work with Seurat v4 and v5. If you encounter any issues, please check your Seurat version and ensure Python dependencies are properly installed.
