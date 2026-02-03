# cellGeometry <img src="man/figures/logo.png" align="right" height="139" alt="" />

<!-- badges: start -->
[![Documentation](https://img.shields.io/badge/docs-pkgdown-blue.svg)](https://zaoqu-liu.github.io/cellGeometry/)
[![R-universe](https://zaoqu-liu.r-universe.dev/badges/cellGeometry)](https://zaoqu-liu.r-universe.dev/cellGeometry)
[![CRAN](https://www.r-pkg.org/badges/version/cellGeometry)](https://cran.r-project.org/package=cellGeometry)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
<!-- badges: end -->

## Overview

**cellGeometry** is an R package for ultrafast deconvolution of bulk RNA-Sequencing data using single-cell RNA-Sequencing reference datasets. The package implements a novel geometric approach based on high-dimensional vector projection, enabling accurate estimation of cell type proportions in heterogeneous tissue samples.

📖 **Documentation**: https://zaoqu-liu.github.io/cellGeometry/

### Key Features

- **High-dimensional geometric deconvolution**: Utilizes vector dot product projection in gene expression space
- **Spillover compensation**: Automatically adjusts for cross-talk between similar cell types
- **Scalable architecture**: Efficiently handles large-scale scRNA-Seq references (>1M cells) via HDF5 backend
- **Multi-pass outlier detection**: Robust estimation through iterative removal of problematic genes
- **Comprehensive diagnostics**: Built-in tools for signature quality assessment and result validation

## Installation

### From R-universe (Recommended)

```r
install.packages("cellGeometry", repos = "https://zaoqu-liu.r-universe.dev")
```

### From GitHub

```r
# Install devtools if not available
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")

devtools::install_github("Zaoqu-Liu/cellGeometry")
```

### Dependencies

Bioconductor (version ≥3.20) packages are required:

```r
# Core dependencies
BiocManager::install(c("ensembldb", "DelayedArray"))

# For HDF5/h5ad file support (recommended for large datasets)
BiocManager::install(c("zellkonverter", "rhdf5", "HDF5Array"))

# For gene ID conversion
BiocManager::install("AnnotationHub")

# Optional: Seurat integration
install.packages("Seurat")
```

## Methodology

### Algorithm Overview

The cellGeometry algorithm operates in two stages:

**Stage 1: Marker Gene Selection**

For each cell type $c$, genes are represented as vectors in a $k$-dimensional space where $k$ is the number of cell clusters. The optimal markers are selected based on the angular distance from each coordinate axis:

$$\theta_c = \arccos\left(\frac{\mathbf{g} \cdot \mathbf{e}_c}{\|\mathbf{g}\|}\right)$$

where $\mathbf{g}$ is the gene expression vector normalized to the unit hypersphere and $\mathbf{e}_c$ is the unit vector for cell type $c$.

**Stage 2: Deconvolution via Vector Projection**

Bulk RNA-Seq samples are deconvoluted by computing the vector projection onto each cell type signature:

$$\text{proj}_{\mathbf{v}}(\mathbf{u}) = \frac{\mathbf{u} \cdot \mathbf{v}}{\|\mathbf{v}\|^2}$$

A compensation matrix $\mathbf{C}$ is applied to adjust for spillover between similar cell types:

$$\hat{\mathbf{p}} = \mathbf{C}^{-1} \cdot \text{proj}(\mathbf{b})$$

where $\hat{\mathbf{p}}$ represents the estimated cell type proportions and $\mathbf{b}$ is the bulk expression profile.

For detailed mathematical foundations, see the [Algorithm Vignette](https://zaoqu-liu.github.io/cellGeometry/articles/algorithm.html).

## Quick Start

### Basic Workflow

```r
library(cellGeometry)

# Load single-cell reference data
# mat: count matrix (genes × cells)
# subcl: cell subclass annotations
# cellgrp: broader cell group annotations

# Stage 1: Build cell markers
mk <- cellMarkers(mat, subclass = subcl, cellgroup = cellgrp)

# Stage 2: Deconvolute bulk RNA-Seq
fit <- deconvolute(mk, bulk_data)

# Extract results
proportions <- fit$subclass$percent
```

### Working with Large Datasets

For datasets exceeding 1 million cells, we recommend using the HDF5 backend:

```r
library(zellkonverter)

# Load h5ad file (data remains on disk)
sce <- readH5AD("large_dataset.h5ad", use_hdf5 = TRUE, reader = "R")

# Extract count matrix and metadata
mat <- sce@assays@data$X
rownames(mat) <- rownames(sce)
meta <- sce@colData@listData

# Proceed with cellMarkers (supports DelayedMatrix)
mk <- cellMarkers(mat, subclass = meta$cell_type, cores = 4)
```

## Vignettes

- [Quick Start Guide](https://zaoqu-liu.github.io/cellGeometry/articles/intro.html) - Getting started with cellGeometry
- [Mathematical Foundations](https://zaoqu-liu.github.io/cellGeometry/articles/algorithm.html) - Detailed algorithm explanation
- [Visualization Guide](https://zaoqu-liu.github.io/cellGeometry/articles/visualization.html) - Plotting and interpretation
- [Real-World Application](https://zaoqu-liu.github.io/cellGeometry/articles/real_world.html) - Complete analysis workflow
- [Brain Atlas Example](https://zaoqu-liu.github.io/cellGeometry/articles/brain_atlas.html) - Large-scale deconvolution

## Key Functions

| Function | Description |
|----------|-------------|
| `cellMarkers()` | Build cell type signatures from scRNA-Seq reference |
| `deconvolute()` | Perform bulk RNA-Seq deconvolution |
| `updateMarkers()` | Refine signature parameters |
| `tune_deconv()` | Optimize deconvolution settings |
| `mergeMarkers()` | Combine multiple scRNA-Seq references |
| `gene2symbol()` | Convert Ensembl IDs to gene symbols |
| `signature_heatmap()` | Visualize gene signature matrix |
| `spillover_heatmap()` | Visualize cell type cross-talk |
| `cos_similarity()` | Compute signature similarity matrix |
| `diagnose()` | Assess signature quality |

## Citation

If you use cellGeometry in your research, please cite:

> Liu Z, Lewis MJ. cellGeometry: Geometric single-cell deconvolution of bulk RNA-Sequencing data. R package version 0.5.8. https://github.com/Zaoqu-Liu/cellGeometry

## Contributing

Contributions are welcome! Please submit issues and pull requests on [GitHub](https://github.com/Zaoqu-Liu/cellGeometry).

## License

This package is licensed under the [GPL-3.0 License](https://www.gnu.org/licenses/gpl-3.0.en.html).

## Contact

- **Maintainer**: Zaoqu Liu (liuzaoqu@163.com)
- **GitHub**: [https://github.com/Zaoqu-Liu/cellGeometry](https://github.com/Zaoqu-Liu/cellGeometry)
- **Documentation**: [https://zaoqu-liu.github.io/cellGeometry/](https://zaoqu-liu.github.io/cellGeometry/)
- **Issues**: [https://github.com/Zaoqu-Liu/cellGeometry/issues](https://github.com/Zaoqu-Liu/cellGeometry/issues)
