# cellGeometry <img src="man/figures/logo.png" align="right" height="139" />

[![R-universe](https://zaoqu-liu.r-universe.dev/badges/cellGeometry)](https://zaoqu-liu.r-universe.dev/cellGeometry)
[![CRAN](https://www.r-pkg.org/badges/version/cellGeometry)](https://cran.r-project.org/package=cellGeometry)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

## Overview

**cellGeometry** is an R package for ultrafast deconvolution of bulk RNA-Sequencing data using single-cell RNA-Sequencing reference datasets. The package implements a novel geometric approach based on high-dimensional vector projection, enabling accurate estimation of cell type proportions in heterogeneous tissue samples.

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

## Example Analysis

### Cell Typist Immune Cell Dataset

This example uses the Cell Typist Global dataset containing 329,762 immune cells from the [CZ CELLxGENE repository](https://cellxgene.cziscience.com/collections/62ef75e4-cbea-454e-a0ce-998ec40223d3).

```r
library(zellkonverter)
library(cellGeometry)
library(AnnotationHub)

# Load data
typist <- readH5AD("celltypist_immune.h5ad", use_hdf5 = TRUE, reader = "R")
mat <- typist@assays@data$X
rownames(mat) <- rownames(typist)
meta <- typist@colData@listData

# Define cell annotations
subcl <- meta$Majority_voting_CellTypist
cellgrp <- meta$Majority_voting_CellTypist_high

# Build markers with parallelization
mk <- cellMarkers(mat, subclass = subcl, cellgroup = cellgrp, cores = 4)

# Convert Ensembl IDs to gene symbols
ah <- AnnotationHub()
ensDb <- ah[["AH113665"]]
mk <- gene2symbol(mk, ensDb)

# Visualize signature
signature_heatmap(mk)
spillover_heatmap(mk)
```

### Signature Refinement

```r
# Remove overlapping cell types
mk <- updateMarkers(mk,
    remove_subclass = c("Helper T cells", "Cytotoxic T cells")
)

# Diagnose signature quality
diagnose(mk)
```

### Simulation and Validation

```r
# Generate pseudo-bulk samples
set.seed(42)
sim_counts <- generate_samples(mk, n = 25)
sim_bulk <- simulate_bulk(mk, sim_counts)

# Deconvolute and evaluate
fit <- deconvolute(mk, sim_bulk)

# Assess performance
metrics <- metric_set(
    obs = sim_counts / rowSums(sim_counts) * 100,
    pred = fit$subclass$percent
)
print(metrics)

# Visualize results
plot_set(sim_counts, fit$subclass$output)
```

### Parameter Tuning

```r
# Optimize deconvolution parameters
tune_result <- tune_deconv(mk, sim_bulk, sim_counts,
    nsubclass = c(15, 20, 25, 30),
    noisefraction = c(0.2, 0.25, 0.3)
)

# Apply optimal parameters
summary(tune_result)
plot_tune(tune_result)
```

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
- **Issues**: [https://github.com/Zaoqu-Liu/cellGeometry/issues](https://github.com/Zaoqu-Liu/cellGeometry/issues)
