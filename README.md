<h> scRNA-seq Analysis Pipeline using Seurat</h> <br>
<p>This repository contains an R-based workflow for analyzing Single-Cell RNA Sequencing (scRNA-seq) data using the Seurat package. The pipeline takes raw count matrices (in HDF5 format) and performs the standard workflow: Quality Control (QC), Normalization, Dimensionality Reduction (PCA), Clustering, and Visualization (UMAP).</p>
<br><br>
<h>This script performs the following steps:</h>
<p>Data Loading: Imports 10X Genomics HDF5 data.</p><br>
<p>Quality Control: Calculates mitochondrial content and filters low-quality cells.</p><br>
<p>Normalization: Log-normalizes the data.</p><br>
<p>Feature Selection: Identifies highly variable features (genes).</p><br>
<p>Scaling & PCA: Linear dimensionality reduction.</p><br>
<p>Clustering: Graph-based clustering at multiple resolutions.</p><br>
<p>Visualization: Generates UMAP projections.</p><br>
