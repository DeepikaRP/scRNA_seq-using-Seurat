<h2> scRNA-seq Analysis Pipeline using Seurat</h2> 
<p>This repository contains an R-based workflow for analyzing Single-Cell RNA Sequencing (scRNA-seq) data using the Seurat package. The pipeline takes raw count matrices (in HDF5 format) and performs the standard workflow: Quality Control (QC), Normalization, Dimensionality Reduction (PCA), Clustering, and Visualization (UMAP).</p>
<br><br>
<h2>This script performs the following steps:</h2>
<p>Data Loading: Imports 10X Genomics HDF5 data.</p>
<p>Quality Control: Calculates mitochondrial content and filters low-quality cells.</p>
<p>Normalization: Log-normalizes the data.</p>
<p>Feature Selection: Identifies highly variable features (genes).</p>
<p>Scaling & PCA: Linear dimensionality reduction.</p>
<p>Clustering: Graph-based clustering at multiple resolutions.</p>
<p>Visualization: Generates UMAP projections.</p>
