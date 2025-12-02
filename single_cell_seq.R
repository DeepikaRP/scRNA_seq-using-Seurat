library(ggplot2)
library(tidyverse)
library(Seurat)
library(hdf5r)

#read h5 file
count_file <- Read10X_h5("count_matrix.h5")
cts <- count_file$`Gene Expression`

#create Seurat object with min 200 features occurring in at least 3 cells
seurat_obj <- CreateSeuratObject(cts, project = "scrna", min.cells = 3, min.features = 200)

#to check if there are low expressed cells or features (QC)
seurat_obj@meta.data

#to identify MT presence and create additional column for %info
seurat_obj[["mt_count"]]<-PercentageFeatureSet(seurat_obj, pattern = "^MT-")
seurat_obj@meta.data

#to visualize the features and counts of cells
VlnPlot(seurat_obj, features = c("nCount_RNA", "nFeature_RNA", "mt_count"), ncol = 3)+
  FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")+
  geom_smooth(method = 'lm')

#filter out features by creating a subset
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & mt_count < 5)

#normalization
seurat_obj <- NormalizeData(seurat_obj)

#to identify outliers, "vst" finds log mean and log variance to find outliers 
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)

#identify top10 variable features
top10 <- head(VariableFeatures(seurat_obj), 10)

#plot variable features
plot1 <- VariableFeaturePlot(seurat_obj)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

#scale object to avoid biases due to variable features
all_genes <- rownames(seurat_obj)
seurat_obj <- ScaleData(seurat_obj, features = all_genes)

#To perform PCA 
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))

#View different dimensions
print(seurat_obj[["pca"]], dims = 1:3, Features = 5)

#elbow plot to identify dimentions 
ElbowPlot(seurat_obj)

#clustering to find nearest neighbors for 1:15 PCAs
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:15)

#Clustering and refining resolution
seurat_obj <- FindClusters(seurat_obj, resolution = c(0.1, 0.3, 0.5, 0.9, 1))
View(seurat_obj@meta.data)
DimPlot(seurat_obj, group.by = "RNA_snn_res.0.5", label = TRUE)

#performing UMAP clustering
seurat_obj <- RunUMAP(seurat_obj, dims = 1:15)
DimPlot(seurat_obj, reduction = "umap", label = TRUE)