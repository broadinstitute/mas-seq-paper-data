library(Seurat)
library(SeuratData)
library(SeuratDisk)
library(ggplot2)
library(sctransform)

# the following variables must be set accordingly
seurat_input_adata_path <- "path/to/seurat/input/adata/objects/"
seurat_output_adata_path <- "path/to/seurat/output/adata/objects/"

# sctransform of long-reads data
input_h5ad_path <- paste(seurat_input_adata_path, "M132TS_immune.v4.harmonized.long.stringtie.seurat_input.h5ad", sep = "")
input_h5seurat_path <- paste(seurat_input_adata_path, "M132TS_immune.v4.harmonized.long.stringtie.seurat_input.h5seurat", sep = "")
output_h5seurat_path <- paste(seurat_output_adata_path, "M132TS_immune.v4.harmonized.long.stringtie.seurat_output.no_mt_pct_regression.h5seurat", sep = "")
Convert(input_h5ad_path, dest = "h5seurat", overwrite = TRUE)
m132ts_long <- LoadH5Seurat(output_h5seurat_path)
m132ts_long$nCount_RNA <- colSums(x = m132ts_long, slot = "counts")
m132ts_long@assays[["RNA"]]@meta.features[["is_mito"]] <- grepl("^MT-|MTRNR", m132ts_long@assays[["RNA"]]@meta.features[["gencode_overlap_gene_names"]])
m132ts_long <- PercentageFeatureSet(m132ts_long, features = rownames(m132ts_long@assays[["RNA"]])[m132ts_long@assays[["RNA"]]@meta.features[["is_mito"]]], col.name = "percent.mt")
m132ts_long <- SCTransform(m132ts_long, verbose = TRUE, return.only.var.genes = FALSE, variable.features.n = 3000, ncells = NULL, n_genes = 10000)
m132ts_long <- RunPCA(m132ts_long, verbose = TRUE)
m132ts_long <- RunUMAP(m132ts_long, dims = 1:30, verbose = TRUE)
m132ts_long <- FindNeighbors(m132ts_long, dims = 1:30, verbose = TRUE)
m132ts_long <- FindClusters(m132ts_long, verbose = TRUE)
DimPlot(object = m132ts_long, reduction = "umap")
SaveH5Seurat(m132ts_long, filename = output_h5seurat_path)
Convert(output_h5seurat_path, dest = "h5ad")

# sctransform of short-reads data
input_h5ad_path <- paste(seurat_input_adata_path, "M132TS_immune.v4.harmonized.short.stringtie.seurat_input.h5ad", sep = "")
input_h5seurat_path <- paste(seurat_input_adata_path, "M132TS_immune.v4.harmonized.short.stringtie.seurat_input.h5seurat", sep = "")
output_h5seurat_path <- paste(seurat_output_adata_path, "M132TS_immune.v4.harmonized.short.stringtie.seurat_output.no_mt_pct_regression.h5seurat", sep = "")
Convert(input_h5ad_path, dest = "h5seurat", overwrite = TRUE)
m132ts_short <- LoadH5Seurat(input_h5seurat_path)
m132ts_short$nCount_RNA <- colSums(x = m132ts_short, slot = "counts")
m132ts_short <- PercentageFeatureSet(m132ts_short, pattern = "^MT-|MTRNR", col.name = "percent.mt")
m132ts_short <- SCTransform(m132ts_short, verbose = TRUE, return.only.var.genes = FALSE, variable.features.n = 3000, ncells = NULL, n_genes = 10000)
m132ts_short <- RunPCA(m132ts_short, verbose = TRUE)
m132ts_short <- RunUMAP(m132ts_short, dims = 1:30, verbose = TRUE)
m132ts_short <- FindNeighbors(m132ts_short, dims = 1:30, verbose = TRUE)
m132ts_short <- FindClusters(m132ts_short, verbose = TRUE)
DimPlot(object = m132ts_short, reduction = "umap")
SaveH5Seurat(m132ts_short, filename = output_h5seurat_path)
Convert(output_h5seurat_path, dest = "h5ad")
