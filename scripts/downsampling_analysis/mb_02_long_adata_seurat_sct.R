library(Seurat)
library(SeuratData)
library(SeuratDisk)
library(ggplot2)
library(sctransform)

# the following variables must be set accordingly
seurat_input_adata_path <- "path/to/seurat/input/adata/objects/"
seurat_output_adata_path <- "path/to/seurat/output/adata/objects/"

# list of all samples to process
prefix_list <- c(
  "M132TS_immune_gencode_1m_harmonized_long",
  "M132TS_immune_gencode_5m_harmonized_long",
  "M132TS_immune_gencode_10m_harmonized_long",
  "M132TS_immune_gencode_20m_harmonized_long",
  "M132TS_immune_gencode_30m_harmonized_long",
  "M132TS_immune_gencode_isoseq_harmonized_long",
  "M132TS_immune_gencode_masseq_harmonized_long",
  "M132TS_immune_gencode_st2_ds_1m_harmonized_long",
  "M132TS_immune_gencode_st2_ds_5m_harmonized_long",
  "M132TS_immune_gencode_st2_ds_10m_harmonized_long",
  "M132TS_immune_gencode_st2_ds_20m_harmonized_long",
  "M132TS_immune_gencode_st2_ds_30m_harmonized_long",
  "M132TS_immune_gencode_st2_ds_isoseq_harmonized_long",
  "M132TS_immune_gencode_st2_ds_masseq_harmonized_long"
  "M132TS_immune_gencode_st2_full_1m_harmonized_long",
  "M132TS_immune_gencode_st2_full_5m_harmonized_long",
  "M132TS_immune_gencode_st2_full_10m_harmonized_long",
  "M132TS_immune_gencode_st2_full_20m_harmonized_long",
  "M132TS_immune_gencode_st2_full_30m_harmonized_long",
  "M132TS_immune_gencode_st2_full_isoseq_harmonized_long",
  "M132TS_immune_gencode_st2_full_masseq_harmonized_long")

for (i in 1:length(prefix_list)) {
  
  prefix <- prefix_list[i]
  input_h5ad_path <- paste(seurat_input_adata_path, prefix, "_seurat_input.h5ad", sep = "")
  input_h5seurat_path <- paste(seurat_input_adata_path, prefix, "_seurat_input.h5seurat", sep = "")
  output_h5seurat_path <- paste(seurat_output_adata_path, prefix, "_seurat_output.h5seurat", sep = "")
  
  print(output_h5seurat_path)
  
  Convert(input_h5ad_path, dest = "h5seurat", overwrite = TRUE)
  m132ts_long <- LoadH5Seurat(input_h5seurat_path)
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

}
