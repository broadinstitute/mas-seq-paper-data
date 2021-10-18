M132TS Main Analysis
--------------------

The notebooks in this folder generate the processed data and figures related to short-read and long-read scRNA-seq analysis of sample M132TS. The notebooks must be run in sequential order. The required inputs, paths, and the produced outputs are specified in the preamble markdown block at the beginning of each notebook. The following is a brief description of each script:

- `mb_00_short_reads_first_pass.ipynb`: A first look at the short-read sample to isolate tumor and immune cells
  
- `mb_01_short_reads_immune.ipynb`: Cell QC of the tumor component of the short-read sample
  
- `mb_02_harmonize_long_short_adata_to_seurat.ipynb`: Harmonizes short-read and long-read AnnData objects (e.g. mutual barcodes, mutual genes, etc.) and produces the input for Seurat analysis
  
- `mb_03_long_short_adata_suerat_sct.R`: Generating VST count features for both short-read and long-read dataset using Seurat's sctransform (SCT) workflow
  
- `mb_04_finalize_long_short_adata_from_seurat.ipynb`: Finalizes the short-read and long-read output from Seurat, including doublet removal, additional QC, clustering, and embedding
  
- `mb_05_isoform_clustering.ipynb`: A workflow for manual (decision-tree-based) and automated [experimental] determination of isoform assignments from genome aligned long reads
  
- `mb_06_ds_de_computation_fisher.ipynb`: Differential expression and differential splicing tests
  
- `mb_07_pseudotime_analysis.ipynb`: Diffusion Pseudotime (DPT) analysis, and a deep-dive study of PTPRPC isforms in pseutotime and associated splicing factors
  
- `mb_08_generate_additional_plots.ipynb`: Generates PTPRPC UMAP plots together with CITE-seq validation
