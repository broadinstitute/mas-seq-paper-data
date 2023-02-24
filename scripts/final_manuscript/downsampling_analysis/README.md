M132TS Downsamplig Analysis
---------------------------

The notebooks in this folder generate the data and figures for the _in silico_ downsampling analysis of M132TS sample.

The downsampled long-reads AnnData objects are specified in `downsampling_series_sample_metadata.yaml`.

To reproduce the results, run the following notebook on all sample keys to generate input files for Seurat:
- mb_01_harmonize_long_short_to_seurat.ipynb

Next, run the following R script to perform sctransform all samples:
- mb_02_long_adata_seurat_sct.R

Next, run the following three notebooks on all sample keys:
- mb_03_finalize_long_short_adata_from_seurat.ipynb
- mb_04_ds_de_computation_fisher.ipynb
- mb_05_umi_stats.ipynb

Next, run the following notebook to generate the figures:
- mb_06_generate_downsampling_analysis_figures.ipynb

Finally, run the following notebook to produce Supplementary Figure 15:
- mb_07_isoform_assignment_specificity_analysis.ipynb
