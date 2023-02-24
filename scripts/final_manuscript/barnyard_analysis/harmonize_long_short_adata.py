#!/usr/bin/env python
# coding: utf-8

# In[86]:


import matplotlib.pylab as plt
import numpy as np
import pandas as pd
import os
import sys
from time import time
import logging
import pickle
from operator import itemgetter
import scanpy as sc
import argparse
import yaml

SMALL_SIZE = 12
MEDIUM_SIZE = 14
BIGGER_SIZE = 16

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

logger = logging.getLogger()
logger.setLevel(logging.INFO)
log_info = logger.warn

import warnings
warnings.filterwarnings("ignore")

sc.settings.set_figure_params(dpi=80, facecolor='white')


# In[62]:


parser = argparse.ArgumentParser(description='Harmonize short-read and long-reads AnnData objects')

parser.add_argument(
    '-s',
    '--short',
    help='Path to short-reads AnnData objects',
    required=True)

parser.add_argument(
    '-l',
    '--long',
    help='Path to long-reads AnnData objects',
    required=True)

parser.add_argument(
    '-o',
    '--output-path',
    help='Output path',
    required=True)

parser.add_argument(
    '-p',
    '--output-prefix',
    help='Output prefix',
    required=True)


# # input
# short_h5_path = '/home/jupyter/mb-ml-data-disk/MAS-seq-analysis/output/t-cell-vdj-cite-seq/M132TS_both.h5ad'
# long_h5_path = '/home/jupyter/mb-ml-data-disk/MAS-seq-analysis/data/t-cell-vdj/long/quant/v5/M132TS_MAS_15x_gene_tx_expression_count_matrix_tx_gene_counts_adata.h5ad'
# output_prefix = 'M132TS_both_new_pipeline__pilot'
# output_path = '/home/jupyter/mb-ml-dev-disk/MAS-seq-analysis/notebooks/cd8_t_cell_vdj_citeseq__paper__revised/scripts'
# args = f'-s {short_h5_path} -l {long_h5_path} -o {output_path} -p {output_prefix}'.split(' ')

args = vars(parser.parse_args())


# In[63]:


# constants
ADATA_SHORT_GENE_IDS_COL = 'gene_ids'
ADATA_LONG_GENE_IDS_COL = 'gene_ids'
KNOWN_GENE_PREFIX = 'ENSG'
DE_NOVO_GENE_PREFIX = 'MASG'

# traits to propagate from short adata to long adata
adata_short_obs_col_propagate = [
    'CD45_TotalSeqC',
    'CD45R_B220_TotalSeqC',
    'CD45RA_TotalSeqC',
    'CD45RO_TotalSeqC',
    'leiden_crude']

adata_short_obsm_col_propagate = [
    'X_pca',
    'X_umap']

# input
short_h5_path = args['short']
long_h5_path = args['long']

# output
output_prefix = args['output_prefix']
output_path = args['output_path']

# make sure the output path exists
os.makedirs(output_path, exist_ok=True)

harmonized_long_adata_h5_path = os.path.join(
    output_path, f'{output_prefix}.harmonized.barnyard.long.h5ad')
harmonized_short_adata_h5_path = os.path.join(
    output_path, f'{output_prefix}.harmonized.barnyard.short.h5ad')

harmonized_long_adata_mutual_genes_h5_path = os.path.join(
    output_path, f'{output_prefix}.harmonized.mutual_genes.barnyard.long.h5ad')
harmonized_short_adata_mutual_genes_h5_path = os.path.join(
    output_path, f'{output_prefix}.harmonized.mutual_genes.barnyard.short.h5ad')

# metrics container
metrics = dict()


# In[64]:


# load short anndata
log_info('Loading short-reads AnnData object...')
adata_short = sc.read(short_h5_path).raw.to_adata()
assert 'leiden_crude' in adata_short.obs.columns

# load long anndata
log_info('Loading long-reads AnnData object...')
adata_long = sc.read(long_h5_path)
adata_long.var_names_make_unique()
adata_long.obs = adata_long.obs.drop('Cell Barcode', axis=1)


# In[87]:


log_info('Inferring mutual barcodes and genes...')


# In[66]:


n_barcodes_short = adata_short.X.shape[0]
n_barcodes_long = adata_long.X.shape[0]
metrics['n_barcodes_short'] = int(n_barcodes_short)
metrics['n_barcodes_long'] = int(n_barcodes_long)


# In[67]:


adata_gene_info_set = set(zip(
    adata_long.var[ADATA_LONG_GENE_IDS_COL].values,
    adata_long.var['is_de_novo'].values,
    adata_long.var['is_gene_id_ambiguous'].values))

gene_ids_set = set(map(itemgetter(0), adata_gene_info_set))
n_total_genes_long = len(gene_ids_set)
n_total_genes_short = len(adata_short.var)
n_de_novo_genes_long = sum(map(lambda gene_id: gene_id.find(DE_NOVO_GENE_PREFIX) == 0, gene_ids_set))
n_gencode_genes_long = sum(map(lambda gene_id: gene_id.find(KNOWN_GENE_PREFIX) == 0, gene_ids_set))

metrics['n_total_genes_long'] = int(n_total_genes_long)
metrics['n_total_genes_short'] = int(n_total_genes_short)
metrics['n_gencode_genes_long'] = int(n_total_genes_long)
metrics['n_de_novo_genes_long'] = int(n_total_genes_long)


# In[68]:


from collections import Counter

adata_short_gene_id_set = set(adata_short.var[ADATA_SHORT_GENE_IDS_COL].values)
adata_long_gene_id_set = set(adata_long.var[ADATA_LONG_GENE_IDS_COL].values)

# drop gencode version suffix ...
drop_version = lambda entry: entry.split('.')[0] if entry.find('ENS') == 0 else entry

unversioned_adata_short_gene_id_counter = Counter([
    drop_version(entry) for entry in adata_short_gene_id_set])
unversioned_adata_long_gene_id_counter = Counter([
    drop_version(entry) for entry in adata_long_gene_id_set])

ver_unambiguous_adata_short_gene_id_list = [
    gene_id for gene_id in unversioned_adata_short_gene_id_counter.keys()
    if unversioned_adata_short_gene_id_counter[gene_id] == 1]
ver_unambiguous_adata_long_gene_id_list = [
    gene_id for gene_id in unversioned_adata_long_gene_id_counter.keys()
    if unversioned_adata_long_gene_id_counter[gene_id] == 1]

gene_id_ambiguous_adata_long_unversioned_gene_id_set = set(map(
    drop_version,
    adata_long[:, adata_long.var['is_gene_id_ambiguous']].var[ADATA_LONG_GENE_IDS_COL].values))

final_unversioned_unambiguous_mutual_gene_id_set =     set(ver_unambiguous_adata_long_gene_id_list)     .intersection(ver_unambiguous_adata_short_gene_id_list)     .difference(gene_id_ambiguous_adata_long_unversioned_gene_id_set)

metrics['n_adata_short_gene_id_set'] = len(adata_short_gene_id_set)
metrics['n_adata_long_gene_id_set'] = len(adata_long_gene_id_set)
metrics['n_ver_unambiguous_adata_short'] = len(ver_unambiguous_adata_short_gene_id_list)
metrics['n_ver_unambiguous_adata_long'] = len(ver_unambiguous_adata_long_gene_id_list)
metrics['n_gene_id_ambiguous_adata'] = len(gene_id_ambiguous_adata_long_unversioned_gene_id_set)
metrics['n_final_unversioned_unambiguous_mutual_gene_id_set'] = len(final_unversioned_unambiguous_mutual_gene_id_set)


# In[69]:


final_adata_short_mutual_keep_var_indices = [
    var_idx
    for var_idx, gene_id in enumerate(adata_short.var[ADATA_SHORT_GENE_IDS_COL])
    if drop_version(gene_id) in final_unversioned_unambiguous_mutual_gene_id_set]

final_adata_long_mutual_keep_var_indices = [
    var_idx
    for var_idx, gene_id in enumerate(adata_long.var[ADATA_LONG_GENE_IDS_COL])
    if drop_version(gene_id) in final_unversioned_unambiguous_mutual_gene_id_set]

# sort both by gene_ids
final_adata_short_mutual_keep_var_indices = sorted(
    final_adata_short_mutual_keep_var_indices,
    key=lambda idx: drop_version(adata_short.var[ADATA_SHORT_GENE_IDS_COL].values[idx]))

final_adata_long_mutual_keep_var_indices = sorted(
    final_adata_long_mutual_keep_var_indices,
    key=lambda idx: drop_version(adata_long.var[ADATA_LONG_GENE_IDS_COL].values[idx]))


# In[70]:


# subset long adata barcodes to short adata
adata_short_barcodes_set = set(adata_short.obs.index.values)
adata_long_keep_indices = []
found_barcodes_set = set()
for idx, bc in enumerate(adata_long.obs.index.values):
    if bc in adata_short_barcodes_set:
        adata_long_keep_indices.append(idx)
        found_barcodes_set.add(bc)
not_found_barcodes_set = adata_short_barcodes_set.difference(found_barcodes_set)

found_barcodes_list = sorted(list(found_barcodes_set))

adata_short_barcode_index_map = {
    bc: idx for idx, bc in enumerate(adata_short.obs.index.values)}
final_adata_short_keep_obs_indices = [
    adata_short_barcode_index_map[barcode]
    for barcode in found_barcodes_list]

adata_long_barcode_index_map = {
    bc: idx for idx, bc in enumerate(adata_long.obs.index.values)}
final_adata_long_keep_obs_indices = [
    adata_long_barcode_index_map[barcode]
    for barcode in found_barcodes_list]
final_adata_long_not_keep_obs_indices = sorted(list(set(adata_long_barcode_index_map.values()).difference(
    final_adata_long_keep_obs_indices)))


# In[71]:


# finally, slice
log_info('Generating harmonized AnnData objects...')
adata_short_mutual_barcodes = adata_short[final_adata_short_keep_obs_indices]
adata_short_mutual_barcodes_genes = adata_short_mutual_barcodes[:, final_adata_short_mutual_keep_var_indices]
adata_long_empty_barcodes = adata_long[final_adata_long_not_keep_obs_indices]
adata_long_mutual_barcodes = adata_long[final_adata_long_keep_obs_indices]
adata_long_mutual_genes = adata_long[:, final_adata_long_mutual_keep_var_indices]
adata_long_mutual_barcodes_genes = adata_long_mutual_barcodes[:, final_adata_long_mutual_keep_var_indices]
adata_long_empty_barcodes_mutual_genes = adata_long_empty_barcodes[:, final_adata_long_mutual_keep_var_indices]


# In[72]:


n_mutual_barcodes = len(final_adata_long_keep_obs_indices)
metrics['n_mutual_barcodes'] = int(n_mutual_barcodes)
metrics['pct_mutual_barcodes'] = float(100. * (n_mutual_barcodes / n_barcodes_short))

# UMI statistics
metrics['mean_umi_per_barcode_short'] = float(np.mean(
    np.asarray(adata_short.X.sum(-1)).flatten()))
metrics['median_umi_per_barcode_short'] = float(np.median(
    np.asarray(adata_short.X.sum(-1)).flatten()))

metrics['mean_umi_per_barcode_long'] = float(np.mean(
    np.asarray(adata_long.X.sum(-1)).flatten()))
metrics['median_umi_per_barcode_long'] = float(np.median(
    np.asarray(adata_long.X.sum(-1)).flatten()))

metrics['mean_umi_per_mutual_barcode_short'] = float(np.mean(
    np.asarray(adata_short_mutual_barcodes.X.sum(-1)).flatten()))
metrics['median_umi_per_mutual_barcode_short'] = float(np.median(
    np.asarray(adata_short_mutual_barcodes.X.sum(-1)).flatten()))

metrics['mean_umi_per_mutual_barcode_long'] = float(np.mean(
    np.asarray(adata_long_mutual_barcodes.X.sum(-1)).flatten()))
metrics['median_umi_per_mutual_barcode_long'] = float(np.median(
    np.asarray(adata_long_mutual_barcodes.X.sum(-1)).flatten()))

metrics['mean_umi_per_empty_barcode_long'] = float(np.mean(
    np.asarray(adata_long_empty_barcodes.X.sum(-1)).flatten()))
metrics['median_umi_per_empty_barcode_long'] = float(np.median(
    np.asarray(adata_long_empty_barcodes.X.sum(-1)).flatten()))

metrics['pct_umis_in_empty_barcodes_long'] = float(100. * adata_long_empty_barcodes.X.sum() / adata_long.X.sum())
metrics['pct_umis_in_empty_barcodes_mutual_genes_long'] = float(100. * adata_long_empty_barcodes_mutual_genes.X.sum() / adata_long_mutual_genes.X.sum())


# In[73]:


for col in adata_short_obs_col_propagate:
    try:
        adata_long_mutual_barcodes.obs[col] = adata_short_mutual_barcodes.obs[col].values.copy()
        adata_long_mutual_barcodes_genes.obs[col] = adata_short_mutual_barcodes_genes.obs[col].values.copy()
    except:
        log_info(f'WARNING: Could not propagate {col}!')
    
for col in adata_short_obsm_col_propagate:
    try:
        adata_long_mutual_barcodes.obsm[col] = adata_short_mutual_barcodes.obsm[col].copy()
        adata_long_mutual_barcodes_genes.obsm[col] = adata_short_mutual_barcodes_genes.obsm[col].copy()
    except:
        log_info(f'WARNING: Could not propagate {col}!')


# In[74]:


# save
log_info('Saving harmonized AnnData objects...')
adata_long_mutual_barcodes.write(harmonized_long_adata_h5_path)
adata_short_mutual_barcodes.write(harmonized_short_adata_h5_path)
adata_long_mutual_barcodes_genes.write(harmonized_long_adata_mutual_genes_h5_path)
adata_short_mutual_barcodes_genes.write(harmonized_short_adata_mutual_genes_h5_path)


# In[ ]:


log_info('Generating plots...')


# In[75]:


# highest expressed genes
with plt.rc_context():
    sc.pl.highest_expr_genes(adata_short_mutual_barcodes_genes, n_top=50, show=False)
    plt.savefig(
        os.path.join(
            output_path,
            f'{output_prefix}.highly_expressed_genes.short.png'),
        bbox_inches="tight")


# In[76]:


new_adata_long_index = list(f'{gs} ({teq})' for gs, teq in zip(
    adata_long_mutual_barcodes_genes.var['gene_names'],
    adata_long_mutual_barcodes_genes.var['transcript_eq_classes']))
adata_long_mutual_barcodes_genes.var['genes_names_teq'] = new_adata_long_index
adata_long_mutual_barcodes_genes.var.set_index('genes_names_teq', inplace=True, drop=True)


# In[77]:


# highest expressed isoforms
with plt.rc_context():
    sc.pl.highest_expr_genes(adata_long_mutual_barcodes_genes, n_top=50, show=False)
    plt.savefig(
        os.path.join(
            output_path,
            f'{output_prefix}.highly_expressed_genes.long.png'),
        bbox_inches="tight")


# In[78]:


total_tx_expr_long = np.asarray(adata_long_mutual_barcodes_genes.X.sum(0)).flatten()
total_gene_expr_short = np.asarray(adata_short_mutual_barcodes_genes.X.sum(0)).flatten()


# In[79]:


short_gene_ids = list(map(drop_version, adata_short_mutual_barcodes_genes.var[ADATA_SHORT_GENE_IDS_COL].values))
long_gene_ids = list(map(drop_version, adata_long_mutual_barcodes_genes.var[ADATA_LONG_GENE_IDS_COL].values))


# In[80]:


from itertools import groupby
from operator import itemgetter

gene_id_to_tx_indices_map = dict()
for g in groupby(enumerate(long_gene_ids), key=itemgetter(1)):
    gene_id = g[0]
    tx_indices = list(map(itemgetter(0), g[1]))
    gene_id_to_tx_indices_map[gene_id] = tx_indices
    
total_gene_expr_long = []
for gene_id in short_gene_ids:
    total_gene_expr_long.append(np.sum(total_tx_expr_long[gene_id_to_tx_indices_map[gene_id]]))
total_gene_expr_long = np.asarray(total_gene_expr_long)


# In[81]:


total_gene_expr_short_tpm = 1_000_000 * total_gene_expr_short / np.sum(total_gene_expr_short)
total_gene_expr_long_tpm = 1_000_000 * total_gene_expr_long / np.sum(total_gene_expr_long)


# In[82]:


import matplotlib as mpl
mpl.rcParams.update(mpl.rcParamsDefault)


# In[83]:


import matplotlib.ticker as tck
from sklearn.metrics import r2_score

fig, ax = plt.subplots(figsize=(4, 4))

ax.plot([1e-1, 1e5], [1e-1, 1e5], '--', lw=1, color='black')
ax.scatter(total_gene_expr_short_tpm, total_gene_expr_long_tpm, s=1, alpha=0.2, color='gray')
r2 = r2_score(np.log1p(total_gene_expr_short_tpm), np.log1p(total_gene_expr_long_tpm))
ax.text(0.15, 3e4, f'$R^2$ = {r2:.2f}', fontsize=10)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xticks([1e-1, 1e1, 1e3, 1e5])
ax.set_yticks([1e-1, 1e1, 1e3, 1e5])
ax.xaxis.set_minor_locator(tck.AutoMinorLocator())
ax.yaxis.set_minor_locator(tck.AutoMinorLocator())
ax.set_xlim((1e-1, 1e5))
ax.set_ylim((1e-1, 1e5))
ax.set_aspect('equal')
# ax.set_title('M132TS')
ax.set_xlabel('NGS total expression (TPM)')
ax.set_ylabel('MAS-seq total expression (TPM)')

fig.tight_layout()

plt.savefig(
    os.path.join(
        output_path,
        f'{output_prefix}.gex.short.long.concordance.png'),
    bbox_inches="tight")


# In[84]:


metrics['gex_concordance_r2'] = float(r2)


# In[85]:


log_info('Saving metrics...')
with open(os.path.join(output_path, f'{output_prefix}.metrics.yaml'), 'w+') as f:
    yaml.dump(metrics, f)


# In[ ]:


log_info('Done!')

