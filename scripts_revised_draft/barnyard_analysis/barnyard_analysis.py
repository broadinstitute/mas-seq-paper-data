#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import sys

import matplotlib.pylab as plt
import numpy as np
import pandas as pd
from time import time
import logging
import pickle
from operator import itemgetter
import json, pprint

import argparse
import yaml

import scanpy as sc
import anndata

from collections import defaultdict
from itertools import groupby
from operator import itemgetter
from typing import List, Dict, Union, Any

from time import time

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
log_info = logger.warning

import warnings
warnings.filterwarnings("ignore")

sc.settings.set_figure_params(dpi=80, facecolor='white')


# In[20]:


parser = argparse.ArgumentParser(description='Harmonize short-read and long-reads AnnData objects')

parser.add_argument(
    '-s',
    '--short',
    help='Path to harmonized short-reads AnnData objects',
    required=True)

parser.add_argument(
    '-l',
    '--long',
    help='Path to harmonized long-reads AnnData objects',
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

parser.add_argument(
    '--min-cells-per-transcript',
    default=0,
    required=False)

parser.add_argument(
    '--min-cells-per-gene',
    default=0,
    required=False)

parser.add_argument(
    '--discovery-adata',
    default='short',
    required=False)

parser.add_argument(
    '--barnyard-obs-column-name',
    default='leiden_crude',
    required=False)

parser.add_argument(
    '--lo-expr-threshold',
    default=1e-5,
    type=float,
    required=False)

parser.add_argument(
    '--hi-expr-threshold',
    default=1e-4,
    type=float,
    required=False)

parser.add_argument(
    '--lo-expr-threshold-sensitive',
    default=1e-6,
    type=float,
    required=False)

parser.add_argument(
    '--hi-expr-threshold-sensitive',
    default=1e-4,
    type=float,
    required=False)

parser.add_argument(
    '--cell-purification-threshold',
    default=1.,
    type=float,
    required=False)

parser.add_argument(
    '--contamination-threshold',
    default=0.1,
    type=float,
    required=False)

parser.add_argument(
    '--min-umi-final',
    default=50,
    type=int,
    required=False)

# # input
# short_h5_path = '/home/jupyter/mb-ml-dev-disk/MAS-seq-analysis/notebooks/cd8_t_cell_vdj_citeseq__paper__revised/scripts/test.harmonized.barnyard.short.h5ad'
# long_h5_path = '/home/jupyter/mb-ml-dev-disk/MAS-seq-analysis/notebooks/cd8_t_cell_vdj_citeseq__paper__revised/scripts/test.harmonized.barnyard.long.h5ad'
# output_prefix = 'test'
# output_path = '/home/jupyter/mb-ml-dev-disk/MAS-seq-analysis/notebooks/cd8_t_cell_vdj_citeseq__paper__revised/scripts'
# args = f'-s {short_h5_path} -l {long_h5_path} -o {output_path} -p {output_prefix}'.split(' ')

# args = vars(parser.parse_args(args))

args = vars(parser.parse_args())


# In[3]:


# input path
final_long_adata_raw_h5_path = args['long']
final_short_adata_raw_h5_path = args['short']
output_path = args['output_path']
output_prefix = args['output_prefix']

# configuration
min_cells_per_transcript = int(args['min_cells_per_transcript'])
min_cells_per_gene = int(args['min_cells_per_gene'])
group_cells_by_obs_key = args['barnyard_obs_column_name']

# consts
GENE_IDS_KEY = 'gene_ids'
GENE_NAMES_KEY = 'gene_names'

# containers
metrics = dict()


# In[4]:


log_info('Loading long-reads AnnData object...')
adata_long = sc.read(final_long_adata_raw_h5_path)
assert group_cells_by_obs_key in adata_long.obs.columns


# In[5]:


log_info('Filtering long-reads AnnData object...')
# remove genes that are lowly expressed
from collections import defaultdict
gene_id_to_tx_indices_map = defaultdict(list)
for i, gid in enumerate(adata_long.var[GENE_IDS_KEY].values):
    gene_id_to_tx_indices_map[gid].append(i)

included_gene_ids = []
tx_counts_i = np.asarray(adata_long.X.sum(0)).flatten()
for gid, tx_indices in gene_id_to_tx_indices_map.items():
    if np.sum(tx_counts_i[tx_indices]) >= min_cells_per_gene:
        included_gene_ids.append(gid)

adata_long = adata_long[:, adata_long.var[GENE_IDS_KEY].values.isin(included_gene_ids)]

# remove transcript that are very lowly expressed
sc.pp.filter_genes(adata_long, min_cells=min_cells_per_transcript)
total_umis = adata_long.X.sum()
tpm_threshold = 1_000_000 * min_cells_per_transcript / total_umis

log_info(f'Removing isoforms with TPM < {tpm_threshold:.2f}')


# In[6]:


# mapping from gene id to spanning tx icatces
from collections import defaultdict
gene_id_to_tx_indices_map = defaultdict(list)
for i, gid in enumerate(adata_long.var[GENE_IDS_KEY].values):
    gene_id_to_tx_indices_map[gid].append(i)

# useful auxiliary data structures    
gene_ids = sorted(list(gene_id_to_tx_indices_map.keys()))
n_genes = len(gene_ids)
n_transcripts = adata_long.shape[1]
gene_id_to_gene_name_map = {
    gene_id: gene_name for gene_id, gene_name in zip(adata_long.var[GENE_IDS_KEY], adata_long.var[GENE_NAMES_KEY])}
gene_name_to_gene_id_map = {
    gene_name: gene_id for gene_id, gene_name in zip(adata_long.var[GENE_IDS_KEY], adata_long.var[GENE_NAMES_KEY])}
gene_names = list(map(gene_id_to_gene_name_map.get, gene_ids))

# mapping from gene id to spanning tx indices
group_ids = adata_long.obs[group_cells_by_obs_key].values.categories.values
group_id_to_obs_indices_map = defaultdict(list)
for group_id in group_ids:
    group_id_to_obs_indices_map[group_id] = [
        idx for idx in range(len(adata_long))
        if adata_long.obs[group_cells_by_obs_key].values[idx] == group_id]


# In[7]:


import scipy

# get gene expression from isoform expression
row_indices = []
col_indices = []
values = []
for j, gene_id in enumerate(gene_ids):
    tx_indices = gene_id_to_tx_indices_map[gene_id]
    row_indices += tx_indices
    col_indices += [j] * len(tx_indices)
    values += [1] * len(tx_indices)
Y_ij = scipy.sparse.coo_matrix((values, (row_indices, col_indices)), shape=(n_transcripts, n_genes)).tocsr()
gex_X_nj = adata_long.X @ Y_ij

# normalize
adata_long_gex = sc.AnnData(
    X=gex_X_nj,
    obs=adata_long.obs,
    var=pd.DataFrame(index=pd.Index(list(map(gene_id_to_gene_name_map.get, gene_ids)))))

adata_long_gex.var_names_make_unique()


# In[8]:


log_info('Loading long-reads AnnData object...')
adata_short = sc.read(final_short_adata_raw_h5_path)


# In[9]:


adata_short_final = adata_short[:, adata_short.var.index.isin(adata_long_gex.var.index.values)]
adata_long_final = adata_long_gex[:, adata_short_final.var.index]


# In[10]:


discovery_adata_name = args['discovery_adata']
log_info(f'Identifying pure cells of each class from the discovery AnnData object ({discovery_adata_name})...')


def get_grouped_expression(adata, group_cells_by_obs_key):
    # mapping from gene id to spanning tx indices
    group_ids = adata.obs[group_cells_by_obs_key].values.categories.values
    group_id_to_obs_indices_map = defaultdict(list)
    for group_id in group_ids:
        group_id_to_obs_indices_map[group_id] = [
            idx for idx in range(len(adata))
            if adata.obs[group_cells_by_obs_key].values[idx] == group_id]
    
    n_genes = adata.shape[1]
    n_groups = len(group_id_to_obs_indices_map)
    group_expr_gi = np.zeros((n_groups, n_genes), dtype=np.int)
    for i_group, group_id in enumerate(group_ids):
        group_expr_gi[i_group, :] = np.asarray(adata.X[group_id_to_obs_indices_map[group_id], :].sum(0)).flatten()
        
    return group_expr_gi


# In[11]:


log_info('[First pass] identifying putative marker genes...')
discovery_adata = {
    'short': adata_short_final,
    'long': adata_long_final}[discovery_adata_name]

metacell_mg = get_grouped_expression(discovery_adata, group_cells_by_obs_key)
normed_metacell_mg = metacell_mg / np.sum(metacell_mg, -1, keepdims=True)
lo_expr_threshold = args['lo_expr_threshold']
hi_expr_threshold = args['hi_expr_threshold']

barnyard_gene_indices_list = []
for group_a in range(2):
    for group_b in range(2):
        lo_in_a = normed_metacell_mg[group_a, :] < lo_expr_threshold
        hi_in_b = normed_metacell_mg[group_b, :] > hi_expr_threshold
        barnyard_mask_g = lo_in_a & hi_in_b
        barnyard_gene_indices = np.where(barnyard_mask_g)[0]
        for idx in barnyard_gene_indices:
            barnyard_gene_indices_list.append((group_a, group_b, idx, discovery_adata.var.index.values[idx]))
        metrics[f'n_genes_lo_group_{group_a}_hi_group_{group_b}__first_pass'] = int(barnyard_mask_g.sum()) 
        log_info(f'[First pass] Number of genes low in group {group_a}, high in group {group_b}: {barnyard_mask_g.sum()}')


# In[12]:


putative_group_1_gene_indices = [t[2] for t in barnyard_gene_indices_list if t[0] == 0 and t[1] == 1]
putative_group_0_gene_indices = [t[2] for t in barnyard_gene_indices_list if t[0] == 1 and t[1] == 0]
putative_group_1_gene_names = [t[3] for t in barnyard_gene_indices_list if t[0] == 0 and t[1] == 1]
putative_group_0_gene_names = [t[3] for t in barnyard_gene_indices_list if t[0] == 1 and t[1] == 0]

putative_group_1_gene_expr_in_group_0_n = np.asarray(
    discovery_adata[discovery_adata.obs[group_cells_by_obs_key] == '0'][:, putative_group_1_gene_indices].X.sum(-1)).flatten()
putative_group_0_gene_expr_in_group_1_n = np.asarray(
    discovery_adata[discovery_adata.obs[group_cells_by_obs_key] == '1'][:, putative_group_0_gene_indices].X.sum(-1)).flatten()


# In[13]:


log_info('Saving putative marker genes diagnostics plot...')

fig, ax = plt.subplots()
ax.hist(np.log1p(putative_group_0_gene_expr_in_group_1_n), bins=100);
ax.set_xlabel('log1p(gex)')
ax.set_ylabel('Number of cells')
ax.set_title('Putative group 0 marker genes in group 1 cells')
fig.tight_layout()
plt.savefig(
    os.path.join(
        output_path,
        f'{output_prefix}.barnyard.putative_01_gex.png'),
    bbox_inches="tight")

fig, ax = plt.subplots()
ax.hist(np.log1p(putative_group_1_gene_expr_in_group_0_n), bins=100);
ax.set_xlabel('log1p(gex)')
ax.set_ylabel('Number of cells')
ax.set_title('Putative group 1 marker genes in group 0 cells')
fig.tight_layout()
plt.savefig(
    os.path.join(
        output_path,
        f'{output_prefix}.barnyard.putative_10_gex.png'),
    bbox_inches="tight")


# In[14]:


log_info('[First pass] Purifying cell types based on putative marker genes...')
cell_purification_threshold = args['cell_purification_threshold']

pure_group_0_cells_mask = putative_group_1_gene_expr_in_group_0_n < cell_purification_threshold
pure_group_1_cells_mask = putative_group_0_gene_expr_in_group_1_n < cell_purification_threshold
pure_group_0_cell_indices = np.where((discovery_adata.obs[group_cells_by_obs_key] == '0').values)[0][pure_group_0_cells_mask]
pure_group_1_cell_indices = np.where((discovery_adata.obs[group_cells_by_obs_key] == '1').values)[0][pure_group_1_cells_mask]
pure_both_indices = pure_group_0_cell_indices.tolist() + pure_group_1_cell_indices.tolist()

metrics[f'n_0_cells__first_pass'] = int((discovery_adata.obs[group_cells_by_obs_key] == '0').values.sum())
metrics[f'n_1_cells__first_pass'] = int((discovery_adata.obs[group_cells_by_obs_key] == '1').values.sum())
metrics[f'n_0_cells__second_passs'] = len(pure_group_0_cell_indices)
metrics[f'n_1_cells__second_pass'] = len(pure_group_1_cell_indices)


# In[15]:


adata_short_final_pure = adata_short_final[pure_both_indices]
adata_long_final_pure = adata_long_final[pure_both_indices]


# In[16]:


log_info('[Second pass] identifying highly specific marker genes...')
discovery_adata = {
    'short': adata_short_final_pure,
    'long': adata_long_final_pure}[discovery_adata_name]

metacell_mg = get_grouped_expression(adata_short_final_pure, group_cells_by_obs_key)
normed_metacell_mg = metacell_mg / np.sum(metacell_mg, -1, keepdims=True)
lo_expr_threshold = args['lo_expr_threshold_sensitive']
hi_expr_threshold = args['hi_expr_threshold_sensitive']

barnyard_gene_indices_list = []
for group_a in range(2):
    for group_b in range(2):
        lo_in_a = normed_metacell_mg[group_a, :] < lo_expr_threshold
        hi_in_b = normed_metacell_mg[group_b, :] > hi_expr_threshold
        barnyard_mask_g = lo_in_a & hi_in_b
        barnyard_gene_indices = np.where(barnyard_mask_g)[0]
        for idx in barnyard_gene_indices:
            barnyard_gene_indices_list.append((group_a, group_b, idx, discovery_adata.var.index.values[idx]))
        metrics[f'n_genes_lo_group_{group_a}_hi_group_{group_b}__second_pass'] = int(barnyard_mask_g.sum())
        log_info(f'[Second pass] Number of genes low in group {group_a}, high in group {group_b}: {barnyard_mask_g.sum()}')


# In[17]:


final_group_1_gene_indices = [t[2] for t in barnyard_gene_indices_list if t[0] == 0 and t[1] == 1]
final_group_0_gene_indices = [t[2] for t in barnyard_gene_indices_list if t[0] == 1 and t[1] == 0]
final_group_1_gene_names = [t[3] for t in barnyard_gene_indices_list if t[0] == 0 and t[1] == 1]
final_group_0_gene_names = [t[3] for t in barnyard_gene_indices_list if t[0] == 1 and t[1] == 0]


# In[18]:


fig, ax = plt.subplots(figsize=(5, 4))

adata = adata_short_final_pure.copy()

contamination_threshold = args['contamination_threshold']
min_umi_final = args['min_umi_final']

final_group_1_gene_expr_in_group_0_n = np.asarray(
    adata[adata.obs['leiden_crude'] == '0'][:, final_group_1_gene_indices].X.sum(-1)).flatten()
final_group_0_gene_expr_in_group_1_n = np.asarray(
    adata[adata.obs['leiden_crude'] == '1'][:, final_group_0_gene_indices].X.sum(-1)).flatten()
final_group_1_gene_expr_in_group_1_n = np.asarray(
    adata[adata.obs['leiden_crude'] == '1'][:, final_group_1_gene_indices].X.sum(-1)).flatten()
final_group_0_gene_expr_in_group_0_n = np.asarray(
    adata[adata.obs['leiden_crude'] == '0'][:, final_group_0_gene_indices].X.sum(-1)).flatten()

outlier_group_1_cells_n = final_group_0_gene_expr_in_group_1_n > (contamination_threshold * final_group_1_gene_expr_in_group_1_n)
outlier_group_1_cells_n = outlier_group_1_cells_n & (final_group_1_gene_expr_in_group_1_n > min_umi_final)
outlier_group_0_cells_n = final_group_1_gene_expr_in_group_0_n > (contamination_threshold * final_group_0_gene_expr_in_group_0_n)
outlier_group_0_cells_n = outlier_group_0_cells_n & (final_group_0_gene_expr_in_group_0_n > min_umi_final)

# other statistics
median_umi_per_cell = np.median(np.asarray(adata_short.X.sum(-1)).flat)
group_1_purity = 100. * final_group_1_gene_expr_in_group_1_n.sum() / (final_group_0_gene_expr_in_group_1_n.sum() + final_group_1_gene_expr_in_group_1_n.sum())
group_0_purity = 100. * final_group_0_gene_expr_in_group_0_n.sum() / (final_group_1_gene_expr_in_group_0_n.sum() + final_group_0_gene_expr_in_group_0_n.sum())
n_group_0_cells = (adata.obs['leiden_crude'] == '0').sum()
n_group_1_cells = (adata.obs['leiden_crude'] == '1').sum()

metrics['n_outlier_group_0_cells__short'] = int(np.sum(outlier_group_0_cells_n))
metrics['n_outlier_group_1_cells__short'] = int(np.sum(outlier_group_1_cells_n))
metrics['median_umi_per_cell__short'] = int(median_umi_per_cell)
metrics['group_0_purity__short'] = float(group_0_purity)
metrics['group_1_purity__short'] = float(group_1_purity)
metrics['n_group_0_cells__short'] = int(n_group_0_cells)
metrics['n_group_1_cells__short'] = int(n_group_1_cells)

# all points
ax.scatter(
    final_group_0_gene_expr_in_group_0_n,
    final_group_1_gene_expr_in_group_0_n,
    s=1,
    label=f'Group 0 (N={n_group_0_cells})')

ax.scatter(
    final_group_0_gene_expr_in_group_1_n,
    final_group_1_gene_expr_in_group_1_n,
    s=1,
    label=f'Group 1 (N={n_group_1_cells})')

# outliers
ax.scatter(
    final_group_0_gene_expr_in_group_0_n[outlier_group_0_cells_n],
    final_group_1_gene_expr_in_group_0_n[outlier_group_0_cells_n],
    s=50,
    facecolor='none',
    edgecolor='red',
    marker='o',
    lw=0.5,
    label=f'Group 1 in Group 0 > {int(100. * contamination_threshold)}% (N={outlier_group_0_cells_n.sum()})')

ax.scatter(
    final_group_0_gene_expr_in_group_1_n[outlier_group_1_cells_n],
    final_group_1_gene_expr_in_group_1_n[outlier_group_1_cells_n],
    s=50,
    facecolor='none',
    edgecolor='black',
    marker='o',
    lw=0.5,
    label=f'Group 0 in Group 1 > {int(100. * contamination_threshold)}% (N={outlier_group_1_cells_n.sum()})')


plt.plot(
    [], [], ' ',
    label=f"Median UMIs per cell: {int(median_umi_per_cell)}")
plt.plot(
    [], [], ' ',
    label=f"Group 1 purity: {group_1_purity:.1f}%")
plt.plot(
    [], [], ' ',
    label=f"Group 0 purity: {group_0_purity:.1f}%")

ax.set_xlim((-40, 2000))
ax.set_ylim((-40, 2000))

ax.set_xlabel('Group 0-specific total GEX')
ax.set_ylabel('Group 1-specific total GEX')

ax.set_title(f'Short-reads Dataset')
ax.legend(fontsize=10)

fig.tight_layout()

log_info('Saving short-reads barnyard plot...')
plt.savefig(
    os.path.join(
        output_path,
        f'{output_prefix}.barnyard.short_reads.png'),
    bbox_inches="tight")


# In[19]:


fig, ax = plt.subplots(figsize=(5, 4))

adata = adata_long_final_pure.copy()

contamination_threshold = args['contamination_threshold']
min_umi_final = args['min_umi_final']

final_group_1_gene_expr_in_group_0_n = np.asarray(
    adata[adata.obs['leiden_crude'] == '0'][:, final_group_1_gene_indices].X.sum(-1)).flatten()
final_group_0_gene_expr_in_group_1_n = np.asarray(
    adata[adata.obs['leiden_crude'] == '1'][:, final_group_0_gene_indices].X.sum(-1)).flatten()
final_group_1_gene_expr_in_group_1_n = np.asarray(
    adata[adata.obs['leiden_crude'] == '1'][:, final_group_1_gene_indices].X.sum(-1)).flatten()
final_group_0_gene_expr_in_group_0_n = np.asarray(
    adata[adata.obs['leiden_crude'] == '0'][:, final_group_0_gene_indices].X.sum(-1)).flatten()

outlier_group_1_cells_n = final_group_0_gene_expr_in_group_1_n > (contamination_threshold * final_group_1_gene_expr_in_group_1_n)
outlier_group_1_cells_n = outlier_group_1_cells_n & (final_group_1_gene_expr_in_group_1_n > min_umi_final)
outlier_group_0_cells_n = final_group_1_gene_expr_in_group_0_n > (contamination_threshold * final_group_0_gene_expr_in_group_0_n)
outlier_group_0_cells_n = outlier_group_0_cells_n & (final_group_0_gene_expr_in_group_0_n > min_umi_final)

# other statistics
median_umi_per_cell = np.median(np.asarray(adata_long.X.sum(-1)).flat)
group_1_purity = 100. * final_group_1_gene_expr_in_group_1_n.sum() / (final_group_0_gene_expr_in_group_1_n.sum() + final_group_1_gene_expr_in_group_1_n.sum())
group_0_purity = 100. * final_group_0_gene_expr_in_group_0_n.sum() / (final_group_1_gene_expr_in_group_0_n.sum() + final_group_0_gene_expr_in_group_0_n.sum())
n_group_0_cells = (adata.obs['leiden_crude'] == '0').sum()
n_group_1_cells = (adata.obs['leiden_crude'] == '1').sum()

metrics['n_outlier_group_0_cells__long'] = int(np.sum(outlier_group_0_cells_n))
metrics['n_outlier_group_1_cells__long'] = int(np.sum(outlier_group_1_cells_n))
metrics['median_umi_per_cell__long'] = int(median_umi_per_cell)
metrics['group_0_purity__long'] = float(group_0_purity)
metrics['group_1_purity__long'] = float(group_1_purity)
metrics['n_group_0_cells__long'] = int(n_group_0_cells)
metrics['n_group_1_cells__long'] = int(n_group_1_cells)

# all points
ax.scatter(
    final_group_0_gene_expr_in_group_0_n,
    final_group_1_gene_expr_in_group_0_n,
    s=1,
    label=f'Group 0 (N={n_group_0_cells})')

ax.scatter(
    final_group_0_gene_expr_in_group_1_n,
    final_group_1_gene_expr_in_group_1_n,
    s=1,
    label=f'Group 1 (N={n_group_1_cells})')

# outliers
ax.scatter(
    final_group_0_gene_expr_in_group_0_n[outlier_group_0_cells_n],
    final_group_1_gene_expr_in_group_0_n[outlier_group_0_cells_n],
    s=50,
    facecolor='none',
    edgecolor='red',
    marker='o',
    lw=0.5,
    label=f'Group 1 in Group 0 > {int(100. * contamination_threshold)}% (N={outlier_group_0_cells_n.sum()})')

ax.scatter(
    final_group_0_gene_expr_in_group_1_n[outlier_group_1_cells_n],
    final_group_1_gene_expr_in_group_1_n[outlier_group_1_cells_n],
    s=50,
    facecolor='none',
    edgecolor='black',
    marker='o',
    lw=0.5,
    label=f'Group 0 in Group 1 > {int(100. * contamination_threshold)}% (N={outlier_group_1_cells_n.sum()})')


plt.plot(
    [], [], ' ',
    label=f"Median UMIs per cell: {int(median_umi_per_cell)}")
plt.plot(
    [], [], ' ',
    label=f"Group 1 purity: {group_1_purity:.1f}%")
plt.plot(
    [], [], ' ',
    label=f"Group 0 purity: {group_0_purity:.1f}%")

# ax.set_xscale('log')
# ax.set_yscale('log')


ax.set_xlim((-20, 1000))
ax.set_ylim((-20, 1000))

ax.set_xlabel('Group 0-specific total GEX')
ax.set_ylabel('Group 1-specific total GEX')

ax.set_title(f'Long-reads Dataset ({output_prefix})')
ax.legend(fontsize=10)

fig.tight_layout()

log_info('Saving long-reads barnyard plot...')
plt.savefig(
    os.path.join(
        output_path,
        f'{output_prefix}.barnyard.long_reads.png'),
    bbox_inches="tight")


# In[ ]:


import yaml

log_info('Saving metrics')
print(os.path.join(output_path, f'{output_prefix}.barnyard.metrics.yaml'))
with open(os.path.join(output_path, f'{output_prefix}.barnyard.metrics.yaml'), 'w+') as f:
    yaml.dump(metrics, f)


# In[ ]:


log_info('Done!')


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




