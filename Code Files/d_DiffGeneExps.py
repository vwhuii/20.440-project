import scipy as sp
from pylab import *
from matplotlib import pyplot as plt
import seaborn as sns
import sklearn
import pandas as pd
import scanpy as sc
import anndata as ad
import os
import scipy.stats as stats
from scipy.stats import chi2_contingency
import umap
import celltypist
from celltypist import models

from c_UMAP_MarkerGene_analysis import AnnData_UMAP
adata_nl, adata_tumor = AnnData_UMAP()
#folder to save figures
folder_path = r"C:\Users\veron\Documents\GitHub\20.440-project\Figures"

#region #### Differential Gene Expression Analysis TUMOR ######### 
#label datasets, combine datasets with labels
adata_combined = adata_nl.concatenate(adata_tumor, batch_key='condition', batch_categories=['normal', 'tumor'])
adata_combined.obs['condition_cell_group'] = adata_combined.obs['condition'].astype(str) + "_" + adata_combined.obs['cell_group'].astype(str)

# Perform differential expression analysis comparing normal vs tumor conditions
# group by cell type and run DGE per cell type for tumor vs ctl
cell_types = ['Naive CD4+ T cell','Helper CD4+ T cell', 'Regulatory T cell']
fig, axs = plt.subplots(2, 3, figsize=(15, 8))

for i, cell_type in enumerate(cell_types):
    adata_sub = adata_combined[adata_combined.obs['cell_group'] == cell_type, :]
    if adata_sub.n_obs == 0:
        continue
    sc.tl.rank_genes_groups(adata_sub, groupby='condition', method='wilcoxon')
    markers = sc.get.rank_genes_groups_df(adata_sub, group='tumor')
    markers['condition'] = 'tumor'
    markers['neg_log10_padj'] = -np.log10(markers['pvals_adj'])
    # filter for significant genes
    markers_filtered = markers[(markers.pvals_adj < 0.05) & ((markers.logfoldchanges > 0.5)
        | (markers.logfoldchanges < -0.5))]

    top_genes = markers_filtered.sort_values('pvals_adj').head(20)
    print('Top 10 genes are', top_genes['names'].tolist())

    # Create pivot table: genes as rows, groups as columns
    pivot = top_genes.pivot_table(index='names',columns='condition',
        values='logfoldchanges',aggfunc='mean')

    # Sort pivot table by log fold change (absolute value)
    pivot['abs_logfoldchange'] = pivot.abs().mean(axis=1)  # Get mean absolute logFC for sorting
    pivot = pivot.sort_values(by='abs_logfoldchange', ascending=False)
    pivot = pivot.drop(columns=['abs_logfoldchange'])  # Remove the temporary sorting column

#sort up and downregulated genes
    conditions = []
    for idx, row in markers.iterrows():
        if (row['pvals_adj'] < 0.05) and (row['logfoldchanges'] > 0.5):
            conditions.append('Upregulated')
        elif (row['pvals_adj'] < 0.05) and (row['logfoldchanges'] < -0.5):
            conditions.append('Downregulated')
        else:
            conditions.append('Not Significant')
    markers['category'] = conditions

    sns.scatterplot(data=markers, x='logfoldchanges', y='neg_log10_padj',
                  hue='category',
                  palette={'Upregulated': 'red', 'Downregulated': 'blue', 'Not Significant': 'grey'},
                    edgecolor=None, alpha=0.7, ax=axs[0,i])

    axs[0,i].axhline(-np.log10(0.05), color='black', linestyle='--', lw=1)  # Significance threshold line
    axs[0,i].axvline(0.5, color='black', linestyle='--', lw=1)   # Fold-change threshold lines
    axs[0,i].axvline(-0.5, color='black', linestyle='--', lw=1)

    axs[0,i].set_xlabel('Log2 Fold Change')
    axs[0,i].set_ylabel('-Log10 Adjusted p-value')
    axs[0,i].set_title(cell_type)
    axs[0,i].legend([],[], frameon=False)  # Hide legend for hue (optional)

#plot heatmap, shows log fold changes of the most significant genes (based off pval_adj)
    sns.heatmap(pivot, annot=False, cmap='viridis', center=0,
                     annot_kws={'size': 12, 'rotation': 0}, ax=axs[1,i])
    axs[1,i].set_aspect(0.05)
    axs[1,i].set_title(cell_type)
    axs[1, i].set_yticks(np.arange(pivot.shape[0]))
    axs[1, i].set_yticklabels(pivot.index, rotation=0, fontsize=10)
    axs[1,i].set_xticks([])
    axs[1,i].set_ylabel('genes')
    axs[1,i].set_xlabel('')

file_path = os.path.join(folder_path, "DGEA_tumor.png")
plt.savefig(file_path)
plt.tight_layout()
plt.show()

## for presentation figure

# Perform differential expression analysis comparing normal vs tumor conditions
# group by cell type and run DGE per cell type for tumor vs ctl
import matplotlib.colors as mcolors


cell_types = ['Naive CD4+ T cell','Helper CD4+ T cell', 'Regulatory T cell']
fig, axs = plt.subplots(1, 3, figsize=(10,4))

for i, cell_type in enumerate(cell_types):
    adata_sub = adata_combined[adata_combined.obs['cell_group'] == cell_type, :]
    if adata_sub.n_obs == 0:
        continue
    sc.tl.rank_genes_groups(adata_sub, groupby='condition', method='wilcoxon')
    markers = sc.get.rank_genes_groups_df(adata_sub, group='tumor')
    markers['condition'] = 'tumor'
    markers['neg_log10_padj'] = -np.log10(markers['pvals_adj'])
    # filter for significant genes
    markers_filtered = markers[(markers.pvals_adj < 0.05) & ((markers.logfoldchanges > 0.5)
        | (markers.logfoldchanges < -0.5))]

    top_genes = markers_filtered.sort_values('pvals_adj').head(20)
    print('Top 10 genes are', top_genes['names'].tolist())

    # Create pivot table: genes as rows, groups as columns
    pivot = top_genes.pivot_table(index='names',columns='condition',
        values='logfoldchanges',aggfunc='mean')

    # Sort pivot table by log fold change (absolute value)
    pivot['abs_logfoldchange'] = pivot.abs().mean(axis=1)  # Get mean absolute logFC for sorting
    pivot = pivot.sort_values(by='abs_logfoldchange', ascending=False)
    pivot = pivot.drop(columns=['abs_logfoldchange'])  # Remove the temporary sorting column

#sort up and downregulated genes
    conditions = []
    for idx, row in markers.iterrows():
        if (row['pvals_adj'] < 0.05) and (row['logfoldchanges'] > 0.5):
            conditions.append('Upregulated')
        elif (row['pvals_adj'] < 0.05) and (row['logfoldchanges'] < -0.5):
            conditions.append('Downregulated')
        else:
            conditions.append('Not Significant')
    markers['category'] = conditions

    cmap = sns.color_palette("vlag", as_cmap=True)
# Normalize around 0
    norm = mcolors.TwoSlopeNorm(vmin=-1.5, vcenter=0, vmax=1.5)

    sns.heatmap(pivot, annot=False, cmap=cmap, norm=norm,
                     annot_kws={'size': 12, 'rotation': 0}, ax=axs[i])
  #  axs[i].set_aspect(0.3)
    axs[i].set_title(cell_type)
    axs[i].set_yticks(np.arange(0.5, pivot.shape[0], 1))  # Offset by 0.5 for centering
    axs[i].set_yticklabels(pivot.index, rotation=0, fontsize=10, va='center')

    axs[i].set_xticks([])
    axs[i].set_ylabel('genes')
    axs[i].set_xlabel('')


plt.tight_layout()
file_path = os.path.join(folder_path, "DGEA_heatmap.png")
plt.savefig(file_path)
plt.show()

#endregion ################ DIFF GENE EXPRESSION ###################

#region #### DGEA WITHIN GROUPS ######### 
#normalize all values and graph Z score to compare between conditions.
#comparing EXPRESSION VALUES

condition = 'tumor'
cell_types = ['Naive CD4+ T cell', 'Helper CD4+ T cell', 'Regulatory T cell']
adata_sub = adata_combined[adata_combined.obs['condition'] == condition, :]
adata_sub = adata_sub[adata_sub.obs['cell_group'].isin(cell_types), :]

# Normalize and log1p transform
sc.pp.normalize_total(adata_sub, target_sum=1e4)
sc.pp.log1p(adata_sub)

# Get top genes from your earlier DEGs
sc.tl.rank_genes_groups(adata_sub, groupby='cell_group', groups=cell_types, method='wilcoxon')

all_dge = []
for ct in cell_types:
    df = sc.get.rank_genes_groups_df(adata_sub, group=ct)
    df['cell_group'] = ct
    df['condition'] = condition
    all_dge.append(df)
all_dge_df = pd.concat(all_dge)

top_genes = (
    all_dge_df.sort_values('pvals_adj')
    .groupby('cell_group')
    .head(10)['names']
    .unique()
)

# Average expression per cell type
mean_exp = pd.DataFrame(
    adata_sub[:, top_genes].X.toarray(),
    columns=top_genes,
    index=adata_sub.obs['cell_group']
).groupby(level=0).mean().T

# Z-score normalization across cell types
mean_exp_z = mean_exp.sub(mean_exp.mean(axis=1), axis=0)
mean_exp_z = mean_exp_z.div(mean_exp.std(axis=1), axis=0)

mean_exp_z = mean_exp_z[cell_types]

# Plot heatmap
plt.figure(figsize=(12, 10))
norm = mcolors.TwoSlopeNorm(vmin=-1.5, vcenter=0, vmax=1.5)
sns.heatmap(mean_exp_z, cmap='vlag', norm=norm, center=0, annot=True)
plt.title(f'Top 10 DE Genes per T Cell Type (Tumor) - Normalized Expression (Z-score)')
plt.xlabel('Cell Type')
plt.ylabel('Gene')
plt.tight_layout()
file_path = os.path.join(folder_path, "DGEA_cellgroup_heatmap.png")
plt.savefig(file_path)
plt.show()
#endregion ################ DGEA WITHIN GROUPS ###################
