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

#region #### GSEA ######### 
tumor_expr = pd.DataFrame(adata_tumor.X.toarray(), columns=adata_tumor.var_names).mean(axis=0)
normal_expr = pd.DataFrame(adata_nl.X.toarray(), columns=adata_nl.var_names).mean(axis=0)

def running_ES(phenotype_correlation_values, genes, gene_set, p):
    # Handle empty gene sets or missing values gracefully
    phenotype_correlation_values = phenotype_correlation_values[genes]  # Ensure correct indexing
    abs_values = abs(phenotype_correlation_values)
    Nr = sum(abs_values[gene_set])  # Sum of absolute values for genes in the gene set
    running_score = [0]

    for gene in genes:
        if gene in gene_set:
            P_hit = (abs_values[gene]**p) / Nr  # Compute contribution from genes in the gene set
            running_score.append(float(running_score[-1] + P_hit))
        else:
            P_miss = 1 / (len(genes) - len(gene_set))  # Compute contribution from genes not in the gene set
            running_score.append(float(running_score[-1] - P_miss))

    running_score.pop(0)  # Remove initial 0 value
    max_ES = max(running_score, key=abs)
    return running_score, max_ES

def calculate_p_value(real_ES, null_ES):
    """
    Calculate p-value by comparing real ES to null ES values of the same sign.
    """
    if real_ES >= 0:
        # Only count positive null ES values for positive real ES
        p_val = np.sum(null_ES >= real_ES) / np.sum(null_ES >= 0)
    else:
        # Only count negative null ES values for negative real ES
        p_val = np.sum(null_ES <= real_ES) / np.sum(null_ES <= 0)
    if p_val == 0:
      p_val = 1 / 1000
    elif p_val == 1:
      p_val = 999 / 1000
    return p_val

def calculate_NES(real_ES, null_ES):
    same_sign_nulls = null_ES[null_ES >= 0] if real_ES >= 0 else null_ES[null_ES < 0]
    mean_null = np.mean(same_sign_nulls)
    std_null = np.std(same_sign_nulls)
    NES = (real_ES - mean_null) / std_null if std_null != 0 else 0
    return NES

gene_sets = {
    "activation": ["Il2ra", "Ifng", "Tnf", "Nme1", "Ifit3"],
    "exhaustion": ["Pdcd1", "Ctla4", "Lag3", "Havcr2", "Tigit"],
    "M2 polarization-driving": ["Stat6", "Pparg", "Klf4", "Irf4"]
}

n_permutations = 1000  # number of random scrambles
results_per_gene_set = {}

fig, axs = plt.subplots(1, 3, figsize=(18, 5), sharey=True)
gene_sets_ordered = list(gene_sets.keys())  # Ensure consistent order

for idx, gene_set_name in enumerate(gene_sets_ordered):
    ax = axs[idx]
    results_per_gene_set[gene_set_name] = {}

    for cell_group in ['Naive CD4+ T cell', 'Helper CD4+ T cell', 'Regulatory T cell']:
        adata_tumor_subset = adata_tumor[adata_tumor.obs['cell_group'] == cell_group]
        adata_nl_subset = adata_nl[adata_nl.obs['cell_group'] == cell_group]

        # Ensure unique var names and intersection
        adata_tumor_subset.var_names_make_unique()
        adata_nl_subset.var_names_make_unique()
        common_genes = adata_tumor_subset.var_names.intersection(adata_nl_subset.var_names)

        adata_tumor_subset = adata_tumor_subset[:, common_genes]
        adata_nl_subset = adata_nl_subset[:, common_genes]

        # Mean expression
        tumor_expr = np.array(adata_tumor_subset.X.mean(axis=0)).flatten()
        normal_expr = np.array(adata_nl_subset.X.mean(axis=0)).flatten()

        df = pd.DataFrame({
            'tumor': tumor_expr,
            'normal': normal_expr
        }, index=common_genes)

        epsilon = 1e-6
        df['log2FC'] = np.log2((df['tumor'] + epsilon) / (df['normal'] + epsilon))
        df = df.sort_values(by='log2FC', ascending=False)

        # GSEA
        gene_list = gene_sets[gene_set_name]
        running_ES_scores, max_ES = running_ES(df['log2FC'], df.index, gene_list, 1)

        null_ES = []
        for _ in range(n_permutations):
            scrambled_genes = np.random.permutation(df.index)
            _, fake_ES = running_ES(df['log2FC'], scrambled_genes, gene_list, 1)
            null_ES.append(fake_ES)
        null_ES = np.array(null_ES)

        p_value = calculate_p_value(max_ES, null_ES)
        NES = calculate_NES(max_ES, null_ES)

        # Store results
        results_per_gene_set[gene_set_name][cell_group] = {
            'ES': max_ES,
            'NES': NES,
            'p_value': p_value,
            'running_ES': running_ES_scores
        }

        # Plot
        ax.plot(running_ES_scores, label=f"{cell_group}\nES={max_ES:.2f}, NES={NES:.2f}, p={p_value:.3g}")

    ax.set_title(f"{gene_set_name} gene set")
    ax.set_xlabel("Gene Rank")
    ax.set_ylabel("Running ES")
    ax.axhline(0, color='gray', linestyle='--')
    ax.legend()

fig.suptitle("GSEA Running Enrichment Score by Gene Set and T Cell Type", fontsize=16)
plt.tight_layout()
file_path = os.path.join(folder_path, "GSEA Running Enrichment Score by Gene Set and T Cell Type.png")
plt.savefig(file_path)
plt.show()

#endregion ################ GSEA ###################

plt.show() #keeps figures open at end