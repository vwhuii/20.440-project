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

from c_UMAP_analysis import AnnData_UMAP
adata_nl, adata_tumor = AnnData()

#region #### CELL-TYPE FREQ ANALYSIS ######### 
#label datasets, combine datasets with labels
adata_combined = adata_nl.concatenate(adata_tumor, batch_key='condition', batch_categories=['normal', 'tumor'])
adata_combined.obs['condition_cell_type'] = adata_combined.obs['condition'].astype(str) + "_" + adata_combined.obs['cell_type'].astype(str)

#group total number of cells in a given subset
num_cells = adata_combined.obs.groupby(['condition','cell_type']).count()
num_cells_total = adata_combined.obs.groupby(['condition']).count()

# Get total cell counts for 'normal' and 'tumor'
normal_count = num_cells_total.loc['normal', 'barcode']
tumor_count = num_cells_total.loc['tumor', 'barcode']

#make a dict with 'Condition' and count
cell_count_total = {'normal': normal_count,
                    'tumor': tumor_count
                    }

#cell_count_total = pd.DataFrame.from_dict(cell_count_total, orient='index', columns=['total_cells'])

#total number of cells per condition
cell_type_counts = adata_combined.obs.groupby(['condition', 'cell_type']).count().reset_index()
cell_type_counts['total_cells']= cell_type_counts['condition'].map(cell_count_total).astype(int)
cell_type_counts['frequency'] = cell_type_counts.barcode / cell_type_counts.total_cells

#calculate frequency: count/total count
plt.figure(figsize = (10,4))
ax = sns.barplot(data = cell_type_counts, x = 'cell_type', y = 'frequency', hue = 'condition')
plt.xticks(rotation = 35, rotation_mode = 'anchor', ha = 'right')
plt.title('Frequency of CD4+ T cell Subsets')
plt.show()

#endregion ################ CELL FREQ ANALYSIS ###################

#region #### DIFFERENTIAL GENE EXPRESSION ANALYSIS ######### 
# Perform differential expression analysis comparing normal vs tumor conditions
sc.tl.rank_genes_groups(adata_combined, groupby='condition_cell_type', method='wilcoxon')
markers = sc.get.rank_genes_groups_df(adata_combined, group=None)
#filter for significant genes
markers_filtered = markers[(markers.pvals_adj < 0.05) & (markers.logfoldchanges > 0.5)]

#plot top 20 differentially expressed genes FOR TUMOR V NORMAL
top_genes = markers_filtered['names'].head(20)
top_gene_expression = markers_filtered['logfoldchanges'].head(20)
top_gene_expression_df = pd.DataFrame({
    'Gene': top_genes,
    'LogFoldChange': top_gene_expression
})
plt.figure(figsize=(10, 5))
top_gene_expression_df.set_index('Gene', inplace=True)
sns.heatmap(top_gene_expression_df.T, annot=True, cmap='viridis', annot_kws={'size': 13, 'rotation':90}, cbar_kws={'label': 'Log Fold Change'})
plt.title("Top 20 Differentially Expressed Genes - Tumor vs Normal")
plt.xlabel("Genes")

#endregion ################ DIFFERENTIAL GENE EXPRESSION ANALYSIS ###################


#region #### GSEA ######### 
tumor_expr = pd.DataFrame(adata_tumor.X.toarray(), columns=adata_tumor.var_names).mean(axis=0)
normal_expr = pd.DataFrame(adata_nl.X.toarray(), columns=adata_nl.var_names).mean(axis=0)

df = pd.DataFrame({
    'tumor': tumor_expr,
    'normal': normal_expr
})

# Step 2: Calculate log2 fold change (log2FC)
epsilon = 1e-6
df['log2FC'] = np.log2((df['tumor'] + epsilon) / (df['normal'] + epsilon))

# Step 3: Sort by log2FC values
df = df.sort_values(by='log2FC', ascending=False)

# Step 4: Define gene sets
gene_sets = {
    "activation": ["Il2ra", "Ifng", "Tnf", "Nme1", "Ifit3"],
    "exhaustion": ["Pdcd1", "Ctla4", "Lag3", "Havcr2", "Tigit"],
    "M2 polarization-driving": ["Stat6", "Pparg", "Klf4", "Irf4"]
}

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
    return running_score, max(running_score, key=abs)

# Step 5: Compute enrichment scores for each gene set
results = {}
for name, genes in gene_sets.items():
    running_ES_scores, max_ES = running_ES(df['log2FC'], df.index, genes, 1)
    results[name] = max_ES

    # Plot the running score for each gene set
    plt.figure(figsize=(10, 6))
    plt.plot(df.index, running_ES_scores, label=f'Running ES for {name}', color='blue')
    plt.axhline(0, color='gray', linestyle='--')  # Add a line at zero for reference
    plt.title(f'Running Enrichment Score (ES) for Gene Set: {name}')
    plt.xlabel('Gene Rank')
    plt.ylabel('Running Enrichment Score (ES)')
    plt.legend()
    plt.gca().get_xaxis().set_visible(False) #hides x axis, it is so many genes
    plt.tight_layout()  # Ensure the plot fits within the figure
    plt.show()

# Print the results for each gene set
print(results)

#endregion ################ GSEA ###################

plt.show() #keeps figures open at end