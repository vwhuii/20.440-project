#region ########### IMPORT######################
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

#endregion ################ END IMPORT###################

#region #### IMPORT DATA #########   

#normal data
base_path = r"C:\Users\veron\Documents\2025 Spr\20.440\20.440_pset6\20.440 Project Data\\"
matrix_file_nl = base_path + "matrix_normal.mtx"
features_file_nl = base_path + "features_normal.tsv"
barcodes_folder_nl = base_path + "GSM4598899_normal_gex_barcodes.tsv/"  # Assuming a folder with multiple binary files
# Load the matrix file
matrix = sp.io.mmread(matrix_file_nl).tocsc()
# Load features
features = pd.read_csv(features_file_nl, sep='\t', header=None)
features.columns = ['gene_id', 'gene_name', 'feature_type']
features.set_index('gene_name', inplace=True)
# Load barcodes from binary files
barcodes = []
for barcode_file in sorted(os.listdir(barcodes_folder_nl)):  # Ensure correct order
    with open(os.path.join(barcodes_folder_nl, barcode_file), 'rb') as f:
        barcodes.append(f.read().decode().strip())
barcodes = pd.DataFrame(barcodes, columns=['barcode'])
# Create AnnData object
adata_nl = sc.AnnData(X=matrix.T, obs=barcodes, var=features)
# Save the AnnData object
adata_nl.write(base_path + "adata_nl.h5ad")

#TUMOR DATA
#create an anndata object using the data from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152022
base_path = r"C:\Users\veron\Documents\2025 Spr\20.440\20.440_pset6\20.440 Project Data\\"
matrix_file = base_path + "matrix.mtx"
features_file = base_path + "features.tsv"
barcodes_folder = base_path + "GSM4598898_tumor_gex_barcodes.tsv/"  # Assuming a folder with multiple binary files

# Load the matrix file
matrix = sp.io.mmread(matrix_file).tocsc()
# Load features
features = pd.read_csv(features_file, sep='\t', header=None)
features.columns = ['gene_id', 'gene_name', 'feature_type']  # Adjust column names as needed
features.set_index('gene_name', inplace=True)
# Load barcodes from binary files
barcodes = []
for barcode_file in sorted(os.listdir(barcodes_folder)):  # Ensure correct order
    with open(os.path.join(barcodes_folder, barcode_file), 'rb') as f:
        barcodes.append(f.read().decode().strip())
barcodes = pd.DataFrame(barcodes, columns=['barcode'])
# Create AnnData object
adata_tumor = sc.AnnData(X=matrix.T, obs=barcodes, var=features)
# Save the AnnData object
adata_tumor.write(base_path + "adata.h5ad")

print("AnnData object created and saved successfully!")
#endregion ################ END DATA IMPORT###################

#region #### PROCESSING DATA #########   
#normal data
adata_nl.raw = adata_nl.copy() #make a copy of the raw data
sc.pp.filter_cells(adata_nl, min_genes=200)  # Remove low-quality cells
sc.pp.filter_genes(adata_nl, min_cells=3)  # Remove rarely expressed genes
sc.pp.normalize_total(adata_nl, target_sum=1e4)  # Normalize counts per cell
sc.pp.log1p(adata_nl)  # Apply log transformation
sc.pp.highly_variable_genes(adata_nl, flavor="seurat", n_top_genes=2000) #top 2000 highly var genes
sc.pl.highly_variable_genes(adata_nl)

#add qc metrics to adata
sc.pp.calculate_qc_metrics(adata_nl, inplace=True)

#label all mitochondrial genes as 'true' for mitochondrial
adata_nl.var["mt"] = adata_nl.var.index.str.startswith("mt-")

#tumor data
#filtering and preprocessing of tumor data
adata_tumor.raw = adata_tumor.copy() #make a copy of the raw data
sc.pp.filter_cells(adata_tumor, min_genes=200)  # Remove low-quality cells
sc.pp.filter_genes(adata_tumor, min_cells=3)  # Remove rarely expressed genes
sc.pp.normalize_total(adata_tumor, target_sum=1e4)  # Normalize counts per cell
sc.pp.log1p(adata_tumor)  # Apply log transformation
sc.pp.highly_variable_genes(adata_tumor, flavor="seurat", n_top_genes=2000) #top 2000 highly var genes
sc.pl.highly_variable_genes(adata_tumor)

#add qc metrics to adata
sc.pp.calculate_qc_metrics(adata_tumor, inplace=True)

#label all mitochondrial genes as 'true' for mitochondrial
adata_tumor.var["mt"] = adata_tumor.var.index.str.startswith("mt-")

#endregion ################ END PROCESS DATA###################


#region #### NORMAL CLUSTERING AND UMAP #########   
#normal
#perform PCA
sc.pp.pca(adata_nl)
sc.pl.pca_variance_ratio(adata_nl, log=True)
sc.tl.pca(adata_nl, n_comps=50)

# Compute neighbors if not already done
sc.pp.neighbors(adata_nl, n_neighbors=15)

# Run leiden clustering with higher resolution for more clustering
#sc.tl.leiden(adata_nl, resolution=1.5)
sc.tl.leiden(adata_nl, resolution=1.5)

# Plot UMAP with leiden clusters
sc.tl.umap(adata_nl, n_components=2)
sc.pl.umap(adata_nl, color='leiden')

leiden_nl = adata_nl.obs['leiden'].values

#Within each cluster, find the gene that is most expressed.
#used to identify marker genes in each cluster
sc.tl.rank_genes_groups(adata_nl, groupby='leiden', method='wilcoxon')
# View the top 20 marker genes for each cluster
sc.pl.rank_genes_groups(adata_nl, n_genes=20, sharey=False)
top_genes_indices = adata_nl.uns['rank_genes_groups']['names']

# Using zip to transpose the lists (create 15 lists with 20 items each)
transposed_lists = list(zip(*top_genes_indices))

# Convert the tuple output of zip into lists
transposed_lists = [list(t) for t in transposed_lists]

# Output the result
for i, transposed_list in enumerate(transposed_lists):
    print(f"List {i+1}: {transposed_list}")

#assign T cell groups to clusters
#create a dictionary with the cluster to cell type
cluster_to_cell_nl = {
    0: 'Naive CD4+ T cell',
    1: 'Memory CD4+ T cell',
    2:'Th17',
    3: 'Helper CD4+ T cell',
    4: 'Helper CD4+ T cell',
    5: 'Helper CD4+ T cell',
    6: 'Memory CD4+ T cell',
    7: 'Th2',
    8: 'Helper CD4+ T cell',
    9: 'Naive CD4+ T cell',
    10: 'Naive CD4+ T cell',
    11:'Helper CD4+ T cell',
    12: 'Regulatory T cell',
    13: 'Regulatory T cell',
    14: 'Th1',
    15: 'Regulatory T cell',
    16: 'Th2',
    17: 'Helper CD4+ T cell',
    18: 'T Follicular Helper Cell',
    19: 'Other'
}

#add column to adata_nl
adata_nl.obs['leiden'] = adata_nl.obs['leiden'].astype(int)
adata_nl.obs['cell_type'] = adata_nl.obs['leiden'].map(cluster_to_cell_nl)
# sc.pl.umap(adata_nl, color='cell_type')
sc.pp.neighbors(adata_nl, n_neighbors=15)

# Run leiden clustering with higher resolution for more clustering
sc.tl.leiden(adata_nl, resolution=1.5)

# Plot UMAP with leiden clusters
sc.tl.umap(adata_nl, n_components=2)
sc.pl.umap(adata_nl, color='leiden', show=True)
sc.pl.umap(adata_nl, color='cell_type', show=True,title='Normal UMAP Cell Type Clustering')
#endregion ################ END CLUSTERING AND UMAP###################

#visualize most expressed gene in that cluster
#plot top most expressed genes
top_genes = []
for i in transposed_lists:
  top_genes.extend([i[0]])
print(top_genes)
#sc.pl.umap(adata_nl, color=top_genes, color_map='plasma')

#Assign the most expressed gene based on max expression in each cluster to each cell
adata_nl.obs['leiden'] = adata_nl.obs['leiden'].astype(int)
adata_nl.obs['most_expressed_gene_max'] = adata_nl.obs['leiden'].map(lambda x: top_genes[x])

#Plot UMAP colored by the most expressed gene for each cluster (based on max expression)
sc.pl.umap(adata_nl, color='most_expressed_gene_max', title="Max Gene Expressed Per Cluster- Normal Lung")

#umap for the dominant TF for a CD4 subtype
adata_nl = adata_nl.copy()

# Compute neighbors if not already done
sc.pp.neighbors(adata_nl)

# Run leiden clustering
sc.tl.leiden(adata_nl, resolution=0.8)

#Plot plain UMAP
sc.pl.umap(adata_nl, color = 'leiden', title = 'Low Cluster Resolution Normal Lung UMAP')

# Plot UMAP with leiden clusters
sc.pl.umap(adata_nl, color=['Tbx21','Gata3','Rorc','Foxp3','Bcl6', 'Cd3d'], cmap = 'plasma')

# Violin plots
adata_nl.var_names_make_unique()
adata_nl.raw = adata_nl
sc.pl.violin(adata_nl,
             keys=['Tbx21','Gata3','Rorc','Foxp3','Bcl6', 'Cd3d'],
             groupby="leiden",
             jitter=0.4,
             scale="width",
             multi_panel=True)



#TUMOR PCA reduction
adata_tumor = adata_tumor.copy()

#perform PCA
sc.pp.pca(adata_tumor)
sc.pl.pca_variance_ratio(adata_tumor, log=True)
sc.tl.pca(adata_tumor, n_comps=50)

# Compute neighbors if not already done
sc.pp.neighbors(adata_tumor, n_neighbors=24)

# Run leiden clustering with higher resolution for more clustering
sc.tl.leiden(adata_tumor, resolution=1.5)

# Plot UMAP with leiden clusters
sc.tl.umap(adata_tumor, n_components=2)
sc.pl.umap(adata_tumor, color='leiden')
#endregion ################ END NORMAL UMAP###################

#region #### TUMOR CLUSTERING AND UMAP #########  
#TUMOR PCA reduction
adata_tumor = adata_tumor.copy()

#perform PCA
sc.pp.pca(adata_tumor)
sc.pl.pca_variance_ratio(adata_tumor, log=True)
sc.tl.pca(adata_tumor, n_comps=50)

# Compute neighbors if not already done
sc.pp.neighbors(adata_tumor, n_neighbors=24)

# Run leiden clustering with higher resolution for more clustering
sc.tl.leiden(adata_tumor, resolution=1.5)

# Plot UMAP with leiden clusters
sc.tl.umap(adata_tumor, n_components=2)
sc.pl.umap(adata_tumor, color='leiden')

#Within each cluster, find the gene that is most expressed.
#used to identify marker genes in each cluster
sc.tl.rank_genes_groups(adata_tumor, groupby='leiden', method='wilcoxon')
# View the top 20 marker genes for each cluster
sc.pl.rank_genes_groups(adata_tumor, n_genes=20, sharey=False)
top_genes_indices_tumor = adata_tumor.uns['rank_genes_groups']['names']

# Using zip to transpose the lists (create 15 lists with 20 items each)
transposed_lists_tumor = list(zip(*top_genes_indices_tumor))

# Convert the tuple output of zip into lists
transposed_lists_tumor = [list(t) for t in transposed_lists_tumor]

# Output the result
for i, transposed_list in enumerate(transposed_lists_tumor):
    print(f"List {i+1}: {transposed_list}")


    #create a dictionary with the cluster to cell type
cluster_to_cell_tumor = {
    0: 'Regulatory T cell',
    1: 'Helper CD4+ T cell',
    2:'T Follicular Helper Cell',
    3: 'Th17',
    4: 'Naive CD4+ T cell',
    5: 'Helper CD4+ T cell',
    6: 'Memory CD4+ T cell',
    7: 'Th1',
    8: 'Naive CD4+ T cell',
    9: 'Memory CD4+ T cell',
    10: 'Regulatory T cell',
    11:'Regulatory T cell'
}

#add column to adata_nl
adata_tumor.obs['leiden'] = adata_tumor.obs['leiden'].astype(int)
adata_tumor.obs['cell_type'] = adata_tumor.obs['leiden'].map(cluster_to_cell_tumor)
sc.pp.neighbors(adata_tumor, n_neighbors=24)

# Run leiden clustering with higher resolution for more clustering
sc.tl.leiden(adata_tumor, resolution=1.5)

# Plot UMAP with leiden clusters
sc.tl.umap(adata_tumor, n_components=2)
sc.pl.umap(adata_tumor, color='leiden', show=True)
sc.pl.umap(adata_tumor, color='cell_type', show=True,title='Tumor UMAP Cell Type Clustering')

#plot top most expressed genes
top_genes = []
for i in transposed_lists:
  top_genes.extend([i[0]])
print(top_genes)

#Assign the most expressed gene based on max expression in each cluster to each cell
adata_tumor.obs['leiden'] = adata_tumor.obs['leiden'].astype(int)
adata_tumor.obs['most_expressed_gene_max'] = adata_tumor.obs['leiden'].map(lambda x: top_genes[x])

#Plot UMAP colored by the most expressed gene for each cluster (based on max expression)
sc.pl.umap(adata_tumor, color='most_expressed_gene_max', title="Max Gene Expressed Per Cluster- Tumor Lung")

adata_tumor = adata_tumor.copy()

# Compute neighbors if not already done
sc.pp.neighbors(adata_tumor)

# Run leiden clustering
sc.tl.leiden(adata_tumor, resolution=0.8)

#Plot plain UMAP
sc.pl.umap(adata_tumor, color = 'leiden', title = 'Low Cluster Resolution Tumor Lung UMAP')

# Plot UMAP with leiden clusters
sc.pl.umap(adata_tumor, color=['Tbx21','Gata3','Rorc','Foxp3','Bcl6', 'Cd3d'], cmap = 'plasma')

# Violin plots
adata_tumor.var_names_make_unique()
adata_tumor.raw = adata_tumor
sc.pl.violin(adata_tumor,
             keys=['Tbx21','Gata3','Rorc','Foxp3','Bcl6', 'Cd3d'],
             groupby="leiden",
             jitter=0.4,
             scale="width",
             multi_panel=True)

#endregion ################ END TUMOR UMAP###################

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
print('done running!')