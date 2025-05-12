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
from matplotlib.colors import LinearSegmentedColormap

from b_Data_processing import AnnData_Processed
adata_nl, adata_tumor = AnnData_Processed()

#folder to save figures
folder_path = r"C:\Users\veron\Documents\GitHub\20.440-project\Figures"

#region #### NORMAL CLUSTERING AND UMAP #########   
#normal
#perform PCA
sc.pp.pca(adata_nl)
sc.pl.pca_variance_ratio(adata_nl, log=True, show=False)
plt.show(block=False)
sc.tl.pca(adata_nl, n_comps=50)

# Compute neighbors if not already done
sc.pp.neighbors(adata_nl, n_neighbors=15)

# Run louvain clustering with higher resolution for more clustering
sc.tl.louvain(adata_nl, resolution=1.5)

# Plot UMAP with louvain clusters
sc.tl.umap(adata_nl, n_components=2)
sc.pl.umap(adata_nl, color='louvain', show=False)
plt.show(block=False)

            
louvain_nl = adata_nl.obs['louvain'].values

#Within each cluster, find the gene that is most expressed.
#used to identify marker genes in each cluster
sc.tl.rank_genes_groups(adata_nl, groupby='louvain', method='wilcoxon')
# View the top 20 marker genes for each cluster
sc.pl.rank_genes_groups(adata_nl, n_genes=20, sharey=False, show=False)
top_genes_indices = adata_nl.uns['rank_genes_groups']['names']
plt.show(block=False)

# Using zip to transpose the lists (create 15 lists with 20 items each)
transposed_lists = list(zip(*top_genes_indices))

# Convert the tuple output of zip into lists
transposed_lists = [list(t) for t in transposed_lists]

# Output the result
#for i, transposed_list in enumerate(transposed_lists):
#    print(f"List {i+1}: {transposed_list}")

#assign T cell groups to clusters
#create a dictionary with the cluster to cell type
cluster_to_cell_nl = {
    0: 'Naive CD4+ T cell',
    1: 'Helper CD4+ T cell', #m
    2:'Th17',
    3: 'Helper CD4+ T cell',
    4: 'Helper CD4+ T cell',
    5: 'Helper CD4+ T cell',
    6: 'Helper CD4+ T cell', #m
    7: 'Helper CD4+ T cell',
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
adata_nl.obs['louvain'] = adata_nl.obs['louvain'].astype(int)
adata_nl.obs['cell_type'] = adata_nl.obs['louvain'].map(cluster_to_cell_nl)


# Let's say you want to remove clusters "3" and "5"
clusters_to_remove = ['19']

# Filter out the unwanted cluster
adata_nl = adata_nl[~adata_nl.obs['louvain'].isin(clusters_to_remove)].copy()


##create another classification for analysis groups, with category called "cell_group"
#create a dictionary with the cluster to cell type
cluster_to_cell_nl = {
    0: 'Naive CD4+ T cell',
    1: 'Helper CD4+ T cell', #m
    2:'Helper CD4+ T cell',
    3: 'Helper CD4+ T cell',
    4: 'Helper CD4+ T cell',
    5: 'Helper CD4+ T cell',
    6: 'Helper CD4+ T cell', #m
    7: 'Helper CD4+ T cell',
    8: 'Helper CD4+ T cell',
    9: 'Naive CD4+ T cell',
    10: 'Naive CD4+ T cell',
    11:'Helper CD4+ T cell',
    12: 'Regulatory T cell',
    13: 'Regulatory T cell',
    14: 'Helper CD4+ T cell',
    15: 'Regulatory T cell',
    16: 'Helper CD4+ T cell',
    17: 'Helper CD4+ T cell',
    18: 'Helper CD4+ T cell',
}

#add column to adata_nl
adata_nl.obs['louvain'] = adata_nl.obs['louvain'].astype(int)
adata_nl.obs['cell_group'] = adata_nl.obs['louvain'].map(cluster_to_cell_nl)


#visualize most expressed gene in that cluster
#plot top most expressed genes
top_genes = []
for i in transposed_lists:
  top_genes.extend([i[0]])
#print(top_genes)

#Assign the most expressed gene based on max expression in each cluster to each cell
adata_nl.obs['louvain'] = adata_nl.obs['louvain'].astype(int)
adata_nl.obs['most_expressed_gene_max'] = adata_nl.obs['louvain'].map(lambda x: top_genes[x])

#Plot UMAP colored by the most expressed gene for each cluster (based on max expression)
sc.pl.umap(adata_nl, color='most_expressed_gene_max', title="Max Gene Expressed Per Cluster- Normal Lung", show=False,return_fig=True)
file_path = os.path.join(folder_path, "Max Gene Expressed Per Cluster- Normal Lung.png")
plt.savefig(file_path)
plt.show(block=False)

##plot for a given TF subtype that corresponds to CD4
#umap for the dominant TF for a CD4 subtype
adata_nl = adata_nl.copy()

# Compute neighbors if not already done
sc.pp.neighbors(adata_nl)

# Compute neighbors if not already done
sc.pp.neighbors(adata_nl)


#Plot plain UMAP
sc.pl.umap(adata_nl, color = 'cell_group', title = 'Low Cluster Resolution Normal Lung UMAP',return_fig=True)
file_path = os.path.join(folder_path, "Low Cluster Resolution Normal Lung UMAP.png")
plt.savefig(file_path)

# Plot UMAP with louvain clusters
grey_to_purple = LinearSegmentedColormap.from_list("grey_purple", ["#D3D3D3", "#6A0DAD"])
sc.pl.umap(adata_nl, color=['Tbx21','Gata3','Rorc','Foxp3','Bcl6', 'Cd3d'], cmap = grey_to_purple, vmin = 0, vmax = 4)

# Violin plots
adata_nl.var_names_make_unique()
adata_nl.raw = adata_nl
sc.pl.violin(adata_nl,
             keys=['Tbx21','Gata3','Rorc','Foxp3','Bcl6', 'Cd3d'],
             groupby="cell_group",
             jitter=0.4,
             scale="width",
             multi_panel=True, show=False)
file_path = os.path.join(folder_path, "Violin Plot Normal.png")
plt.savefig(file_path)
plt.show()
plt.close()
#endregion ################ END NORMAL UMAP###################

#region #### TUMOR CLUSTERING AND UMAP #########  
#TUMOR PCA reduction
adata_tumor = adata_tumor.copy()

#perform PCA
sc.pp.pca(adata_tumor)
sc.pl.pca_variance_ratio(adata_tumor, log=True, show=False)
plt.show(block=False)
sc.tl.pca(adata_tumor, n_comps=50)

# Compute neighbors if not already done
sc.pp.neighbors(adata_tumor, n_neighbors=24)

# Run louvain clustering with higher resolution for more clustering
sc.tl.louvain(adata_tumor, resolution=1.5)

# Plot UMAP with louvain clusters
sc.tl.umap(adata_tumor, n_components=2)
sc.pl.umap(adata_tumor, color='louvain', show=False)
plt.show(block=False)

louvain_tumor = adata_tumor.obs['louvain'].values

#TUMOR PCA reduction

# Plot UMAP with louvain clusters

#Within each cluster, find the gene that is most expressed.
#used to identify marker genes in each cluster
sc.tl.rank_genes_groups(adata_tumor, groupby='louvain', method='wilcoxon')
# View the top 20 marker genes for each cluster
sc.pl.rank_genes_groups(adata_tumor, n_genes=20, sharey=False,show=False)
plt.show(block=False)
top_genes_indices_tumor = adata_tumor.uns['rank_genes_groups']['names']

# Using zip to transpose the lists (create 15 lists with 20 items each)
transposed_lists_tumor = list(zip(*top_genes_indices_tumor))

# Convert the tuple output of zip into lists
transposed_lists_tumor = [list(t) for t in transposed_lists_tumor]

# Output the result
#for i, transposed_list in enumerate(transposed_lists_tumor):
#    print(f"List {i+1}: {transposed_list}")


   #create a dictionary with the cluster to cell type
cluster_to_cell_tumor = {
    0: 'Regulatory T cell',
    1: 'Helper CD4+ T cell',
    2:'T Follicular Helper cell',
    3: 'Helper CD4+ T cell',
    4: 'Naive CD4+ T cell',
    5: 'Helper CD4+ T cell',
    6: 'Helper CD4+ T cell',
    7: 'Th1',
    8: 'Naive CD4+ T cell',
    9: 'Helper CD4+ T cell',
    10: 'Regulatory T cell',
    11:'Regulatory T cell'
}

#add column to adata_nl
adata_tumor.obs['louvain'] = adata_tumor.obs['louvain'].astype(int)
adata_tumor.obs['cell_type'] = adata_tumor.obs['louvain'].map(cluster_to_cell_tumor)


#create a dictionary with the cluster to cell group (3 total)
cluster_to_cell_tumor = {
    0: 'Regulatory T cell',
    1: 'Helper CD4+ T cell',
    2:'Helper CD4+ T cell',
    3: 'Helper CD4+ T cell',
    4: 'Naive CD4+ T cell',
    5: 'Helper CD4+ T cell',
    6: 'Helper CD4+ T cell',
    7: 'Helper CD4+ T cell',
    8: 'Naive CD4+ T cell',
    9: 'Helper CD4+ T cell',
    10: 'Regulatory T cell',
    11:'Regulatory T cell'
}

#add column to adata_nl
adata_tumor.obs['louvain'] = adata_tumor.obs['louvain'].astype(int)
adata_tumor.obs['cell_group'] = adata_tumor.obs['louvain'].map(cluster_to_cell_tumor)

# Plot UMAP with louvain clusters
sc.pl.umap(adata_tumor, color='cell_group', show=True,title='Tumor UMAP Cell Group for Analysis',return_fig=True)
file_path = os.path.join(folder_path, "Tumor UMAP Cell Group for Analysis.png")
plt.savefig(file_path)
plt.show(block=False)

#plot top most expressed genes
top_genes = []
for i in transposed_lists:
  top_genes.extend([i[0]])
#print(top_genes)

#Assign the most expressed gene based on max expression in each cluster to each cell
adata_tumor.obs['louvain'] = adata_tumor.obs['louvain'].astype(int)
adata_tumor.obs['most_expressed_gene_max'] = adata_tumor.obs['louvain'].map(lambda x: top_genes[x])

#Plot UMAP colored by the most expressed gene for each cluster (based on max expression)
sc.pl.umap(adata_tumor, color='most_expressed_gene_max', title="Max Gene Expressed Per Cluster- Tumor Lung")

adata_tumor = adata_tumor.copy()

# Compute neighbors if not already done
sc.pp.neighbors(adata_tumor)


# Plot UMAP with louvain clusters
grey_to_purple = LinearSegmentedColormap.from_list("grey_purple", ["#D3D3D3", "#6A0DAD"])
sc.pl.umap(adata_tumor, color=['Tbx21','Gata3','Rorc','Foxp3','Bcl6', 'Cd3d'], cmap = grey_to_purple, vmin = 0, vmax = 4)


# Violin plots
adata_tumor.var_names_make_unique()
adata_tumor.raw = adata_tumor
sc.pl.violin(adata_tumor,
             keys=['Tbx21','Gata3','Rorc','Foxp3','Bcl6', 'Cd3d'],
             groupby="cell_group",
             jitter=0.4,
             scale="width",
             multi_panel=True,show=False)
file_path = os.path.join(folder_path, "Violin_tumor.png")
plt.savefig(file_path)
plt.show(block=False)
plt.close()

#endregion ################ END TUMOR UMAP###################

def AnnData_UMAP():
    normal = adata_nl
    tumor = adata_tumor
    return normal, tumor
