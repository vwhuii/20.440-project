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

from b_Data_processing import AnnData_Processed
adata_nl, adata_tumor = AnnData_Processed()


##plot for a given TF subtype that corresponds to CD4
#umap for the dominant TF for a CD4 subtype
adata_nl = adata_nl.copy()

# Compute neighbors if not already done
sc.pp.neighbors(adata_nl)

# Run leiden clustering
sc.tl.leiden(adata_nl, resolution=0.8)

#Plot plain UMAP
sc.pl.umap(adata_nl, color = 'leiden', title = 'Low Cluster Resolution Normal Lung UMAP', show=False)
plt.show(block=False)

# Plot UMAP with leiden clusters
sc.pl.umap(adata_nl, color=['Tbx21','Gata3','Rorc','Foxp3','Bcl6', 'Cd3d'], cmap = 'plasma', show=False)
plt.show(block=False)

# Violin plots
adata_nl.var_names_make_unique()
adata_nl.raw = adata_nl
sc.pl.violin(adata_nl,
             keys=['Tbx21','Gata3','Rorc','Foxp3','Bcl6', 'Cd3d'],
             groupby="leiden",
             jitter=0.4,
             scale="width",
             multi_panel=True, show=False)
plt.show(block=False)
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

# Run leiden clustering with higher resolution for more clustering
sc.tl.leiden(adata_tumor, resolution=1.5)

# Plot UMAP with leiden clusters
sc.tl.umap(adata_tumor, n_components=2)
sc.pl.umap(adata_tumor, color='leiden', show=False)
plt.show(block=False)

leiden_tumor = adata_tumor.obs['leiden'].values

#TUMOR PCA reduction


adata_tumor = adata_tumor.copy()

# Compute neighbors if not already done
sc.pp.neighbors(adata_tumor)

# Run leiden clustering
sc.tl.leiden(adata_tumor, resolution=0.8)

#Plot plain UMAP
sc.pl.umap(adata_tumor, color = 'leiden', title = 'Low Cluster Resolution Tumor Lung UMAP', show=False)
plt.show(block=False)

# Plot UMAP with leiden clusters
sc.pl.umap(adata_tumor, color=['Tbx21','Gata3','Rorc','Foxp3','Bcl6', 'Cd3d'], cmap = 'plasma', show=False)
plt.show(block=False)

# Violin plots
adata_tumor.var_names_make_unique()
adata_tumor.raw = adata_tumor
sc.pl.violin(adata_tumor,
             keys=['Tbx21','Gata3','Rorc','Foxp3','Bcl6', 'Cd3d'],
             groupby="leiden",
             jitter=0.4,
             scale="width",
             multi_panel=True, show=False)
plt.show(block=False)

#endregion ################ END TUMOR UMAP###################

def AnnData_UMAP():
    normal = adata_nl
    tumor = adata_tumor
    return normal, tumor
