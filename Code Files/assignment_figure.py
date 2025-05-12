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

#folder to save figures
folder_path = r"C:\Users\veron\Documents\GitHub\20.440_pset6\Figures"

#region #### NORMAL CLUSTERING AND UMAP #########   
#normal
#perform PCA
sc.pp.pca(adata_nl)
sc.pl.pca_variance_ratio(adata_nl, log=True, show=False)
plt.show(block=False)
sc.tl.pca(adata_nl, n_comps=50)

# Compute neighbors if not already done
sc.pp.neighbors(adata_nl, n_neighbors=15)

# Run leiden clustering with higher resolution for more clustering
sc.tl.leiden(adata_nl, resolution=1.5)

# Plot UMAP with leiden clusters
sc.tl.umap(adata_nl, n_components=2)
sc.pl.umap(adata_nl, color='leiden', show=False)
file_path = os.path.join(folder_path, "normal UMAP with leiden clusters.png")
plt.savefig(file_path)
plt.show(block=False)
