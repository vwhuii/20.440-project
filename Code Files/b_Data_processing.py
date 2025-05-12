import scipy as sp
from matplotlib import pyplot as plt
import scanpy as sc
import anndata as ad

from a_Import_Data import AnnData_Object
adata_nl, adata_tumor = AnnData_Object()

#region #### PROCESSING DATA #########   
#normal data
adata_nl.raw = adata_nl.copy() #make a copy of the raw data
sc.pp.filter_cells(adata_nl, min_genes=200)  # Remove low-quality cells
sc.pp.filter_genes(adata_nl, min_cells=3)  # Remove rarely expressed genes
sc.pp.normalize_total(adata_nl, target_sum=1e4)  # Normalize counts per cell
sc.pp.log1p(adata_nl)  # Apply log transformation
sc.pp.highly_variable_genes(adata_nl, flavor="seurat", n_top_genes=2000) #top 2000 highly var genes

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


#add qc metrics to adata
sc.pp.calculate_qc_metrics(adata_tumor, inplace=True)

#label all mitochondrial genes as 'true' for mitochondrial
adata_tumor.var["mt"] = adata_tumor.var.index.str.startswith("mt-")

#endregion ################ END PROCESS DATA###################

def AnnData_Processed():
    normal = adata_nl
    tumor = adata_tumor
    return normal, tumor
