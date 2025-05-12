import scipy as sp
from pylab import *
import pandas as pd
import scanpy as sc
import anndata as ad
import os

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

def AnnData_Object():
    normal = adata_nl
    tumor = adata_tumor
    return normal, tumor