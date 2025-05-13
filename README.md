# 20.440 Project: The Role of CD4+ Subsets in Lung Adenocarcinoma Progression Using Gene Expression

## Overview

This repo aims to analyze the dataset (GSE152022) to further analyze CD4+ T cell phenotypes in lung adenocarcinoma. First, CD4+ T cell subsets are identified and labeled via PCA, louvain, and UMAP display. Then, cells are grouped into three CD4+ subsets: naive, helper, and regulatory. Marker genes are evaluated. Differential Gene Expression Analysis is conducted on cells from normal or tumor lungs via log fold change. Then, genes are compared between cell groups from the tumor condition where their logFC values are normalized and displayed via gene expression score (Z-score) to be compared across groups. Lastly, gene set enrichment analysis (GSEA) is completed for a subset of genes specific to CD4+ activation, CD4+ exhaustion, and CD4+ secretion of M2-polarizing factors. 

This repo contains all necessary code to run the analysis for the project, "Uncovering the Role of CD4+ Subsets in Lung Adenocarcinoma Progression Using Gene Expression".

All code blocks were written by Logan and Veronica with exception to publically available libraries that were imported in the beginning. 
Assignment of CD4+ T cell subset using top most expressed genes were performed with help of Enrichr (https://maayanlab.cloud/Enrichr/)


## Data

The authors of this dataset performed scRNA-seq and TCR sequencing analysis of enriched CD4+ T cells from pooled lung tissues of normal mice (n=10) and murine lung adenocarcinoma tumor bearing mice (n=7) (GSE152022). The data can be processed to better understand various CD4+ subtypes in the presence of lung adenocarcinoma in murine models. 
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152022 

## Folder Structure

20.440 Project Data: Contains all the data files extracted from the dataset linked above needed for analysis. It includes the barcodes, features, and gene expression matrix for each tumor and normal cell. Raw data can also be downloaded from the GEO link directly.

Figures: Contains the saved figures that are generated from this report. 

Code: Contains the codes file to generate the analysis in the order of a, b, c, d, e. Relevant variables are carried over via functions to ensure that they are connected and updated across sub .py files 

a_Import_Data: imports the datafiles from a local folder into the variables adata_nl and adata_tumor for both normal and tumor samples. 

b_Data_process: processes the data through filtering and quality control metrics

c_UMAP_MarkerGene_analysis: performs dimensionality reduction, clustering via leiden, and UMAP. Marker Gene Analysis is also embedded into this file

d_DiffGeneExps_analysis: performs all differential gene expression analysis between normal v tumor, as well as between three CD4+ cell groups.

e_GSEA_analysis: conducts GSEA analysis


## Installation

You will have to install the libraries needed to complete the data analysis. First, if not downlaoded already, download the newest version of Python 3.13.3. Then, you can download the necessary libraries required to run this script by using pip install library-name. 
Libraries that are used and the versions used include: scipy (1.15.2), matplotlib (3.10.1), seaborn (0.13.2), pandas (2.2.3), scanpy (1.11.1), anndata (0.11.4), scipy (1.5.2), celltypist (1.6.3), umap-learn (0.5.7), scikit-learn (1.5.2). <br>

To initialize, follow these steps prior to running code in VS Code or any similar Python processor: <br>
1. Make sure that pip is installed and that Python is integrated in the PATH. The instructions in steps 2-4 are done in the terminal.<br>

2. Install conda software (https://www.anaconda.com/docs/getting-started/anaconda/install#windows-installation) <br>
3. Clone the repository <br>
git clone https://github.mit.edu/lbeatty/20.440-pset6-done <br>
cd 20.440-finalproject <br>
*Extract files from a zipped folder if needed <br>
4. Set up the conda environment (using the environment.yml file) <br>
conda env create -f environment.yml <br>
conda activate 440env #440env is the environment name used in the original setup <br>
5. Install packages noted in requirements.txt (scanpy, scipy, numpy, matplotlib, pandas) <br>
pip install -r requirements.txt <br>
6. Run .py code files <br>

## References:
Son J, Cho JW, Park HJ, Moon J et al. Tumor-Infiltrating Regulatory T-cell Accumulation in the Tumor Microenvironment Is Mediated by IL33/ST2 Signaling. Cancer Immunol Res 2020 Nov;8(11):1393-1406. PMID: 32878747

Xie Z, Bailey A, Kuleshov MV, Clarke DJB., Evangelista JE, Jenkins SL, Lachmann A, Wojciechowicz ML, Kropiwnicki E, Jagodnik KM, Jeon M, & Ma’ayan A. Gene set knowledge discovery with Enrichr. Current Protocols, 1, e90. 2021. doi: 10.1002/cpz1.90
