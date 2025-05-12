# 20.440 Project: Uncovering the Role of CD4+ T Cells in Lung Adenocarcinoma Progression

## Overview

This repo aims to analyze the dataset (GSE152022) to further analyze CD4+ T cell phenotypes in lung adenocarcinoma. First CD4+ T cell subsets are identified and labeled via PCA, leiden, and UMAP display. Then, basic analysis on CD4+ T cell subset frequency and gene expression is analyzed. Gene set enrichment analysis (GSEA) is completed for a subset of genes specific to CD4+ activation, CD4+ exhaustion, and CD4+ secretion of M2-polarizing factors. 

This repo contains all necessary code to run the analysis for the project, "Uncovering the Role of CD4+ T Cells in Lung Adenocarcinoma Progression".

All code blocks were written by Logan and Veronica with exception to publically avaiable libraries that were imported in the beginning. 
Assignment of CD4+ T cell subset using top most expressed genes were performed with help of Enrichr (https://maayanlab.cloud/Enrichr/)

The one figure that is produced and separately saved as a .png for this assignment (pset6) has been isolated to its own .py file and can be found under Code -> assignment_figure.py

## Data

The authors of this dataset performed scRNA-seq and TCR sequencing analysis of enriched CD4+ T cells from pooled lung tissues of normal mice (n=10) and murine lung adenocarcinoma tumor bearing mice (n=7) (GSE152022). The data can be processed to better understand various CD4+ subtypes in the presence of lung adenocarcinoma in murine models. 
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152022 

## Folder Structure

20.440 Project Data: Contains all the data files extracted from the dataset linked above needed for analysis. It includes the barcodes, features, and gene expression matrix for each tumor and normal cell. Raw data can also be downloaded from the GEO link directly.

Figures: Contains the saved figures that are generated from this report. 

Code: Contains the codes file to generate the analysis in the order of a, b, c, d, e. Relevant variables are carried over via functions to ensure that they are connected and updated across sub .py files 

a_Import_Data: imports the datafiles from a local folder into the variables adata_nl and adata_tumor for both normal and tumor samples. 

b_Data_process: processes the data through filtering and quality control metrics

c_UMAP_analysis: performs dimensionality reduction, clustering via leiden, and UMAP.

d_TF_analysis: conducts analysis looking at the main TF for each CD4+ subset

e_Gene_analysis: conducts analysis with the top 20 most variable genes. 

assignment_figure: contains the one figure that is produced and saved in the "Figures" folder for the purposes of Pset6. It was isolated from the rest of the code found in c_UMAP_analysis for ease of running purposes. 


440projectcode.py in the main folder is the entire .py file unseparated for tracking purposes. The code is then broken down into sub sections as seen in the "Code" folder.

## Installation

You will have to install the libraries needed to complete the data analysis. First, if not downlaoded already, download the newest version of Python 3.13.3. Then, you can download the necessary libraries required to run this script by using pip install library-name. 
Libraries that are used and the versions used include: scipy (1.15.2), matplotlib (3.10.1), seaborn (0.13.2), pandas (2.2.3), scanpy (1.11.1), anndata (0.11.4), scipy (1.5.2), celltypist (1.6.3), umap-learn (0.5.7), scikit-learn (1.5.2), leidenalg (0.10.2).

## References:
Son J, Cho JW, Park HJ, Moon J et al. Tumor-Infiltrating Regulatory T-cell Accumulation in the Tumor Microenvironment Is Mediated by IL33/ST2 Signaling. Cancer Immunol Res 2020 Nov;8(11):1393-1406. PMID: 32878747

Xie Z, Bailey A, Kuleshov MV, Clarke DJB., Evangelista JE, Jenkins SL, Lachmann A, Wojciechowicz ML, Kropiwnicki E, Jagodnik KM, Jeon M, & Ma’ayan A. Gene set knowledge discovery with Enrichr. Current Protocols, 1, e90. 2021. doi: 10.1002/cpz1.90
