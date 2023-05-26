Our Project

Group 1 aims to improve and automating single cell analysis, including adding additional methods for the qc filtering process as well as normalization aspects. We want to compare our Scanpy implementation additions
against default filtering options we've learned in lab as well as against Seurat, comparing accuracy and precision, what gets filtered out, and if there's improvement in filtering speed/runtime. 



--
Function Signature: add_MAD(anndata, vars)

Parameters

anndata: An AnnData object which is a data structure for annotated data that contains observed data (.obs), variables (.var), and unspecific annotations (.uns).

vars: A list of variables present in anndata.obs for which the MAD based filtering will be performed.

The function goes through the list of variables specified in vars and for each variable, it computes the median,subtracts the median from each observation, resulting 
in the deviation of each observation from the median, and then adds a new column to the anndata.obs DataFrame with the prefix 'MAD_' followed by the variable name, 
which stores the deviation of each observation from the median. The function returns the AnnData object with additional columns in anndata.obs for each variable processed, 
containing the deviation of each observation from the median for the respective variable. This function can be used for preprocessing scRNA-seq data in order to identify and 
filter outliers based on the MAD of specified variables.

Running the add_MAD function

Pre-requisites
Before running the add_MAD function, make sure you have the following installed:
Python 3.6 or higher.
Packages: scanpy and anndata. 

Data Preparation
Prepare your single-cell RNA sequencing data in an AnnData format. Ensure the AnnData object has observed data .obs which should contain the variables that you wish to run the add_MAD function on.

Execution Steps

1. Import the necessary libraries and the function:
2. Load your data into an AnnData object. Here's an example with a hypothetical data file:
3.Define the variables you want to use for MAD based filtering:
4.Run the add_MAD function:
5. The modified adata object will now contain new columns in its .obs attribute with the prefix 'MAD_', representing the deviation of each observation from the median for the respective variable.

---------

Files

filter_cells.py - Contains add_MAD, used for adding filtering obs to a given anndata instance for usage in filtering

filter_cells_checks.py - Contains test information for add_MAD, and an example about how it could be used for filtering.

pytransform.py - Contains (currently unfinished) \_step1, used to filter and remove genes based on dispersal ratio, and \_get_model_pars, which fits the anndata to a quasi-Poisson model for distribution.
