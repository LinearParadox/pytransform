Pytransform is a python port of SCTransform, a variance stabilization method for single cell RNA sequencing data.  
Instead of taking the log of the data and scaling it pearson residuals are generated using a generalized linear model.  
For more information, please see the following:  
[Normalization and variance stabilization of single-cell RNA-seq data using regularized negative binomial regression](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1874-1)  
In addition, we have implemented an automated filtering function which scores cells based on absolute deviation from  
the median.  

#Installation:  
To install, please use the following:  

pip install git+https://github.com/LinearParadox/pytransform.git
  
  
## Pytransform  

from pytransform import pytransform

useage: pytransform(anndata, min_cells=5, num_genes=2000, num_cells=5000, workers=os.cpu_count()-1, inplace=True,
                verbose=False)  

anndata -- The anndata object to apply pearson residuals to  
min_cells -- minimum number of cells a gene must be expressed in to be retained  
num_genes -- number of genes to use in the initial fitting step  
num_cells -- number of cells to use in the initial fitting step  
workers -- number of cpu cores to use  
inplace -- default is true. Whether to return a new anndata object or to modify the provided object.    
**Save the counts  
in adata.raw, this will overwrite them**  
verbose -- Default false. Currently does nothing, in the future will print information if True.
---------------------------------------------------------------------------------------------------------  
  
  


### Function Signature: add_MAD(anndata, vars)  


Parameters

anndata: An AnnData object

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
Data should be loaded in using the standard scanpy workflow. Observations should be added, such as mitochondrial percent.
Execution Steps

1. Import the necessary libraries and the function:
2. Load your data into an AnnData object.
3. Compute percent MT, as well as read count distributions
4.Define the variables you want to use for MAD based filtering:
5.Run the add_MAD function:
6. The modified adata object will now contain new columns in its .obs attribute with the prefix 'MAD_', representing the deviation of each observation from the median for the respective variable.

---------
