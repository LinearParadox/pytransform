import anndata as ad

def add_MAD(anndata, vars):
    for n in vars:
        median = anndata.obs[n].median()
        anndata.obs["MAD_"+n] = anndata.obs[n] - median
    return anndata


#For mitochondrias specifically, we only care about the right end of the tail (positive values)

