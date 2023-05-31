import numpy as np
import pandas as pd
import anndata as ad
from scipy.sparse import csr_matrix
import scanpy as sc
print(ad.__version__)

counts = csr_matrix(np.random.poisson(1, size=(10, 10)), dtype=np.float32)

adata = ad.AnnData(counts)
adata.obs_names = [f"Cell_{i:d}" for i in range(adata.n_obs)]
adata.var_names = [f"Gene_{i:d}" for i in range(adata.n_vars)]
print(adata.obs_names[:10])

adata.obs["num"] = [0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0]


def add_MAD(anndata, vars):
    for n in vars:
        median = anndata.obs[n].median()
        anndata.obs["MAD_"+n] = anndata.obs[n] - median
    return anndata

v = ["num","rand"]
#adata = add_MAD(adata, v)
#adata = adata[adata.obs["MAD_num"]<3]
adata.obs["rand"] = np.random.poisson(size=10)
adata = add_MAD(adata, v)
print(adata.obs)
adata = adata[(adata.obs["MAD_num"]<3) & (adata.obs["MAD_rand"]>=0.0)]
print(adata.obs)
