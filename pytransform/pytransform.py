import anndata as ad
import pandas as pd
import numpy as np
from math import inf
import scanpy as sc

def _step1(anndata, sample, min_cells):
    sc.pp.filter_genes(anndata, min_cells=min_cells)
    down_sample = sc.pp.subsample(anndata, n_obs=min(sample, anndata.shape[0]), copy=True)
    gmeans = np.log10(down_sample.X.mean(0))
    amean = down_sample.X.mean(0)
    x_sq = down_sample.X.copy()
    x_sq.data **= 2
    genes_var = x_sq.mean(0)-np.square(amean)
    overdisp_fact = genes_var-amean
    is_overdisp = overdisp_fact > 0
    down_sample_filt = down_sample[:, list(np.asarray(is_overdisp).flatten())]
    gmeans = gmeans[is_overdisp]



def pytransform(anndata, training_cell=inf, min_cells=5):
    pass