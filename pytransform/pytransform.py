import anndata as ad
import pandas as pd
import numpy as np
from math import inf
import scanpy as sc
import scipy



def _step1(anndata, sample, min_cells, num_genes):
    down_sample = sc.pp.subsample(anndata, n_obs=min(sample, anndata.shape[0]), copy=True)
    sc.pp.filter_genes(down_sample, min_cells=min_cells)
    gmeans = np.log10(down_sample.X.mean(0))
    amean = down_sample.X.mean(0)
    x_sq = down_sample.X.copy()
    x_sq.data **= 2
    genes_var = x_sq.mean(0)-np.square(amean)
    overdisp_fact = genes_var-amean
    is_overdisp = overdisp_fact > 0
    down_sample_filt = down_sample[:, list(np.asarray(is_overdisp).flatten())]
    gmeans = gmeans[is_overdisp]
    kde= scipy.stats.gaussian_kde(np.asarray(gmeans).flatten())
    probs=kde.pdf(np.asarray(gmeans).flatten())
    down_sample_filt=down_sample_filt.var[down_sample_filt.var.sample(2000, weights=probs, axis=0).index]







def pytransform(anndata, training_cell=inf, min_cells=5, num_genes=2000):
    pass