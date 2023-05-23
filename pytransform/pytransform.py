import anndata as ad
import pandas as pd
import numpy as np
from math import inf
import scanpy as sc
import scipy


def _get_model_pars(anndata, bins, latent_vars):
    pass
    # Implement poisson fitting
    #input anndata object, already filtered using the step one function
    # should ideally include latent vars too.
    # bins is for multicore support, where you bin the count matrix into m cells by (n/bins) genes


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
    down_sample_filt=down_sample_filt[:, down_sample_filt.var.sample(n=num_genes, weights=probs).index]
    log_gmeans = np.log10(np.expm1(down_sample_filt.X.log1p().mean(0)))
    return down_sample_filt, log_gmeans







def pytransform(anndata, training_cell=inf, min_cells=5, num_genes=2000):
    pass