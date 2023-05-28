import logging

from anndata import AnnData
import pandas as pd
import numpy as np
from math import inf
import scanpy as sc
import scipy
from glm.glm import GLM
from glm.families import QuasiPoisson


def _regularize(anndata, model_pars):
    pass

    




def _get_model_pars(anndata, latent,):
    latent = np.array(np.log1p(anndata.X.mean(1))).flatten()
    models = [GLM(family=QuasiPoisson()) for _ in range(0, anndata.shape[1])]
    for n in range(0, anndata.shape[1]):
        models[n].fit(pd.DataFrame(dict(y=anndata.X[:, n].toarray().flatten(), log_umi=latent)), formula="y~log_umi")
    means = anndata.X.mean(0)
    x_sq = anndata.X.copy()
    x_sq.data **= 2
    genes_var = x_sq.mean(0) - np.square(means)
    predicted_theta = np.square(means)/(genes_var-means)
    actual_theta = np.array([x.dispersion_ for x in models])
    diff_theta = np.array(predicted_theta/actual_theta)
    for n in range(0, len(diff_theta)):
        if diff_theta[n] < 1e-3:
            models[n].is_overdispersed = True
    anndata.vars["step1_models"] = models
    return anndata


def _step1(anndata, min_cells, num_genes, num_cells):
    sc.pp.filter_genes(anndata, min_cells=min_cells)
    down_sample = sc.pp.subsample(anndata, n_obs=min(num_cells, anndata.shape[0]), copy=True)
    for n in [anndata, down_sample]:
        n.var["amean"] = np.asarray(n.X.mean(0)).flatten()
        n.var["gmean"] = np.asarray(np.expm1(down_sample.X.log1p().mean(0))).flatten()
        n.var["log_gmeans"] = np.asarray(np.log10(np.expm1(n.X.log1p().mean(0)))).flatten()
        x_sq = n.X.copy()
        x_sq.data **= 2
        genes_var = np.asarray(x_sq.mean(0)) - np.asarray(np.square(n.var["amean"]))
        n.var["var"] = genes_var.flatten()
        n.var["overdisp_fact"] = n.var["var"] - n.var["amean"]
    down_sample_filt = down_sample[:, down_sample.var["overdisp_fact"]>0]
    gmeans = np.asarray(down_sample_filt.var["gmean"]).flatten()
    if num_genes < down_sample_filt.n_vars: # If you want to downsample genes
        logging.info("""Subsampling to {} genes. Note this will limit the number of genes present in the final product.
                     This is not recommended unless needed to speed up computation.""".format(num_genes))
        kde= scipy.stats.gaussian_kde(np.asarray(gmeans).flatten())
        probs=kde.pdf(np.asarray(gmeans).flatten())
        down_sample_filt=down_sample_filt[:, down_sample_filt.var.sample(n=num_genes, weights=probs).index]
    return down_sample_filt

def pytransform(anndata, min_cells=5, num_genes=inf, num_cells=5000, verbose=False ):
    if verbose:
        logging.basicConfig(level=logging.INFO)
    else:

    pass
    # #Initial filtering and sampling
    sampling_filtering = _step1(anndata, min_cells, min(num_genes, anndata.n_vars), num_cells)
    # #Fit model parameters
    # model_pars = _get_model_pars(sampling_filtering[0], bins=1)
    # #Normalize the anndata object
    # normalized_anndata = normalize_pearson_residuals(sampling_filtering,model_pars)
    # return normalized_anndata


###

# def normalize_pearson_residuals(
#     adata: AnnData,
#     model_pars: np.ndarray,
#     *,
#     clip: Optional[float] = None,
#     check_values: bool = True,
#     layer: Optional[str] = None,
#     inplace: bool = True,
#     copy: bool = False,
# ) -> Optional[Dict[str, np.ndarray]]:
#     if copy:
#         if not inplace:
#             raise ValueError("`copy=True` cannot be used with `inplace=False`.")
#         adata = adata.copy()
#
#
#     X = _get_obs_rep(adata, layer=layer)
#     computed_on = layer if layer else 'adata.X'
#
#     # Use the provided model parameters (e.g., theta estimates) for each gene
#     residuals = _pearson_residuals(X, model_pars, clip, check_values, copy=~inplace)
#     settings_dict = dict(theta=model_pars, clip=clip, computed_on=computed_on)
#
#     if inplace:
#         _set_obs_rep(adata, residuals, layer=layer)
#         adata.uns['pearson_residuals_normalization'] = settings_dict
#     else:
#         results_dict = dict(X=residuals, **settings_dict)
#
#     return adata if copy else (results_dict if not inplace else None)
##

