from anndata import AnnData
import pandas as pd
import numpy as np
from math import inf
import scanpy as sc
import scipy
from glm.glm import GLM
from glm.families import QuasiPoisson


def _get_model_pars(anndata, bins, latent,):
    latent = np.array(np.log1p(anndata.X.mean(1))).flatten()
    models = [GLM(family=QuasiPoisson()) for _ in range(0, anndata.shape[1])]
    for n in range(0, anndata.shape[1]):
        models[n].fit(pd.DataFrame({
            "y": anndata.X[:, n].toarray().flatten(),
            "log_umi": latent
        }), formula="y~log_umi")
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
    return models


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

def pytransform(anndata, training_cell=inf, min_cells=5, num_genes=2000):
    pass
    # #Initial filtering and sampling
    # sampling_filtering = _step1(anndata,min_cells,num_genes)
    # #Fit model parameters
    # model_pars = _get_model_pars(sampling_filtering[0], bins=1)
    # #Normalize the anndata object
    # normalized_anndata = normalize_pearson_residuals(sampling_filtering,model_pars)
    # return normalized_anndata