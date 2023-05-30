import logging

from anndata import AnnData
import pandas as pd
import numpy as np
from math import inf
import scanpy as sc
import scipy
from glm.glm import GLM
from glm.families import QuasiPoisson
from KDEpy import FFTKDE

_EPS = np.finfo(float).eps

def _outliers(y, x, threshold):
    kde = FFTKDE(kernel="gaussian", bw="ISJ").fit(x)
    width = (max(x)-min(x))*kde.bw / 2
    eps = _EPS * 10
    breaks1 = np.arange(min(x),max(x)+ width,width)
    breaks2 = np.arange(min(x) - eps - width/2,max(x)+width,width)
    scores1 = _robust_scale(y, x, breaks1)
    scores2= _robust_scale(y, x, breaks2)
    return np.abs(np.vstack((scores1, scores2))).min(0)>threshold

def _robust_scale(y, x, breaks): ## Possibly change later, implementation confusing
    bins = np.digitize(x, breaks)
    unique_bins = np.unique(bins)
    res = np.zeros(bins.size)
    for i in range(0, unique_bins.size):
        items = y[bins==unique_bins[i]]
        res[bins==unique_bins[i]] = (items-np.median(items))/(scipy.stats.median_abs_deviation(items)+_EPS)
    return res

def _pearson_residual(y, mu, theta, min_var = -inf):
    model_var = mu + (mu**2 / theta)
    return ( (y - mu) / (np.sqrt(model_var)))

def _regularize(anndata, model_pars, bw_adjust=3):
    anndata.var["Poisson"] = np.where((anndata.var["amean"] < 1e-3), True, False)
    model_pars.var["Poisson"] = np.where((model_pars.var["overdisp_fact"] <= 0)
                                         | (model_pars.var["theta"] == inf), True, False)
    poisson = pd.concat([anndata.var["Poisson"], model_pars.var["Poisson"]], axis=0, join="outer")
    poisson = poisson[~poisson.index.duplicated(keep='first')]
    anndata.var["Poisson"] = poisson
    model_pars.var["dispersion_par"] = np.log10(1+(10**model_pars.var["log_gmeans"]/model_pars.var["theta"]))
    log_gmeans = model_pars.var["log_gmeans"].to_numpy()
    dispersion_outliers = _outliers(model_pars.var["dispersion_par"].to_numpy(), log_gmeans, 10)
    intercept_outliers = _outliers(model_pars.var["intercept"].to_numpy(), log_gmeans, 10)
    coefficient_outliers = _outliers(model_pars.var["coef"].to_numpy(), log_gmeans, 10)
    all_outliers =dispersion_outliers | coefficient_outliers | intercept_outliers \
                  | np.invert((model_pars.var["theta"] == inf).to_numpy())
    model_pars.var["outliers"] = all_outliers
    model_pars = model_pars[:, model_pars.var["outliers"] == False]
    overdispersed_models = model_pars[model_pars.var["Poisson"] == False]
    kde = FFTKDE(kernel="gaussian", bw="ISJ").fit(overdispersed_models.var["log_gmean"])
    bw = kde.bw*bw_adjust
    x_points = np.maximum(anndata.var["log_gmeans"].to_numpy(), min(overdispersed_models.var["log_gmeans"]))
    x_points = np.minimum(x_points, max(overdispersed_models.var["log_gmeans"]))
    o = np.sort(x_points)
    fit_mtx = pd.DataFrame(data=x_points, columns=["x_points"])





    




def _get_model_pars(anndata):
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
    diff_theta = np.array(predicted_theta/actual_theta).flatten()
    for n in range(0, len(diff_theta)):
        if diff_theta[n] < 1e-3:
            models[n].is_overdispersed = True
    theta = np.array([x.dispersion_ for x in models])
    coefficient = np.array([x.coef_[1] for x in models])
    intercept = np.array([x.coef_[0] for x in models])
    anndata.var["step1_models"] = models
    anndata.var["coef"] = coefficient
    anndata.var["intercept"] = intercept
    anndata.var["theta"] = theta
    return anndata


def _step1(anndata, min_cells,  num_cells, num_genes=inf):
    sc.pp.filter_genes(anndata, min_cells=min_cells)
    down_sample = sc.pp.subsample(anndata, n_obs=min(num_cells, anndata.shape[0]), copy=True)
    sc.pp.filter_genes(down_sample, min_cells=1)
    for n in [anndata, down_sample]:
        print(n.shape)
        n.var["amean"] = np.asarray(n.X.mean(0)).flatten()
        n.var["gmean"] = np.asarray(np.expm1(n.X.log1p().mean(0))).flatten()
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
    pass
    # #Initial filtering and sampling
    #sampling_filtering = _step1(anndata, min_cells, min(num_genes, anndata.n_vars), num_cells)
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

