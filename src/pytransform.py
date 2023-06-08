import logging
import os
import xml.dom

from statsmodels.nonparametric.kernel_regression import KernelReg
from anndata import AnnData
import pandas as pd
import numpy as np
from math import inf
import scanpy as sc
import scipy
from glm.glm import GLM
from glm.families import QuasiPoisson
from KDEpy import FFTKDE
import multiprocessing

_EPS = np.finfo(float).eps


def _calc_model(anndata, latent):
    """
    Calculates the poisson regression model per gene. Implemented here for multiprocessing
    :param anndata: anndata object
    :param latent: latent variable to use. The latent variable is per cell. This is log(umi) by default/
    :return: the overdispersion, slope, and intercept for the regression model
    """
    models = [GLM(family=QuasiPoisson()) for _ in range(0, anndata.shape[1])]
    for n in range(0, anndata.shape[1]):
        logging.critical(n)
        models[n].fit(pd.DataFrame(dict(y=anndata.X[:, n].toarray().flatten(), log_umi=latent.flatten())), formula="y~log_umi")
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
    theta = [float(x.dispersion_) for x in models]
    coefficient = [float(x.coef_[1]) for x in models]
    intercept = [float(x.coef_[0]) for x in models]
    return theta, coefficient, intercept

def _outliers(y, x, threshold):
    """
    Calculates if any of the parameters are outliets
    :param y:  the variable to calculate the outlier against
    :param x: the target to fit against
    :param threshold: the outlier threshold
    :return: THe parameters that are outliers
    """
    kde = FFTKDE(kernel="gaussian", bw="ISJ").fit(x)
    width = (max(x)-min(x))*kde.bw / 2
    eps = _EPS * 10
    breaks1 = np.arange(min(x),max(x)+ width,width)
    breaks2 = np.arange(min(x) - eps - width/2,max(x)+width,width)
    scores1 = _robust_scale(y, x, breaks1)
    scores2= _robust_scale(y, x, breaks2)
    return np.abs(np.vstack((scores1, scores2))).min(0)>threshold

def _robust_scale(y, x, breaks): ## Possibly change later, implementation confusing
    """
    :param y: Metric to calculate scale for
    :param x: Metric to scale against
    :param breaks: number of bins
    :return: outlier scale
    """
    bins = np.digitize(x, breaks)
    unique_bins = np.unique(bins)
    res = np.zeros(bins.size)
    for i in range(0, unique_bins.size):
        items = y[bins==unique_bins[i]]
        res[bins==unique_bins[i]] = (items-np.median(items))/(scipy.stats.median_abs_deviation(items)+_EPS)
    return res


def _regularize(anndata, model_pars, bw_adjust=3):
    """
    Regularizes all genes in the anndata object based on the previously fitted outlier models.
    :param anndata: AN anndata object
    :param model_pars: theta, slope, intercept
    :param bw_adjust: The factor to scale the kernel density estimator by
    :return: The anndata matrix with model pars estimated for all genes.
    """
    anndata.var["Poisson"] = np.where((anndata.var["amean"] < 1e-3), True, False)
    model_pars.var["Poisson"] = np.where((model_pars.var["overdisp_fact"] <= 0)
                                         | (model_pars.var["theta"] == inf), True, False)
    poisson = pd.concat([anndata.var["Poisson"], model_pars.var["Poisson"]], axis=0, join="inner")
    poisson = poisson[~poisson.index.duplicated(keep='first')]
    anndata.var["Poisson"] = poisson
    model_pars.var["dispersion_par"] = np.log10(1+10**model_pars.var["log_gmeans"]/(model_pars.var["theta"]))
    log_gmeans = model_pars.var["log_gmeans"].to_numpy()
    dispersion_outliers = _outliers(model_pars.var["dispersion_par"].to_numpy(), log_gmeans, 10)
    intercept_outliers = _outliers(model_pars.var["intercept"].to_numpy(), log_gmeans, 10)
    coefficient_outliers = _outliers(model_pars.var["coef"].to_numpy(), log_gmeans, 10)
    all_outliers =dispersion_outliers | coefficient_outliers | intercept_outliers \
                  | (model_pars.var["theta"] == inf).to_numpy()
    model_pars.var["outliers"] = all_outliers
    poisson = pd.concat([anndata.var["Poisson"], model_pars.var["outliers"]], axis=0, join="inner")
    poisson = poisson[~poisson.index.duplicated(keep='first')]
    anndata.var["Poisson"] = poisson
    model_pars = model_pars[:, model_pars.var["outliers"] == False]
    anndata = anndata[:, anndata.var["Poisson"]==False]

    overdispersed_models = model_pars[:, model_pars.var["Poisson"] == False]
    kde = FFTKDE(kernel="gaussian", bw="ISJ").fit(overdispersed_models.var["log_gmeans"].to_numpy())
    bw = kde.bw*bw_adjust
    x_points = np.maximum(anndata.var["log_gmeans"].to_numpy(), min(overdispersed_models.var["log_gmeans"]))
    x_points = np.minimum(x_points, max(overdispersed_models.var["log_gmeans"]))
    o = np.sort(x_points)
    fit_mtx = pd.DataFrame(index =anndata.var.index, data=x_points, columns=["x_points"])
    for n in ["dispersion_par", "intercept", "coef"]:
        ks = KernelReg(overdispersed_models.var[n].to_numpy(), overdispersed_models.var["log_gmeans"].to_numpy(),var_type="c", bw=[bw])
        fit_mtx[n] = ks.fit(o)[0]
    fit_mtx["theta"] = (10**anndata.var["log_gmeans"])/(10**fit_mtx["dispersion_par"].to_numpy()-1)
    sum_mean = anndata.X.mean(1).sum()
    for n in anndata[:, anndata.var["Poisson"]].var_names:
        fit_mtx.at[n, "dispersion_par"] = 0
        fit_mtx.at[n, "theta"] = inf
        fit_mtx.at[n, "coef"] = 1
        fit_mtx.at[n, "intercept"] = anndata.var.loc[n, "amean"] - sum_mean

    return fit_mtx




def _get_residuals(anndata, model_pars):
    """
    :param anndata: Anndata object
    :param model_pars: slope, coefficient, overdispersion pars for the anndata
    :return: Pearson residuals for the anndata object
    """
    median = np.apply_along_axis(lambda v: np.median(v[np.nonzero(v)]), 0, anndata.X.toarray())
    min_var = (median/5)**2
    latent = np.array(np.log10(anndata.X.sum(1))).flatten()
    X = anndata.X.copy()
    params = model_pars[["intercept", "coef"]]
    d = X.data.copy()
    x,y = X.nonzero()
    mu = np.exp(params.values[:,0][y]+params.values[:,1][y]*latent[x])
    var = mu + (mu**2/model_pars["theta"].values.flatten()[y])
    for ind, n in enumerate(var):
        if n < min_var[y[ind]]:
            var[ind] = min_var[y[ind]]
    X.data[:]=d-(mu/var**0.5)
    X.data[X.data<0]=0
    X.eliminate_zeros()
    clip = np.sqrt(X.shape[0]/30)
    X.data[X.data>clip]=clip
    return X



    




def _get_model_pars(anndata, num_workers):
    """

    :param anndata: anndata object
    :param num_workers: number of cores to uses
    :return: anndata with dispersion pars estimated for each gene in the object
    """
    latent = np.array(np.log10(anndata.X.sum(1))).flatten()
    split = np.array_split(anndata.X.toarray(), num_workers, axis=1)
    pars_list = []
    z = 0
    for n in split:
        pars_list.append((anndata[:, z:z + n.shape[1]], latent))
    with multiprocessing.Pool(num_workers) as p:
        result = p.starmap(_calc_model, pars_list)
    coefficient = []
    intercept = []
    theta = []
    for n in result:
        theta += n[0]
        coefficient += n[1]
        intercept += n[2]
    anndata.var["coef"] = coefficient
    anndata.var["intercept"] = intercept
    anndata.var["theta"] = theta

    return anndata


def _step1(anndata, min_cells,  num_cells, num_genes=2000):
    """
    :param anndata: anndata object
    :param min_cells: minimum number of cells for a gene to be expressed in to be retained
    :param num_cells: number of cells to use for the initial regression step
    :param num_genes: The number of genes to estimate parameters for in the initial regression step
    :return: an anndata object downsampled from the original anndata to num_cells cells and num_genes genes.
    """
    sc.pp.filter_genes(anndata, min_cells=min_cells)
    down_sample = sc.pp.subsample(anndata, n_obs=min(num_cells, anndata.shape[0]), copy=True)
    sc.pp.filter_genes(down_sample, min_cells=1)
    for n in [anndata, down_sample]:
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
        kde= scipy.stats.gaussian_kde(np.asarray(gmeans).flatten())
        probs=kde.pdf(np.asarray(gmeans).flatten())
        down_sample_filt=down_sample_filt[:, down_sample_filt.var.sample(n=num_genes, weights=probs).index]
    return down_sample_filt

def pytransform(anndata, min_cells=5, num_genes=2000, num_cells=5000, workers=os.cpu_count()-1, inplace=True,
                verbose=False, highly_variable=3000):
    """
    Performs pearson residual estimation on the anndata object. Please note that the initial regression estimation is
    very compute intensive, therefore it is advised to not select too many genes for the num_genes parameter
    :param anndata: an anndata object **NOTE: BATCH EFFECTS ARE CURRENTLY UNSUPPORTED. PLEASE SPLIT THE ANNDATA BEFORE
    YOU APPLY PYTRANSFORM**
    :param min_cells: The min number of cells a gene should be expressed in to not be filtered out. Default is 5
    :param num_genes: Number of genes to use in the initial regression step. Default is 2000
    :param num_cells: Number of cells to use in hte initial regression step. Default is 5000
    :param workers: Number of CPU cores to use. Default is available_cores-1
    :param inplace: Whether to modify the anndata object or return a copy.
    :param verbose: Currently, unimplemented. In the future, will set the level of messages to output
    :param highly_variable: The number of highly variable genes to retain in the final object. Set to 0 to retain all
    genes. Note: scanpy's default highly variable gene selection is not designed to take in pearson residuals, and will
    output strange results if used with these. If you intend to select highly variable genes, please use this function
    instead.
    :return: An anndata object if inplace=False, else nothing.
    """
    if verbose:
        logging.basicConfig(level=logging.INFO)
    sub_samp = _step1(anndata, min_cells=min_cells, num_genes=num_genes, num_cells=num_cells)
    models = _get_model_pars(sub_samp, workers)
    params = _regularize(anndata, models)
    anndata = anndata[:, anndata.var["Poisson"]==False]
    residuals = _get_residuals(anndata, params)
    if inplace:
        anndata.X = residuals
        if highly_variable > 0:
            x_sq = anndata.X.copy()
            x_mean = anndata.X.mean(0)
            x_sq.data **=2
            pearson_var = np.asarray(x_sq.mean(0)) - np.asarray(np.square(x_mean))
            anndata.var["Pearson_variance"] = pearson_var.flatten()
            variable_indeces = anndata.var["Pearson_variance"].sort_values(ascending=False).iloc[0:highly_variable].index
            anndata = anndata[:, variable_indeces]
            anndata.X = scipy.sparse.csr_matrix(anndata.X)
    else:
        return_val = anndata.copy()
        return_val.X = residuals
        if highly_variable > 0:
            x_sq = return_val.X.copy()
            x_mean = return_val.X.mean(0)
            x_sq.data **= 2
            pearson_var = np.asarray(x_sq.mean(0)) - np.asarray(np.square(x_mean))
            return_val.var["Pearson_variance"] = pearson_var.flatten()
            variable_indeces = return_val.var["Pearson_variance"].sort_values(ascending=False).iloc[0:highly_variable].index
            return_val = return_val[:, variable_indeces]
            return_val.X = scipy.sparse.csr_matrix(return_val.X)
        return return_val



