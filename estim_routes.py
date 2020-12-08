"""
Estimate frequencies of alternative allele in locations and clusters
"""
import numpy as np
import os
from itertools import product
from pop_hist_opt_wright_fisher import HistOpt, PopData, ClusterCov

import warnings

# warnings.filterwarnings("ignore", category=RuntimeWarning)

def run_freq_bayes(cl_tag, dist_tag, pop_tag):
    """
    Run optimisation in particular cluster and species type, under particular hypothesis
    :param cl_tag: tag of cluster
    :param dist_tag: tag of hypothesis
    :param pop_tag: desi/kabuli tag
    :return:
    """

    seeds = [239, 30, 123]  # Seeds for three chains

    for iseed in seeds:

        np.random.seed(seed=iseed)  # Set random seed

        file_ndata = path_data + 'n_data_{}.txt'.format(pop_tag)  # file with counts of alternative allele
        file_qdata = path_data + 'q_data_{}.txt'.format(pop_tag)  # file with total allele counts

        file_cov = path_cov + '{}_cov_{}.txt'.format(cl_tag, dist_tag)  # file with covariance matrix

        data = PopData(file_qdata, file_ndata)  # Dataset
        mx_cov = ClusterCov(file_cov)  # Covariance matrix
        opt_hist = HistOpt(mx_cov, data, ab_flag=True, ilr_flag=True)  # Optimizer

        path_res = path_mcmc + '{}_{}_{}_{}/'.format(cl_tag, dist_tag, pop_tag, iseed)  # Path for MCMC
        print(path_res)

        opt_hist.optimise(path_res, n_thr=30)  # Run optimiser


# ------------------------------------------------
# PATHS
# ------------------------------------------------
path_pref = 'data/'
path_data = path_pref + 'samples/'
path_cov = path_pref  + 'cov_mx/'

path_mcmc = path_pref + 'mcmc/'

if not os.path.exists(path_mcmc):
    os.mkdir(path_mcmc)

# ------------------------------------------------
# TAGS
# ------------------------------------------------
cluster_tags = ['ETHI', 'IND', 'LEB', 'MOR', 'TUR', 'UZB']
dist_tags = ['routes', 'lin']
pop_tags = ['desi', 'kabuli']

# Combinations to skip: we do not have Kabulis in Ethiopid and India
skip_combinations = [('kabuli', 'ETHI', 'lin'),
                    ('kabuli', 'ETHI', 'routes'),
                    ('kabuli', 'IND', 'lin'),
                    ('kabuli', 'IND', 'routes')]

# ------------------------------------------------
# OPTIMISATION
# ------------------------------------------------
for pop_tag, cl_tag, dist_tag in product(pop_tags, cluster_tags, dist_tags):

    if (pop_tag, cl_tag, dist_tag) in skip_combinations:
        continue
    run_freq_bayes(cl_tag, dist_tag, pop_tag)


