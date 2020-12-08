"""
Get statistics
"""
import numpy as np
from itertools import product
# from pop_hist_opt import HistOpt, PopData, ClusterCov
from pop_hist_opt_wright_fisher import HistOpt, PopData, ClusterCov


def run_freq_bayes(cl_tag, dist_tag, pop_tag):
    """
    Run optimisation in particular cluster and species type, under particular hypothesis
    :param cl_tag: tag of cluster
    :param dist_tag: tag of hypothesis
    :param pop_tag: desi/kabuli tag
    :return:
    """

    seeds = [239, 30, 123]  # Seeds for two chains

    for iseed in seeds:

        np.random.seed(seed=iseed)  # Set random seed

        file_ndata = path_data + 'n_data_{}.txt'.format(pop_tag)  # file with counts of alternative allele
        file_qdata = path_data + 'q_data_{}.txt'.format(pop_tag)  # file with total allele counts

        file_cov = path_cov + '{}_cov_{}.txt'.format(cl_tag, dist_tag)  # file with covariance matrix

        data = PopData(file_qdata, file_ndata)  # Dataset
        mx_cov = ClusterCov(file_cov)  # Covariance matrix
        opt_hist = HistOpt(mx_cov, data, ab_flag=True, ilr_flag=True)  # Optimizer

        path_res = path_mcmc + '{}_{}_{}_{}/'.format(cl_tag, dist_tag, pop_tag, iseed)  # Path for MCMC

        result_sum = 0
        result_ll = 0
        mean_f_a = []
        for snp_id in range(opt_hist.n_snp):
            # print(snp_id)

            k1 = opt_hist.q[:, snp_id]
            n1 = opt_hist.n[:, snp_id]

            # In the number of alternative alleles in all locations is zero - return
            if sum(k1) == 0:
                mean_f_a += [0]
                continue

            if sum(k1) == sum(n1):
                mean_f_a += [1]
                continue

            if (opt_hist.n_loc - sum(k1 == 0)) <= 2:
                mean_f_a += [k1.sum() / n1.sum()]
                continue

            if (opt_hist.n_loc - sum(k1 == n1)) <= 2:
                mean_f_a += [k1.sum() / n1.sum()]
                continue


            x0 = np.loadtxt(path_res + 'mean' + str(snp_id) + '.txt')
            if x0[0] == 0:
                print('No optimisation', snp_id)

            f_a = 1 / (1 + np.exp(x0[0]))
            mean_f_a += [f_a]

            res = opt_hist.logprob(x0, snp_id)
            result_ll += res[0]


            x_pop = x0[0] + x0[1:(opt_hist.n_loc + 1)]
            res = opt_hist.likelihood_snp(snp_id, x_pop)
            result_sum += res

        np.savetxt(fname=path_estim + '{}_{}_{}_{}'.format(cl_tag, dist_tag, pop_tag, iseed) + '_freq_a.txt',
                   X=mean_f_a, fmt='%.10f')
        #print(len(mean_f_a))
        print(cl_tag, dist_tag, pop_tag, result_sum, result_ll, iseed)
        # opt_hist.optimise(path_res, n_thr=30)  # Run optimiser


# ------------------------------------------------
# PATHS
# ------------------------------------------------
path_pref = 'data/'
path_data = path_pref + 'samples/'
path_cov = path_pref  + 'cov_mx/'

path_mcmc = path_pref + 'mcmc/'
path_estim = path_pref + 'estim/'  # Path for estimates

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

k = 0
pop_tag = pop_tags[k]
cl_tag = cluster_tags[k]
dist_tag = dist_tags[k]


# ------------------------------------------------
# OPTIMISATION
# ------------------------------------------------
for pop_tag, cl_tag, dist_tag in product(pop_tags, cluster_tags, dist_tags):

    if (pop_tag, cl_tag, dist_tag) in skip_combinations:
        continue
    run_freq_bayes(cl_tag, dist_tag, pop_tag)


