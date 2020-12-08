"""
Tools to test a species spread hypothesis within a cluster
Input data:
        1) covariance matrix between locations -
            this covariance matrix is obtained inline with the hypothesis
        2.1) counts of alternative allele in each location
        2.2) total number number of alleles in each location
"""

import os
import numpy as np
import scipy.stats as st

from functools import partial
from pandas import read_csv
from multiprocess import Pool

from pyhmc import hmc


def ilr(f):
    """
    ilr transformation of frequencies
    :param f: frequency
    :return: ilr-transformation
    """
    return np.log((1 - f) / f)


def inv_ilr(x):
    """
    Inverce ilr transformation
    :param x: a coordinate within ilr space
    :return: a corresponding frequency
    """
    return 1 / (1 + np.exp(x))


class PopData:
    """
    The dataset class: everything about SNPs in all locations
    """
    def __init__(self, file_qdata, file_ndata):
        """
        :param file_qdata: file with counts of alternative allele in each location
        :param file_ndata: file with the total number of alternative allele in each location
        """
        self.q = read_csv(file_qdata, sep='\t',index_col=0)  # Counts
        self.n = read_csv(file_ndata, sep='\t',index_col=0)  # Total number
        self.loc_name = list(self.q.index)  # names of loci
        self.n_loc, self.n_snp = self.q.shape  # number of locations and SNPs


class ClusterCov:
    """
    Covariance matrix between locations based on the hypothesis
    """
    def __init__(self, file_cov):
        """
        :param file_cov: File with covariance matrix
        """
        self.cov = read_csv(file_cov, sep='\t',index_col=0)  # Covariance matrix
        self.loc_name = list(self.cov.index)  # names of locations
        self.n_loc = self.cov.shape[0]  # number of loci

        # We normalised covariance matrix to the mean of diagonal elements.
        # The diagonal elements are variances of each locations according to the ancestral state,
        # which, in our liner/routes model, can be considered as
        # lengths of branches in the binary tree from the root to leaves.
        # To normalise these distances, we normalise all variances and covariances
        # to the mean distance from the root to leaves.
        self.cov /= np.mean(np.diag(self.cov))


class HistOpt:
    """
    Estimator of alternative allele frequencies within each location and in the ancestral state
    """
    def __init__(self, mx_cov: ClusterCov, data: PopData, ab_flag=False, ilr_flag=True):
        """
        :param mx_cov: covariance matrix
        :param data: the dataset
        :param ab_flag: whether the prior is used for frequencies in ancestral state
        :param ilr_flag: whether to use the ilr transformation
        """

        self.ab_flag = ab_flag
        self.ilr_flag = ilr_flag

        # shared names of locations in covariance matrix and the dataset
        self.loc_name = [y for y in mx_cov.loc_name
                         if y in data.loc_name]


        # Dataset
        self.q = data.q.loc[self.loc_name, :]  # alt.allele counts
        self.n = data.n.loc[self.loc_name, :]  # total allele counts
        self.q = self.q.as_matrix()
        self.n = self.n.as_matrix()

        self.n_loc = len(self.loc_name)  # number of locations
        self.n_snp = data.n_snp  # number of SNPs

        # freq_anc, freq_pop = self.set_init_freqs()  # sample frequencies of alternative alleles

        # self.x_anc = ilr(freq_anc)  # ilr-transformed frequencies in locations
        # self.x_pop = ilr(freq_pop)  # ilr-transformed frequencies in ancestral state


#       # Covariance matrix
        self.m_c = mx_cov.cov.loc[self.loc_name, self.loc_name]
        self.m_c = self.m_c.as_matrix()
        self.m_c_inv = np.linalg.inv(self.m_c)  # Inverse covariance matrix

        # # Priors for proportionality coefficients
        # self.sigma = 0.01 * np.ones(self.n_snp)
        # self.sigma_mean = 0.001
        # self.sigma_var = 0.0001

        # # Priors for frequency in the ancestral state
        # if self.ab_flag:
        #     self.a = np.ones(self.n_snp)
        #     self.b = np.ones(self.n_snp)
        # self.step = 0.0001

        # mcmc chain
        # self.chain_freq_opt = []
        # self.count = 0


    def optimise(self, path_data, n_thr=2):
        """
        Optimisation of parameters
        :param path_data: directory where to save the data
        :param freqs_init: initial frequencies
        :param sigma_params: parameters of disribution of sigma
        :param n_thr: number of threads
        :return:
        """

        if not os.path.exists(path_data):  # make dir
            os.makedirs(path_data)

        # self.load_init(freqs_init)  # load initial dataset

        func = partial(self.optimize_snp_new, path_chain=path_data)

        with Pool(n_thr) as workers:  # run the optimisation for each SNP separately
            pmap = workers.map
            # pmap(func, range(self.n_snp))
            pmap(func, range(1))




    def optimize_snp_new(self, snp_id, path_chain):
        """
        Optimisation of parameters for a SNPs by Hamiltonian Monte Carlo
        :param snp_id: Id of SNP
        :param path_chain: path to save MCMC chain
        """

        # snp_id = 6

        # print(snp_id)

        k1 = self.q[:, snp_id]
        n1 = self.n[:, snp_id]


        # In the number of alternative alleles in all locations is zero - return
        if sum(k1) == 0:
            return

        if sum(k1) == sum(n1):
            return

        if (self.n_loc - sum(k1 == 0)) <= 2:
            # print('Skip zero', path_chain, snp_id)
            return

        if (self.n_loc - sum(k1 == n1)) <= 2:
            # print('Skip non-zero', path_chain, snp_id)
            return


        # niter = 50000
        # samples = hmc(self.logprob, x0=[0] + list(np.random.uniform(-1, 1, self.n_loc)) + [0.1] + [1] * 2,
        #               args=(snp_id,), n_samples=2000, epsilon=0.1)
        # print(samples[-1])
        # samples = hmc(self.logprob, x0=samples[-1],
        #               args=(snp_id,), n_samples=niter, epsilon=0.2)




        # print(list(zip(k1, n1)))

        # Burn-in
        niter = 5000
        samples = [[k1.sum()/n1.sum()] + list(np.random.uniform(-1, 1, self.n_loc)) + [0.1] + [1] * 2]

        params = np.array(samples[-1])
        print(self.logprob(np.array(samples[-1]), snp_id))

        eps = 0.11
        while len(samples) < (niter / 4):
            eps = eps - 0.01
            samples = hmc(self.logprob, x0=samples[-1],
                          args=(snp_id,), n_samples=niter, epsilon=eps)
            _, idx = np.unique(samples, axis=0, return_index=True)
            samples = samples[np.sort(idx)]


        print(len(samples), eps)


        # Find epsilon after burn-in
        niter = 5000
        eps = 0.18
        x0 = samples[-1]
        samples = []

        while len(samples) < (niter / 4) and eps > 0:
            if eps > 0.019:
                eps -= 0.01
            else:
                eps = eps * 0.7
            samples = hmc(self.logprob,
                          x0=x0,
                          args=(snp_id,), n_samples=niter, epsilon=eps)
            _, idx = np.unique(samples, axis=0, return_index=True)
            samples = samples[np.sort(idx)]

            print(len(samples), eps)


        niter = 100000
        if eps > 0:
            samples = hmc(self.logprob,
                          x0=samples[-1],
                          args=(snp_id,), n_samples=niter, epsilon=eps)
            _, idx = np.unique(samples, axis=0, return_index=True)
            samples = samples[np.sort(idx)]
        else:
            print('Bad epsilon')



        print('Len', path_chain, snp_id, len(samples), eps)

        if len(samples) < niter/4:
            print('Low number of samples', path_chain, snp_id)



        np.savetxt(fname=path_chain +
                   'snp' + str(snp_id) + '.txt', X=samples, fmt='%.10f')

        np.savetxt(fname=path_chain +
                         'eps' + str(snp_id) + '.txt', X=eps, fmt='%.10f')

        n_thin = 10
        samples = samples[::n_thin]
        np.savetxt(fname=path_chain +
                         'mean' + str(snp_id) + '.txt', X=samples.mean(0), fmt='%.10f')

        if 0 in samples.var(0):
            print('Bad News', snp_id)
        np.savetxt(fname=path_chain +
                         'var' + str(snp_id) + '.txt', X=samples.var(0), fmt='%.10f')


    def likelihood_snp(self, snp_id, x_pop):
        """
        Lod-likelihood of the data in one SNP: Binomial distribution
        :param snp_id: SNP id
        :return: Binimoal log-likelood
        """
        # p = 0
        # for loc_id in range(self.n_loc):
        #     if self.n[loc_id, snp_id] > 0:
        #
        #         freq_pop = 1 / (1 + np.exp(self.x_pop[loc_id, snp_id]))
        #         p = st.binom.logpmf(n=self.n[loc_id, snp_id],
        #                              k=self.q[loc_id, snp_id],
        #                              p=freq_pop)
        #         print(p)
        # return p


        freq_pop = 1 / (1 + np.exp(x_pop))
        tmp = self.q[:, snp_id] * np.log(freq_pop) + \
             (self.n[:, snp_id] - self.q[:, snp_id]) * np.log(1 - freq_pop)
        # print(tmp)
        return tmp.sum()

    # ==========================================================================
    # WORK WITH FREQUENCIES and ILR COORDINATES
    # ==========================================================================

    # def set_init_freqs(self):
    #     """
    #
    #     :return:
    #     """
    #     q_sum = reduce(lambda x, y: x + y, self.q)
    #     n_sum = reduce(lambda x, y: x + y, self.n)
    #
    #     freq_mean = q_sum / n_sum * 0 + 0.5
    #     freq_anc = freq_mean
    #     freq_pop = np.tile(freq_mean, (self.n_loc, 1))
    #
    #     return freq_anc, freq_pop


    # ==========================================================================
    # HAMILTONIAN MONTE CARLO: log-posteriors and gradients
    # ==========================================================================

    def logprob(self, params, snp_id):  # Workable

        k1 = self.q[:, snp_id]
        n1 = self.n[:, snp_id]

        n_pop = len(k1)

        x_a = params[0]  # ILR of ancestral state
        r_p = params[1:(n_pop + 1)]  # any an additional component

        z = params[n_pop + 1]
        a = params[n_pop + 2]
        b = params[n_pop + 3]

        x_p = x_a + r_p
        f_p = 1 / (1 + np.exp(x_p))  # frequency in populations
        f_a = 1 / (1 + np.exp(x_a))  # frequency in ancestral state


        q1 = (r_p).T.dot(self.m_c_inv) * z
        quad_form = (q1 * (r_p).T).sum() / 2

        logp = (k1 * np.log(f_p) + (n1 - k1) * np.log(1 - f_p)).sum()
        logp -= quad_form * f_a * (1 - f_a)
        # logp += (a - 1) * np.log(f_a) + (a / m * (1-m) - 1) * np.log(1 - f_a)
        logp += (b - 1 + n_pop / 2) * x_a - (n_pop + a + b - 2) * np.log(1 + np.exp(x_a))

        barrier = 10
        # barrier = 15
        logp += (np.log(barrier - np.abs(r_p))).sum()

        # additional for z param
        z_lambda = 1
        logp += n_pop / 2 * np.log(z) - z_lambda * z

        # additional for "a" param
        logp += st.uniform.logpdf(a / (a + b)) + st.expon.logpdf(a + b)

        # grad_bin = (n1 / (1 + np.exp(x_p))) - k1
        grad_bin = n1 * f_p - k1

        grad = np.concatenate(([sum(grad_bin) +
                                - (-1) * quad_form * (np.exp(x_a) - np.exp(-x_a)) / (
                                            np.exp(x_a) + 2 + np.exp(-x_a)) ** 2 +
                                b - 1 + n_pop / 2 - (n_pop + a + b - 2) * (1 - f_a)],
                               grad_bin - q1 * f_a * (1 - f_a) + (-1) * np.sign(r_p) / (barrier - np.abs(r_p)),
                               [-quad_form * f_a * (1 - f_a) / z + n_pop / 2 / z - z_lambda],
                               [np.log(f_a) - 1],
                               [np.log(1 - f_a) - 1]))

        return logp, grad