__author__ = 'sibirrer'

import numpy as np
from cosmoHammer import ParticleSwarmOptimizer
from cosmoHammer import LikelihoodComputationChain
from cosmoHammer import CosmoHammerSampler
from cosmoHammer.util import InMemoryStorageUtil

import astrofunc.util as util
from astrofunc.util import Util_class
from astrofunc.LensingProfiles.gaussian import Gaussian
from astrofunc.LightProfiles.moffat import Moffat


class Chain(object):
    """
    contains the routines to be fitted by a mcmc, meant for PSF estimations/re-center
    """
    def __init__(self, image, sigma, poisson, sampling_option, deltaPix=1, mask=None, subgrid_res=1,
                 x_grid=None, y_grid=None):
        """
        initializes all the classes needed for the chain
        """
        self.image = util.image2array(image)
        self.numPix = len(image)
        self.deltaPix = deltaPix
        self.subgrid_res = subgrid_res
        self.background = sigma
        self.poisson = poisson
        if x_grid is None or y_grid is None:
            self.x_grid, self.y_grid = util.make_grid(self.numPix, deltaPix, subgrid_res)
        else:
            util_class = Util_class()
            self.x_grid, self.y_grid = util_class.make_subgrid(x_grid, y_grid, subgrid_res)
        self.gaussian = Gaussian()
        self.moffat = Moffat()
        self.sampling_option = sampling_option
        if not mask is None:
            self.mask = util.image2array(mask)
        else:
            self.mask = np.ones((self.numPix, self.numPix))

    def X2_chain_gaussian(self, args):
        """
        routine to compute X2 given variable parameters for a MCMC/PSO chain
        """
        amp = args[0]
        sigma = np.exp(args[1])
        center_x = args[2]
        center_y = args[3]
        model = self.gaussian.function(self.x_grid, self.y_grid, amp, sigma, sigma, center_x, center_y)
        X2 = util.compare(model, self.image, self.background, self.poisson)
        return -X2, None

    def X2_chain_moffat(self, args):
        """
        routine to compute X2 given variable parameters for a MCMC/PSO chain
        """
        amp = args[0]
        alpha = args[1]
        beta = args[2]
        center_x = args[3]
        center_y = args[4]
        model = self.moffat.function(self.x_grid, self.y_grid, amp, alpha, beta, center_x, center_y)
        X2 = util.compare(model, self.image, self.background, self.poisson)
        return -X2, None

    def __call__(self, a):
        if self.sampling_option == 'psf_gaussian':
            return self.X2_chain_gaussian(a)
        if self.sampling_option == 'psf_moffat':
            return self.X2_chain_moffat(a)

    def computeLikelihood(self, ctx):
        if self.sampling_option == 'psf_gaussian':
            likelihood, _ = self.X2_chain_gaussian(ctx.getParams())
        elif self.sampling_option == 'psf_moffat':
            likelihood, _ = self.X2_chain_moffat(ctx.getParams())
        else:
            raise ValueError('sampling option %s not valid!')
        return likelihood

    def setup(self):
        pass


class Fitting(object):
    """
    contains the fitting routines with act on class Chain
    """
    def gaussian_fit(self, image, mean, sigma, poisson, n_walker=100, n_iter=100, threadCount=1):
        """
        fits a gaussian psf to an image
        :param image:
        :return:a
        """
        image_new = image - mean
        # amp, sigma, center_x, center_y
        mean_start = np.array([np.max(image_new), 0, 0, 0])
        lowerLimit = np.array([0, -4, -2, -2])
        upperLimit = np.array([np.max(image_new)*100, 4, 2, 2])
        sigma_start = np.array([np.max(image_new), 1, 0.1,0.1])
        chain = Chain(image_new, sigma, poisson, 'psf_gaussian')
        mean_list, particles = self.pso(n_particles=n_walker, n_iterations=n_iter, lowerLimit=lowerLimit, upperLimit=upperLimit, X2_chain=chain, threadCount=threadCount)
        return mean_list

    def moffat_fit(self, image, mean, sigma, poisson, n_walker=100, n_iter=100, threadCount = 1):
        """
        fits a moffat profile to an image
        :param image:
        :return:
        """
        image_new = image - mean
        #amp, alpha, beta, center_x, center_y
        mean_start = np.array([np.max(image_new), 2, 0, 0, 0])
        lowerLimit = np.array([0, 0, -4, -2, -2])
        upperLimit = np.array([np.max(image_new), 10, 4, 2, 2])
        sigma_start = np.array([np.max(image_new), 1, 1, 0.1, 0.1])
        chain = Chain(image_new, sigma, poisson, 'psf_moffat')
        mean_list, particles = self.pso(n_particles=n_walker, n_iterations=n_iter, lowerLimit=lowerLimit, upperLimit=upperLimit, X2_chain=chain, threadCount=threadCount)
        return mean_list

    def mcmc_CH(self, walkerRatio, n_run, n_burn, mean_start, sigma_start, lowerLimit, upperLimit, X2_chain, threadCount=1, init_pos=None):
        """
        runs mcmc on the parameter space given parameter bounds with CosmoHammerSampler
        returns the chain
        """
        params = np.array([mean_start, lowerLimit, upperLimit, sigma_start]).T

        chain = LikelihoodComputationChain(
            min=lowerLimit,
            max=upperLimit)

        # chain.addCoreModule(CambCoreModule())
        chain.addLikelihoodModule(X2_chain)
        chain.setup()

        store = InMemoryStorageUtil()
        sampler = CosmoHammerSampler(
            params= params,
            likelihoodComputationChain=chain,
            filePrefix="testHammer",
            walkersRatio=walkerRatio,
            burninIterations=n_burn,
            sampleIterations=n_run,
            threadCount=threadCount,
            initPositionGenerator=init_pos,
            storageUtil=store)

        sampler.startSampling()

        if sampler._sampler.pool is not None:
            sampler._sampler.pool.close()

        return store.samples

    def pso(self, n_particles, n_iterations, lowerLimit, upperLimit, X2_chain, threadCount = 1, print_positions=False):
        """
        returns the best fit for the lense model on catalogue basis with particle swarm optimizer
        """


        pso = ParticleSwarmOptimizer(X2_chain, lowerLimit,  upperLimit, n_particles, threads=threadCount)

        X2_list = []
        vel_list = []
        pos_list = []

        for swarm in pso.sample(n_iterations):
            X2_list.append(pso.gbest.fitness*2)
            vel_list.append(pso.gbest.velocity)
            pos_list.append(pso.gbest.position)


        print(pso.gbest.fitness/np.sum(X2_chain.mask), 'reduced X^2 of best position')
        if print_positions == True:
            print(pso.gbest.position)
        result = tuple(pso.gbest.position)
        return result, [X2_list, pos_list, vel_list]

    def fit_sample(self, star_list, mean, sigma, poisson, n_walk=100, n_iter=100, threadCount=1, psf_type='gaussian'):
        """
        routine to fit a sample of several stars and to show the variation
        :param mean:
        :param sigma:
        :param poisson:
        :param walkerRatio:
        :param n_run:
        :param n_burn:
        :param threadCount:
        :return:
        """
        n = len(star_list) # number of different objects
        if psf_type == 'gaussian'or psf_type == 'pixel':
            numParam = 4
        elif psf_type == 'moffat':
            numParam = 5
        else:
            raise ValueError('type %s is not a valid input' %(type))
        mean_list = np.zeros((n, numParam))
        for i in range(n):
            image = star_list[i]
            image = util.cut_edges(image, 31)
            if psf_type == 'gaussian'or psf_type == 'pixel':
                mean_list[i] = self.gaussian_fit(image, mean, sigma, poisson,n_walker=n_walk, n_iter=n_iter, threadCount = threadCount)
            elif psf_type == 'moffat':
                mean_list[i] = self.moffat_fit(image, mean, sigma, poisson, n_walker=n_walk, n_iter=n_iter, threadCount = threadCount)
        return mean_list
