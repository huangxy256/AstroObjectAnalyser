__author__ = 'sibirrer'

import numpy as np


from cosmoHammer import ParticleSwarmOptimizer
from cosmoHammer import LikelihoodComputationChain
from cosmoHammer import CosmoHammerSampler
from cosmoHammer.util import InMemoryStorageUtil


import lenstronomy.util as util
from lenstronomy.util import Util_class
from lenstronomy.FunctionSet.gaussian import Gaussian
from lenstronomy.FunctionSet.moffat import Moffat
from lenstronomy.FunctionSet.sersic import Sersic_elliptic, Sersic, CoreSersic
from lenstronomy.ImSim.make_image import MakeImage
from lenstronomy.Workflow.parameters import Param

class Chain(object):
    """
    contains the routines to be fitted by a mcmc
    """
    def __init__(self, image, sigma, poisson, sampling_option, deltaPix=1, mask=None, subgrid_res=1, psf_kwargs={}, x_grid=None, y_grid=None):
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
        self.sersic_elliptic = Sersic_elliptic()
        self.sersic_core = CoreSersic()
        self.sersic = Sersic()
        self.sampling_option = sampling_option
        if not mask is None:
            self.mask = util.image2array(mask)
        else:
            self.mask = np.ones((self.numPix, self.numPix))
        self.makeImage = MakeImage(kwargs_data={}, kwargs_options={'lens_type': 'GAUSSIAN', 'source_type': 'GAUSSIAN', 'psf_type':'pixel'})
        self.psf_kwargs = psf_kwargs
        kwargs_options = {'lens_type': 'SPEP', 'source_type': 'GAUSSIAN', 'psf_type':'pixel'}
        self.param = Param(kwargs_options)

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

    def X2_chain_moffat(self,args):
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

    def X2_chain_sersic(self, args):
        """
        routineto compute X2 given parameters for MCMC/PSO chain
        """
        model, penalty = self.get_model_sersic(args)
        model = util.image2array(model)
        X2 = util.compare(model*self.mask, self.image*self.mask, self.background, self.poisson)
        return -X2 - penalty, None

    def X2_chain_sersics2(self, args):
        """
        routineto compute X2 given parameters for MCMC/PSO chain
        """
        model, penalty = self.get_model_sersic2(args)
        model = util.image2array(model)
        X2 = util.compare(model*self.mask, self.image*self.mask, self.background, self.poisson)
        return -X2 - penalty, None

    def get_model_sersic2(self, args):
        penalty = 0
        I0_sersic = np.exp(args[0])
        R_sersic = args[1]
        n_sersic = args[2]
        phi_G = args[3]
        q = args[4]
        center_x = args[5]
        center_y = args[6]
        I0_sersic2 = np.exp(args[7])
        R_sersic2 = args[8]
        n_sersic2 = 1
        center_x2 = args[9]
        center_y2 = args[10]
        model1 = self.sersic_elliptic.function(self.x_grid, self.y_grid, I0_sersic, R_sersic, n_sersic, phi_G, q, center_x, center_y)
        model2 = self.sersic.function(self.x_grid, self.y_grid, I0_sersic2, R_sersic2, n_sersic2, center_x2, center_y2)
        model = util.array2image(model1 + model2)
        model = self.makeImage.re_size_convolve(model, self.numPix, self.deltaPix, self.subgrid_res, self.psf_kwargs)
        return model, penalty

    def X2_chain_sersic_core(self, args):
        """
        routine to compute X2 given parameters for MCMC/PSO chain
        """
        model, penalty = self.get_model_sersic_core(args)
        model = util.image2array(model)
        X2 = util.compare(model*self.mask, self.image*self.mask, self.background, self.poisson)
        return -X2 - penalty, None

    def get_model_sersic(self, args):
        penalty = 0
        I0_sersic = args[0]
        R_sersic = args[1]
        n_sersic = args[2]
        phi_G = args[3]
        q = args[4]
        center_x = args[5]
        center_y = args[6]
        model = self.sersic_elliptic.function(self.x_grid, self.y_grid, I0_sersic, R_sersic, n_sersic, phi_G, q, center_x, center_y)
        model = util.array2image(model)
        model = self.makeImage.re_size_convolve(model, self.numPix, self.deltaPix, self.subgrid_res, self.psf_kwargs)
        return model, penalty

    def get_model_sersic_core(self,args):
        # I0_sersic, Rb, R_sersic, n_sersic, gamma, phi_G, q, center_x, center_y, I0_2, R_2, n_2=1, center_x_2, center_y_2
        penalty = 0
        I0_sersic = args[0]
        Rb = args[1]
        R_sersic = args[2]
        n_sersic = args[3]
        gamma = args[4]
        phi_G = args[5]
        q = args[6]
        center_x = args[7]
        center_y = args[8]
        I0_2 = np.exp(args[9])
        R_2 = args[10]
        n_2 = 1
        center_x2 = args[11]
        center_y2 = args[12]
        model1 = self.sersic_core.function(self.x_grid, self.y_grid, I0_sersic, Rb, R_sersic, n_sersic, gamma, phi_G, q, center_x, center_y)
        model2 = self.sersic.function(self.x_grid, self.y_grid, I0_2, R_2, n_2, center_x2, center_y2)
        model = util.array2image(model1 + model2)
        model = self.makeImage.re_size_convolve(model, self.numPix, self.deltaPix, self.subgrid_res, self.psf_kwargs)
        return model, penalty

    def X2_chain_sersic_triple(self, args):
        """
        routine to compute X2 given parameters for MCMC/PSO chain
        """
        model, penalty = self.get_model_sersic_triple(args)
        model = util.image2array(model)
        X2 = util.compare(model*self.mask, self.image*self.mask, self.background, self.poisson)
        return -X2 - penalty, None

    def get_model_sersic_triple(self, args):
        """
        routine to get double sersic + single sersic profile
        """
        penalty = 0
        I0_sersic = np.exp(args[0])
        R_sersic = args[1]
        n_sersic = args[2]
        phi_G = args[3]
        q = args[4]
        center_x = args[5]
        center_y = args[6]
        I0_sersic2 = np.exp(args[7])
        R_sersic2 = args[8]
        center_x_2 = args[9]
        center_y_2 = args[10]
        I0_3 = np.exp(args[11])
        R_3 = args[12]
        n_3 = args[13]
        n_2 = args[14]

        model1 = self.sersic_elliptic.function(self.x_grid, self.y_grid, I0_sersic, R_sersic, n_sersic, phi_G, q, center_x, center_y)
        model2 = self.sersic.function(self.x_grid, self.y_grid, I0_sersic2, R_sersic2, n_2, center_x_2, center_y_2)
        model3 = self.sersic_elliptic.function(self.x_grid, self.y_grid, I0_3, R_3, n_3, phi_G, q, center_x, center_y)
        model = util.array2image(model1 + model2 + model3)
        model = self.makeImage.re_size_convolve(model, self.numPix, self.deltaPix, self.subgrid_res, self.psf_kwargs)
        return model, penalty

    def __call__(self, a):
        if self.sampling_option == 'psf_gaussian':
            return self.X2_chain_gaussian(a)
        if self.sampling_option == 'psf_moffat':
            return self.X2_chain_moffat(a)
        if self.sampling_option == 'sersic_profile':
            return self.X2_chain_sersic(a)
        if self.sampling_option == 'sersic_profile2':
            return self.X2_chain_sersics2(a)
        if self.sampling_option == 'core_sersic':
            return self.X2_chain_sersic_core(a)
        if self.sampling_option == 'sersic_triple':
            return self.X2_chain_sersic_triple(a)

    def computeLikelihood(self, ctx):
        if self.sampling_option == 'psf_gaussian':
            likelihood, _ = self.X2_chain_gaussian(ctx.getParams())
        elif self.sampling_option == 'psf_moffat':
            likelihood, _ = self.X2_chain_moffat(ctx.getParams())
        elif self.sampling_option == 'sersic_profile':
            likelihood, _ = self.X2_chain_sersic(ctx.getParams())
        elif self.sampling_option == 'sersic_profile2':
            likelihood, _ = self.X2_chain_sersics2(ctx.getParams())
        elif self.sampling_option == 'sersic_core':
            likelihood, _ = self.X2_chain_sersic_core(ctx.getParams())
        elif self.sampling_option == 'sersic_triple':
            likelihood, _ = self.X2_chain_sersic_triple(ctx.getParams())
        else:
            raise ValueError('sampling option %s not valid!')
        return likelihood

    def setup(self):
        pass


class Fitting(object):
    """
    contains the fitting routines with act on class Chain
    """
    def gaussian_fit(self, image, mean, sigma, poisson, n_walker=100, n_iter=100, threadCount = 1):
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

    def sersic_fit(self, image, ra_coords, dec_coords, sigma, poisson, mask, deltaPix, psf_kwargs, n_particles, n_iterations, lowerLimit=None, upperLimit=None, threadCount=1):
        """
        fit a sersic profile with PSO method
        :param image:
        :param mean:
        :param sigma:
        :param poisson:
        :param n_walker:
        :param n_iter:
        :param threadCount:
        :return:
        """
        if lowerLimit is None or upperLimit is None:
            #I0_sersic 0,R_sersic 1, n_sersic 2, phi_G 3, q 4, center_x 5, center_y 6
            mean_start = np.array([1,1,2,np.pi, 0.8, 0, 0])
            lowerLimit = np.array([0, 0, 0, 0, 0.2,-0.5,-0.5])
            upperLimit = np.array([100,5,5,2*np.pi,1,0.5,0.5])

        chain = Chain(image, sigma, poisson, 'sersic_profile', deltaPix=deltaPix, mask=mask, subgrid_res=1, psf_kwargs=psf_kwargs, x_grid=ra_coords, y_grid=dec_coords)
        result, particles = self.pso(n_particles, n_iterations, lowerLimit=lowerLimit, upperLimit=upperLimit, X2_chain=chain, threadCount=threadCount)
        kwargs = {'I0_sersic': result[0], 'R_sersic': result[1], 'n_sersic': result[2], 'phi_G': result[3]
            , 'q': result[4], 'center_x': result[5], 'center_y': result[6]}
        model, _ = chain.get_model_sersic(result)
        return kwargs, model, particles

    def sersic2_fit(self, image, ra_coords, dec_coords, sigma, poisson, mask, deltaPix, psf_kwargs, n_particles=120, n_iterations=400, lowerLimit=None, upperLimit=None, threadCount=6):
        """
        fit a sersic profile with PSO method
        :param image:
        :param mean:
        :param sigma:
        :param poisson:
        :param n_walker:
        :param n_iter:
        :param threadCount:
        :return:
        """
        #I0_sersic 0,k_sersic 1, n_sersic 2, phi_G 3, q 4, center_x 5, center_y 6, I0_sersic2 7,k_sersic2 8, center_x2 9, center_y2 10
        x_mins, y_mins, values = util.neighborSelect(util.image2array(-image*mask), ra_coords, dec_coords)
        if lowerLimit is None and upperLimit is None:
            i1 = np.where(values == min(values))[0]
            x_center1 = x_mins[i1]
            y_center1 = y_mins[i1]
            if i1 == 0:
                i2 = 1
            else:
                i2 = 0
            x_center2 = x_mins[i2]
            y_center2 = y_mins[i2]

            mean_start = np.array([3,2,2,np.pi, 0.8, x_center1, y_center1, 1,0.1, x_center2, y_center2])
            sigma_array = np.array([2,1,1,1,0.2,0.01,0.01,2,0.05,0.01,0.01])
            lowerLimit = mean_start-sigma_array
            upperLimit = mean_start+sigma_array
        chain = Chain(image, sigma, poisson, 'sersic_profile2', deltaPix=deltaPix, mask=mask, subgrid_res=2, psf_kwargs=psf_kwargs, x_grid=ra_coords, y_grid=dec_coords)
        result, particles = self.pso(n_particles, n_iterations, lowerLimit=lowerLimit, upperLimit=upperLimit, X2_chain=chain, threadCount=threadCount)
        kwargs_lens_light = {'I0_sersic': np.exp(result[0]), 'R_sersic': result[1], 'n_sersic': result[2], 'phi_G': result[3]
            , 'q': result[4], 'center_x': result[5], 'center_y': result[6], 'I0_2': np.exp(result[7]), 'R_2': result[8], 'n_2': 1
            , 'center_x_2': result[9], 'center_y_2': result[10]}
        model, _ = chain.get_model_sersic2(result)
        return kwargs_lens_light, model, particles

    def core_sersic_fit(self, image, ra_coords, dec_coords, sigma, poisson, mask, deltaPix, psf_kwargs, n_particles=120, n_iterations=400, lowerLimit=None, upperLimit=None, threadCount=6):
        """
        fit a core sersic and a spehrical sersic profile
        :return:
        """
        if lowerLimit is None or upperLimit is None:
            # I0_sersic, Rb, R_sersic, n_sersic, gamma, phi_G, q, center_x, center_y, I0_2, R_2, n_2=1, center_x_2, center_y_2
            x_mins, y_mins, values = util.neighborSelect(util.image2array(-image*mask), ra_coords, dec_coords)
            i1 = np.where(values == min(values))[0]
            x_center1 = x_mins[i1]
            y_center1 = y_mins[i1]
            if i1 == 0:
                i2 = 1
            else:
                i2 = 0
            x_center2 = x_mins[i2]
            y_center2 = y_mins[i2]
            mean_start = np.array([2, 0.1, 1, 1, 1, 0, 0.8, x_center1, y_center1, 1, 0.1, x_center2, y_center2])
            sigma_array = np.array([2, 0.1, 0.5, 0.5, 1, 1, 0.1, 0.05, 0.05, 2, 0.05, 0.01, 0.01])
            lowerLimit = mean_start-sigma_array
            upperLimit = mean_start+sigma_array
        chain = Chain(image, sigma, poisson, 'core_sersic', deltaPix=deltaPix, mask=mask, subgrid_res=2, psf_kwargs=psf_kwargs, x_grid=ra_coords, y_grid=dec_coords)
        result, particles = self.pso(n_particles, n_iterations, lowerLimit=lowerLimit, upperLimit=upperLimit, X2_chain=chain, threadCount=threadCount)
        model,_ = chain.get_model_sersic_core(result)
        kwargs_lens_light = {'I0_sersic': result[0], 'Rb': result[1], 'R_sersic': result[2]
            , 'n_sersic': result[3], 'gamma': result[4], 'phi_G': result[5], 'q': result[6], 'center_x': result[7]
            , 'center_y': result[8], 'I0_2': result[9], 'R_2': result[10], 'n_2': 1., 'center_x_2': result[11], 'center_y_2': result[12]}
        return kwargs_lens_light, model, particles

    def sersic_triple_fit(self, image, ra_coords, dec_coords, sigma, poisson, mask, deltaPix, psf_kwargs, n_particles=120, n_iterations=200, lowerLimit=None, upperLimit=None, threadCount=6):
        """
        fit a core sersic and a spehrical sersic profile
        :return:
        """
        # I0_sersic, R_sersic, n_sersic, phi_G, q, center_x, center_y, I0_2, R_2, n_2=1, center_x_2, center_y_2, I0_3, R_3, n_3
        if lowerLimit is None or upperLimit is None:
            x_mins, y_mins, values = util.neighborSelect(util.image2array(-image*mask), ra_coords, dec_coords)
            i1 = np.where(values == min(values))[0]
            x_center1 = x_mins[i1]
            y_center1 = y_mins[i1]
            if i1 == 0:
                i2 = 1
            else:
                i2 = 0
            try:
                x_center2 = x_mins[i2]
                y_center2 = y_mins[i2]
            except:
                x_center2 = x_center1
                y_center2 = y_center1
                print('no second minima found, choosing the same start position for the center.')
            mean_start = np.array([1,2,3,np.pi, 0.8, x_center1, y_center1, 1,0.1, x_center2, y_center2, 1, 0.5, 3, 1])
            sigma_array = np.array([1,1,1,1,0.2,0.01,0.01,2,0.05,0.01,0.01,1,0.1,1,0.1])
            lowerLimit = mean_start-sigma_array
            upperLimit = mean_start+sigma_array
        chain = Chain(image, sigma, poisson, 'sersic_triple', deltaPix=deltaPix, mask=mask, subgrid_res=2, psf_kwargs=psf_kwargs, x_grid=ra_coords, y_grid=dec_coords)
        result, particles = self.pso(n_particles, n_iterations, lowerLimit=lowerLimit, upperLimit=upperLimit, X2_chain=chain, threadCount=threadCount)
        model, _ = chain.get_model_sersic_triple(result)
        kwargs_lens_light = {'I0_sersic': np.exp(result[0]), 'R_sersic': result[1]
            , 'n_sersic': result[2], 'phi_G': result[3], 'q': result[4], 'center_x': result[5]
            , 'center_y': result[6], 'I0_2': np.exp(result[7]), 'R_2': result[8], 'center_x_2': result[9], 'center_y_2': result[10]
            , 'I0_3': np.exp(result[11]), 'R_3': result[12], 'n_3': result[13], 'n_2': result[14]}
        return kwargs_lens_light, model, particles

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
        result=tuple(pso.gbest.position)
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
        mean_list = np.zeros((n,numParam))
        for i in range(n):
            image = star_list[i]
            image = util.cut_edges(image, 31)
            if psf_type == 'gaussian'or psf_type == 'pixel':
                mean_list[i] = self.gaussian_fit(image, mean, sigma, poisson,n_walker=n_walk, n_iter=n_iter, threadCount = threadCount)
            elif psf_type == 'moffat':
                mean_list[i] = self.moffat_fit(image, mean, sigma, poisson, n_walker=n_walk, n_iter=n_iter, threadCount = threadCount)
        return mean_list

    # def fit_sample_pso(self, star_list, mean, sigma, poisson, walkerRatio=10, n_run= 1000, n_burn=1000, threadCount = 1, psf_type = 'gaussian'):