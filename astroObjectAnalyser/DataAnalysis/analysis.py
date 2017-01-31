__author__ = 'sibirrer'

import numpy as np
import scipy.ndimage.interpolation as interp

import astrofunc.util as util
from astrofunc.LightProfiles.moffat import Moffat
from astrofunc.LensingProfiles.gaussian import Gaussian

from astroObjectAnalyser.DataAnalysis.psf_fitting import Fitting
from astroObjectAnalyser.DataAnalysis.catalogues import Catalogue


class Analysis(Catalogue):
    """
    class which analyses data and fits a psf
    """

    def get_psf(self, image, cat, mean, rms, poisson, psf_type = 'moffat', restrict_psf=None):
        """
        fit a given psf model
        :param image: cutout image to fit a profile on
        :param mean: mean background level
        :param rms: rms value of background
        :param poisson: reduction factor for poisson noise
        :return: parameters of the psf model
        """
        kwargs_cut = self.estimate_star_thresholds(cat)
        mask = self.find_objects(cat, kwargs_cut)
        star_list = self.get_objects_image(image, cat, mask, cut_radius=10, cut_fixed=61)
        fitting = Fitting()
        mean_list, std_list = fitting.fit_sample(star_list, mean, rms, poisson, n_walk=50, n_iter=50, threadCount=1, psf_type=psf_type)
        kernel, mean_list, restrict_psf, star_list_shift = self.stacking(star_list, mean_list, mean, rms, poisson, psf_type, restrict_psf= restrict_psf)
        return kernel, mean_list, restrict_psf

    def get_psf_outside(self, HDUFile, image_no_border, exp_time, psf_type='moffat', restrict_psf=None, kwargs_cut={}):
        """
        fit a given psf model without saving it in the system class, dangerous!
        :param image: cutout image to fit a profile on
        :param mean: mean background level
        :param rms: rms value of background
        :param poisson: reduction factor for poisson noise
        :return: parameters of the psf model
        """
        mean, rms = self.get_background(HDUFile)
        cat = self.get_source_cat(HDUFile)
        if kwargs_cut == {}:
            kwargs_cut = self.estimate_star_thresholds(cat)
        mask = self.find_objects(cat, kwargs_cut)
        star_list = self.get_objects_image(image_no_border, cat, mask, cut_radius=10, cut_fixed=61)
        fitting = Fitting()
        mean_list = fitting.fit_sample(star_list, mean, rms, exp_time, n_walk=100, n_iter=100, threadCount=1, psf_type=psf_type)
        kernel, mean_list, restrict_psf, star_list_shift = self.stacking(star_list, mean_list, mean, rms, exp_time, psf_type, restrict_psf=restrict_psf)
        return kernel, mean_list, restrict_psf, star_list_shift

    def stacking(self, star_list, mean_list, mean, rms, poisson, psf_type, restrict_psf=None):
        """

        :param star_list:
        :return:
        """
        if restrict_psf == None:
            restrict_psf = self.add_criteria_2(star_list, mean_list, mean, rms, poisson)
        shifteds = []
        mean_list_select = []
        for i in range(len(star_list)):
            if restrict_psf[i] is True:
                data = star_list[i]-mean
                if psf_type == 'gaussian' or psf_type == 'pixel':
                    amp, sigma, center_x, center_y = mean_list[i]
                elif psf_type == 'moffat':
                    amp, alpha, beta, center_x, center_y = mean_list[i]
                else:
                    raise ValueError('psf type %s not valid', (psf_type))
                shifted = interp.shift(data, [-center_y, -center_x], order=2)
                shifteds.append(shifted)
                mean_list_select.append(mean_list[i])
                print('=== object ===', i, center_x, center_y)
                import matplotlib.pylab as plt
                fig, ax1 = plt.subplots()
                im=ax1.matshow(np.log10(shifted), origin='lower')
                plt.axes(ax1)
                fig.colorbar(im)
                plt.show()

        combined = sum(shifteds)
        mean_list_select = np.mean(mean_list_select[:])
        new=np.empty_like(combined)
        max_pix = np.max(combined)
        p = combined[combined >= max_pix/10**6]  #in the SIS regime
        new[combined < max_pix/10**6] = 0
        new[combined >= max_pix/10**6] = p
        kernel = util.kernel_norm(new)
        return kernel, mean_list_select, restrict_psf, shifteds

    def add_criteria(self, star_list, mean_list, mean, rms, exp_time):
        """
        additional criteria of whether one should use the object to estimate the psf.
        :param star_list:
        :param mean_list:
        :return: list of bool with length of star_list
        """
        gaussian = Gaussian()
        num = len(star_list)
        X2_list = np.zeros((num))
        restricted_list = []
        for i in range(num):
            data = star_list[i]
            numPix = 21
            data = util.cut_edges(data, numPix)
            amp, sigma, center_x, center_y = mean_list[i]
            x_grid, y_grid = util.make_grid(numPix,1.)
            model_1D = gaussian.function(x_grid, y_grid, amp, sigma, sigma, center_x, center_y)
            model = util.array2image(model_1D)
            X2_list[i] = np.sum((data - model - mean)**2/(model/exp_time+rms**2))
        X2_min = np.min(X2_list)
        for i in range(num):
            if X2_list[i] > 10**2 * X2_min:
                restricted_list.append(False)
            else:
                restricted_list.append(True)
        return restricted_list

    def add_criteria_2(self, star_list, mean_list, mean, rms, exp_time):
        num = len(star_list)
        restricted_list = []
        for i in range(num):
            restricted_list.append(True)
        return restricted_list

    def get_psf_kwargs_update(self, image_name, psf_type, psf_size=None, psf_size_large=91, filter_object=None, kwargs_cut={}):
        """
        does the same as get_psf_kwargs but can also restrict itself to specially chosen objects
        :param image_name:
        :param psf_type:
        :param psf_size:
        :param filter_object:
        :return:
        """
        exp_time = self.system.get_exposure_time(image_name)
        HDUFile, image_no_border = self.system.get_HDUFile(image_name)

        kernel_large, mean_list, restrict_psf, star_list = self.analysis.get_psf_outside(HDUFile, image_no_border, exp_time, psf_type, filter_object, kwargs_cut)
        if psf_type == 'gaussian':
            sigma = mean_list[1]*self.system.get_pixel_scale(image_name)
            psf_kwargs = {'psf_type': psf_type, 'sigma': sigma}
        elif psf_type == 'moffat':
            alpha = mean_list[1]
            beta = mean_list[2]
            alpha *= self.system.get_pixel_scale(image_name)
            psf_kwargs = {'psf_type': psf_type, 'alpha': alpha, 'beta': beta}
        elif psf_type == 'pixel':
            kernel = util.cut_edges(kernel_large, psf_size)
            kernel = util.kernel_norm(kernel)
            kernel_large = util.cut_edges(kernel_large, psf_size_large)
            kernel_large = util.kernel_norm(kernel_large)
            kernel_list = []
            for i in range(len(star_list)):
                if i == 0:
                   kernel_list.append(kernel_large)
                else:
                    star = star_list[i]
                    kernel_star = util.cut_edges(star, psf_size_large)
                    kernel_star = util.kernel_norm(kernel_star)
                    kernel_list.append(kernel_star-kernel_large)
            psf_kwargs = {'psf_type': psf_type, 'kernel': kernel, 'kernel_large': kernel_large, 'kernel_list': kernel_list}
        else:
            raise ValueError('psf type %s not in list' % psf_type)
        return psf_kwargs, restrict_psf, star_list

    def get_psf_from_fits(self, path2fits, psf_type, psf_size, psf_size_large=91):
        """
        ment to import a psf from Tiny Tim
        :param path2fits: path to the fits file
        :return:
        """
        psf_data = pyfits.getdata(path2fits)
        kernel = util.cut_edges_TT(psf_data, psf_size)
        kernel_large = util.cut_edges_TT(psf_data, psf_size_large)
        kernel_large = util.kernel_norm(kernel_large)
        kernel_large_norm = np.copy(kernel_large)
        kernel = util.kernel_norm(kernel)
        psf_kwargs = {'psf_type': psf_type, 'kernel': kernel, 'kernel_large':kernel_large_norm}
        return psf_kwargs

    def get_psf_from_system(self, image_name, psf_size, psf_size_large=91):
        """
        ment to import a psf from Tiny Tim
        :param path2fits: path to the fits file
        :return:
        """
        psf_data = self.system.get_psf_data(image_name)
        psf_kwargs = self.cut_psf(psf_data, psf_size, psf_size_large)
        return psf_kwargs



    def get_psf_errors(self, psf_kwargs, data_kwargs, star_list):
        """
        returns a error map of sigma prop Intensity for a stacked psf estimation
        :param psf_kwargs:
        :param star_list:
        :return:
        """
        psf_size = len(psf_kwargs['kernel_large'])
        kernel_mean = util.image2array(psf_kwargs['kernel_large'])
        weights = np.zeros(len(star_list))
        cov_i = np.zeros((psf_size**2,psf_size**2))
        num_stars = len(star_list)
        for i in range(0,num_stars):
            star_list_i = star_list[i].copy()
            star = util.cut_edges(star_list_i, psf_size)
            weights[i] = np.sum(star)
            rel_array = np.array([util.image2array(star)/weights[i]-kernel_mean])
            a = (rel_array.T).dot(rel_array)
            cov_i += a
        factor = 1./(num_stars -1)
        #weights_sum = sum(weights)
        sigma2_stack = factor*util.array2image(np.diag(cov_i))
        psf_stack = psf_kwargs['kernel_large'].copy()
        sigma2_stack_new = sigma2_stack# - (data_kwargs['sigma_background']**2/weights_sum)
        sigma2_stack_new[np.where(sigma2_stack_new < 0)] = 0
        psf_stack[np.where(psf_stack < data_kwargs['sigma_background'])] = data_kwargs['sigma_background']
        error_map = sigma2_stack_new/(psf_stack)**2
        #error_map[np.where(error_map < psf_stack**2/data_kwargs['reduced_noise'])] = 0
        # n = len(error_map)
        #error_map[(n-1)/2-1:(n-1)/2+2,(n-1)/2-1:(n-1)/2+2] += 0
        error_map = filters.gaussian_filter(error_map, sigma=0.5)
        return error_map
