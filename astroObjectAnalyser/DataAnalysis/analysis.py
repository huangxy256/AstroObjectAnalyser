__author__ = 'sibirrer'

import numpy as np
import scipy.ndimage.interpolation as interp
#import matplotlib.pyplot as plt

from astroObjectAnalyser.DataAnalysis.psf_fitting import Fitting
import astrofunc.util as util
from lenstronomy.FunctionSet.moffat import Moffat
from lenstronomy.FunctionSet.gaussian import Gaussian


class Analysis(object):
    """
    class which analyses data and fits a psf
    """

    def get_source_cat(self, HDUFile):
        """

        :param HDUFile:
        :return: catalogue
        """
        return HDUFile[2]

    def get_background(self, HDUFile):
        """
        filters the mean and rms value of the background computed by sextractor
        :param cat:
        :return: mean, rms
        """
        mean_found, rms_found = False, False
        list = HDUFile[1].data[0][0]
        for line in list:
            line = line.strip()
            line = line.split()
            if line[0] == 'SEXBKGND' or line[0] == 'SEXBKGND=':
                mean = float(line[1])
                mean_found = True
            if line[0] == 'SEXBKDEV' or line[0] == 'SEXBKDEV=':
                rms = float(line[1])
                rms_found = True
        if mean_found == False or rms_found == False:
            raise ValueError('no mean and rms value found in list.')
        return mean, rms

    def estimate_star_thresholds(self, cat):
        """
        estimates the cuts in the different sextractor quantities
        :param cat:
        :return:
        """
        mag = np.array(cat.data['MAG_BEST'],dtype=float)
        size = np.array(cat.data['FLUX_RADIUS'],dtype=float)
        #ellipticity = cat.data['ELLIPTICITY']

        kwargs_cuts = {}
        mag_max = min(np.max(mag), 34)
        mag_min = np.min(mag)
        delta_mag = mag_max - mag_min
        kwargs_cuts['MagMaxThresh'] = mag_max - 0.7*delta_mag
        kwargs_cuts['MagMinThresh'] = mag_min #+ 0.01*delta_mag

        mask = (mag<mag_max-0.5*delta_mag)
        kwargs_cuts['SizeMinThresh'] = max(0, np.min(size[mask]))
        kwargs_cuts['SizeMaxThresh'] = max(0, np.min(size[mask])+4)
        kwargs_cuts['EllipticityThresh'] = 0.1
        kwargs_cuts['ClassStarMax'] = 1.
        kwargs_cuts['ClassStarMin'] = 0.5
        return kwargs_cuts

    def estimate_galaxy_thresholds(self, cat):
        """
        estimates the cuts for selecting some galaxies
        :param selfself:
        :return:
        """
        mag = np.array(cat.data['MAG_BEST'],dtype=float)
        size = np.array(cat.data['FLUX_RADIUS'],dtype=float)
        ellipticity = cat.data['ELLIPTICITY']

        kwargs_cuts = {}
        mag_max = min(np.max(mag), 34)
        mag_min = np.min(mag)
        delta_mag = mag_max - mag_min
        kwargs_cuts['MagMaxThresh'] = mag_max - 0.5*delta_mag
        kwargs_cuts['MagMinThresh'] = mag_min + 0.2*delta_mag

        mask = (mag>mag_max-0.5*delta_mag)
        kwargs_cuts['SizeMinThresh'] = np.min(size[mask])*5
        kwargs_cuts['SizeMaxThresh'] = np.max(size[mask])
        kwargs_cuts['EllipticityThresh'] = 0.9
        kwargs_cuts['ClassStarMax'] = 0.1
        kwargs_cuts['ClassStarMin'] = 0.
        return kwargs_cuts

    def find_objects(self, cat, kwargs_cut):
        """

        :param cat: hdu[2] catalogue objects comming from sextractor
        :return: selected objects in the catalogue data list
        """
        mag = np.array(cat.data['MAG_BEST'],dtype=float)
        size = np.array(cat.data['FLUX_RADIUS'],dtype=float)
        ellipticity = cat.data['ELLIPTICITY']
        classStar = cat.data['CLASS_STAR']
        SizeMaxThresh = kwargs_cut['SizeMaxThresh']
        SizeMinThresh = kwargs_cut['SizeMinThresh']
        EllipticityThresh = kwargs_cut['EllipticityThresh']
        MagMaxThresh = kwargs_cut['MagMaxThresh']
        MagMinThresh = kwargs_cut['MagMinThresh']
        ClassStarMax = kwargs_cut['ClassStarMax']
        ClassStarMin = kwargs_cut['ClassStarMin']

        mask = (size<SizeMaxThresh) & (ellipticity<EllipticityThresh) & (size>SizeMinThresh) & (mag<MagMaxThresh) & (mag>MagMinThresh) & (classStar<ClassStarMax) & (classStar>ClassStarMin)
        return mask

    def get_objects_image(self, image, cat, mask, cut_radius=10, cut_fixed=None):
        """
        returns all the cutouts of the locations of the selected objects
        :param image:
        :param cat:
        :param mask:
        :return:
        """
        nx, ny = image.shape
        x_center = np.array(cat.data['X_IMAGE'],dtype=float)
        y_center = np.array(cat.data['Y_IMAGE'],dtype=float)
        size = np.array(cat.data['FLUX_RADIUS'],dtype=float)
        size_mask = size[mask]
        x_center_mask = x_center[mask]
        y_center_mask = y_center[mask]
        num_objects = len(x_center_mask)
        list = []

        for i in range(num_objects):
            xc, yc = x_center_mask[i], y_center_mask[i]
            if cut_fixed == None:
                numPix = size_mask[i]*cut_radius
            else:
                numPix = cut_fixed
            if (xc-numPix > 0) and (xc+numPix < nx) and (yc-numPix > 0) and (yc+numPix < ny):
                cutout = image[xc-numPix:xc+numPix+1,yc-numPix:yc+numPix+1]
                list.append(cutout)
        return list

    def get_psf(self, image, cat, mean, rms, poisson, psf_type = 'moffat', restrict_psf = None):
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
        star_list = self.get_objects_image(image, cat, mask, cut_radius=10, cut_fixed= 61)
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