__author__ = 'sibirrer'

import numpy as np


class Catalogue(object):
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
        mean, rms = False, False
        list = HDUFile[1].data[0][0]
        for line in list:
            line = line.strip()
            line = line.split()
            if line[0] == 'SEXBKGND' or line[0] == 'SEXBKGND=':
                mean = float(line[1])
            if line[0] == 'SEXBKDEV' or line[0] == 'SEXBKDEV=':
                rms = float(line[1])
        if mean is False or rms is False:
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
        x_center = np.array(cat.data['X_IMAGE'], dtype=float)
        y_center = np.array(cat.data['Y_IMAGE'], dtype=float)
        size = np.array(cat.data['FLUX_RADIUS'], dtype=float)
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