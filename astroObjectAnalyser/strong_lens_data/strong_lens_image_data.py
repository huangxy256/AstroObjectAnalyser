__author__ = 'sibirrer'


#external modules
import astropy.io.fits as pyfits
import astropy.wcs as pywcs
import astropy.coordinates as coords
import numpy as np
import os
import astrofunc.util as util
from astrofunc.util import Util_class



#internal modules
from astroObjectAnalyser.DataAnalysis.analysis import Analysis


class StrongLensImageData(object):
    """
    contains all the information associated with a given band image (e.g.  r-band frame of specific lens)
    """

    def __init__(self, local_filename=None, local_psf_filename=None, local_wht_filename=None, ra=None, dec=None,
                 ra_cutout_cent=None, dec_cutout_cent=None, cutout_scale=None, data_type='cosmos'):

        """
        initialize data class with file names (absolute paths), coordinates and data type specifications
        """
        self.ra = ra
        self.dec = dec
        self.ra_cutout_cent = ra_cutout_cent
        self.dec_cutout_cent = dec_cutout_cent
        self.cutout_scale = cutout_scale

        self.local_filename = local_filename
        self.local_psf_filename = local_psf_filename
        self.local_wht_filename = local_wht_filename

        self.analysis = Analysis()
        self.util_class = Util_class()
        self.data_type = data_type


    @property
    def data_cutout(self):
        if not hasattr(self, '_data_cutout'):
            self.image_cutout(self.ra_cutout_cent, self.dec_cutout_cent, self.cutout_scale)
        return self._data_cutout

    def del_cutout(self):
        if hasattr(self, '_data_cutout'):
            del self._data_cutout
            del self._ra_coords_cutout
            del self._dec_coords_cutout
            del self._header_cutout

    @property
    def header_cutout(self):
        if not hasattr(self, '_header_cutout'):
            print('WARINING: New cutout image is built with default number of pixels.')
            self.image_cutout(self.ra_cutout_cent, self.dec_cutout_cent, self.cutout_scale)
        return self._header_cutout

    @property
    def header_primary(self):
        if not hasattr(self, '_header_primary'):
            self.header_info()
        return self._header_primary

    @property
    def header(self):
        if not hasattr(self, '_header'):
            self.header_info()
        return self._header

    @property
    def cd1(self):
        if not hasattr(self, '_cd1'):
            self.pixel_scale()
        return self._cd1

    @property
    def cd2(self):
        if not hasattr(self, '_cd2'):
            self.pixel_scale()
        return self._cd2

    @property
    def naxis1(self):
        if not hasattr(self, '_naxis1'):
            self.pixel_number()
        return self._naxis1

    @property
    def naxis2(self):
        if not hasattr(self, '_naxis2'):
            self.pixel_number()
        return self._naxis2

    @property
    def exposure_time(self):
        if not hasattr(self, '_exposure_time'):
            self._exposure_time = self.header_primary.get('EXPTIME')
        return self._exposure_time

    @property
    def exposure_map(self):
        if not hasattr(self, '_exposure_map'):
            self.image_cutout(self.ra_cutout_cent, self.dec_cutout_cent, self.cutout_scale, exposure_map=True)
        return self._exposure_map

    @property
    def CCD_gain(self):
        if not hasattr(self, '_CCD_gain'):
            if self.data_type == 'DES':
                self._CCD_gain = self.header.get('GAIN')
            else:
                self._CCD_gain = self.header_primary.get('CCDGAIN')
            print("ccd gain =", self._CCD_gain)
        return self._CCD_gain

    @property
    def background(self):
        if not hasattr(self, '_background'):
            self._background_mean, self._background_rms = self._get_background()
        return self._background_mean, self._background_rms

    @property
    def get_psf_data(self):
        if not hasattr(self, '_psf_data'):
            self._psf_data = self._get_psf_data()
        return self._psf_data

    def psf_fit(self, psf_type):
        private = '_' + psf_type
        if not hasattr(self, private):
            _, psf_variables = self._get_psf(psf_type)
            setattr(self, private, psf_variables)
        return getattr(self, private)

    @property
    def psf_kernel(self):
        if not hasattr(self, '_psf_kernel'):
            kernel, psf_variables = self._get_psf(psf_type='moffat')
            self._psf_kernel = kernel
            if not hasattr(self, '_moffat'):
                self._moffat = psf_variables
        return self._psf_kernel

    @property
    def get_cat(self):
        if not hasattr(self, '_cat'):
            cat = self._get_cat()
            self._cat = cat
        return self._cat

    @property
    def get_cutout_coords(self):
        if not hasattr(self, '_ra_coords_cutout') or not hasattr(self, '_dec_coords_cutout'):
            self.image_cutout(self.ra_cutout_cent, self.dec_cutout_cent, self.cutout_scale)
        return self._ra_coords_cutout, self._dec_coords_cutout

    @property
    def transforms(self):
        if not hasattr(self, '_pix2coord_transform') or not hasattr(self, '_coord2pix_transform'):
            self._transform()
        return self._pix2coord_transform, self._coord2pix_transform

    def get_subgrid(self, subgrid_res=2):
        if not hasattr(self, '_ra_coords_cutout') or not hasattr(self, '_dec_coords_cutout'):
            self.image_cutout(self.ra_cutout_cent, self.dec_cutout_cent, self.cutout_scale)
        cos_dec = np.cos(self.dec / 360 * 2 * np.pi)
        print('test', cos_dec)
        ra_coords_sub, dec_coord_sub = self.util_class.make_subgrid(self._ra_coords_cutout*cos_dec, self._dec_coords_cutout, subgrid_res)

        return ra_coords_sub, dec_coord_sub

    def _get_background(self):
        """

        :return: mean and rms value of background
        """
        HDUFile, _ = self.get_HDUFile()
        mean, rms = self.analysis.get_background(HDUFile)
        return mean, rms

    def _get_cat(self):
        """

        :return: sextractor catalogue
        """
        HDUFile, image_no_boarder = self.get_HDUFile()
        cat = self.analysis.get_source_cat(HDUFile)
        return cat

    def _get_psf(self, psf_type):
        """

        :param psf_type:
        :return:
        """
        exp_time = self.exposure_time
        HDUFile, image_no_boarder = self.get_HDUFile()
        mean, rms = self.analysis.get_background(HDUFile)
        cat = self.analysis.get_source_cat(HDUFile)
        kernel, mean_list, filter_object = self.analysis.get_psf(image_no_boarder, cat, mean, rms, exp_time, psf_type)
        return kernel, mean_list

    def get_HDUFile(self, force=False):
        if not(hasattr(self, '_HDUFile') and hasattr(self, '_image_no_border') and (not force)):
            image = self.image_full()
            exp_time = self.exposure_time
            CCD_gain = self.CCD_gain
            import pyextract.image_config as ImageConfig
            conf_args = ImageConfig.config_arguments(exp_time, CCD_gain)
            self._HDUFile, self._image_no_border = ImageConfig.get_source_cat(image, conf_args)
        return self._HDUFile, self._image_no_border

    def _get_psf_data(self):
        """
        loads in psf data stored in Tiny_Tim folder of the lens system
        :return:
        """
        return pyfits.getdata(self.local_psf_filename)

    def header_info(self):
        """
        reads in the header info and performs checks on whether the header has the right format to deal with
        """
        self._header_primary = pyfits.getheader(self.local_filename) # this is the primary header which does not contain general information
        file = pyfits.open(self.local_filename)
        if self.data_type == 'DES' or self.data_type == 'cosmos':
            self._header = file[0].header
        else:
            self._header = file['SCI'].header  # this is the header of the science image
        file.close()

    def pixel_scale(self):
        """
        returns the pixel scale of the image (units still unclear!!!)
        """
        if not hasattr(self,'header'):
            self.header_info()

        self._cd1 = self.header.get('CDELT1') if self.header.get('CDELT1') else np.sqrt(self.header.get('CD1_1')**2 + self.header.get('CD1_2')**2)
        self._cd2 = self.header.get('CDELT2') if self.header.get('CDELT2') else np.sqrt(self.header.get('CD2_1')**2 + self.header.get('CD2_2')**2)
        if self.cd1 is None or self.cd2 is None:
            raise Exception("Missing CD or CDELT keywords in header")

    def _transform(self, approx=False):
        """
        initializes the the matrix which transforms pixel to ra/dec
        """
        if not hasattr(self, 'header'):
            self.header_info()

        CD1_1 = self.header.get('CD1_1')*3600  # change in arc sec per pixel d(ra)/dx
        CD1_2 = self.header.get('CD1_2')*3600
        CD2_1 = self.header.get('CD2_1')*3600
        CD2_2 = self.header.get('CD2_2')*3600

        self._pix2coord_transform = np.array([[CD1_1, CD1_2], [CD2_1, CD2_2]])
        det = CD1_1*CD2_2 - CD1_2*CD2_1
        self._coord2pix_transform = np.array([[CD2_2, -CD1_2], [-CD2_1, CD1_1]])/det

    def pixel_number(self):
        """
        reads in number of pixel per axis for original image
        """
        if self.header['NAXIS'] > 2:
            raise TypeError("Too many (%i) dimensions!" % self.header['NAXIS'])
        self._naxis1 = self.header['NAXIS1']
        self._naxis2 = self.header['NAXIS2']

    def image_full(self):
        """
        array of one full band, do only use this function when really needed as images can be quite large
        """
        file = pyfits.open(self.local_filename)
        if self.data_type == 'cosmos' or self.data_type == 'DES':
            data_full = file[0].data
        else:
            data_full = file['SCI'].data
        file.close()
        return data_full

    def exposure_full(self):
        """
        array of one full band exposure time. do only use this function when really needed as images can be quite large
        """
        if self.data_type == 'cosmos':
            file = pyfits.open(self.local_wht_filename)
            exp_full = file[0].data
            print("separate exposure map loaded")
        else:
            file = pyfits.open(self.local_filename)
            if self.data_type == 'DES':
                exp_full = file[1].data
            else:
                exp_full = file['WHT'].data
        file.close()
        return exp_full

    def image_cutout(self, xc, yc, cutout_scale, cutout_filename=None, exposure_map=False):
        """
        used to make a cut out image. Can also be called after the instance has been initialised, in which case
        the cutout image will be updated (i.e. changed).
        """
        # if not hasattr(self,'header'):
        #     self.header_info()
        if isinstance(cutout_filename, str):
            file = pyfits.open(cutout_filename)
            if self.data_type == 'cosmos' or self.data_type == 'DES':
                self._header_cutout = file[0].header
                self._data_cutout = file[0].data
            else:
                self._header_cutout = file['SCI'].header
                self._data_cutout = file['SCI'].data
            file.close()
        else:
            if cutout_scale is None:
                cutout_scale = 50
                print("New cutout is generated with default cutout scale.")
            xw, yw = cutout_scale, cutout_scale
            self._data_cutout, self._header_cutout, exp_map, self._ra_coords_cutout, self._dec_coords_cutout =\
                self._cutout(self.local_filename, xc, yc, xw, yw, exposure_map=exposure_map)
            if exposure_map:
                self._exposure_map = exp_map
        pass

    def _cutout(self, fits_filename, xc, yc, xw, yw, units='pixels',
            clobber=True, useMontage=False, coordsys='galactic', verbose=False, exposure_map=False):
        """
        Inputs:
            file  - pyfits HDUList (must be 2D)
            xc,yc - x and y coordinates in the fits files' coordinate system (CTYPE)
            xw,yw - x and y width (pixels or wcs)
            units - specify units to use: either pixels or wcs
            outfile - optional output file
        """
        # file = pyfits.open(fits_filename)
        # head = file['SCI'].header.copy()


        head = self.header.copy()
        wcs = pywcs.WCS(head)
        if units == 'wcs':
            if coordsys == 'celestial' and wcs.wcs.lngtyp == 'GLON':
                xc, yc = coords.Position((xc, yc), system=coordsys).galactic()
            elif coordsys == 'galactic' and wcs.wcs.lngtyp == 'RA':
                xc, yc = coords.Position((xc, yc), system=coordsys).j2000()
        xx, yy = wcs.all_world2pix(xc, yc, 0)
        print('the center of the image is at pixel coordinates %f, %f.' % (xx, yy))
        if units == 'pixels':
            xmin, xmax = np.max([0, xx-xw]), np.min([self.naxis1, xx+xw])
            ymin, ymax = np.max([0, yy-yw]), np.min([self.naxis2, yy+yw])
        elif units == 'wcs':
            xmin, xmax = np.max([0, xx-xw/np.abs(self.cd1)]), np.min([self.naxis1, xx+xw/np.abs(self.cd1)])
            ymin, ymax = np.max([0, yy-yw/np.abs(self.cd2)]), np.min([self.naxis2, yy+yw/np.abs(self.cd2)])
        else:
            raise Exception("Can't use units %s." % units)

        if xmax < 0 or ymax < 0:
            raise ValueError("Max Coordinate is outside of map: %f,%f." % (xmax,ymax))
        if ymin >= head.get('NAXIS2') or xmin >= head.get('NAXIS1'):
            raise ValueError("Min Coordinate is outside of map: %f,%f." % (xmin,ymin))

        head = self.change_header(head, xmin, xmax, ymin, ymax)
        img = self.image_full()[ymin:ymax, xmin:xmax].copy()
        # img = file['SCI'].data[ymin:ymax, xmin:xmax]
        if verbose: print("Cut image %s with dims %s to %s.  xrange: %f:%f, yrange: %f:%f" % (fits_filename, file['SCI'].data.shape, img.shape, xmin, xmax, ymin, ymax))
        if exposure_map:
            exp_map = self.exposure_full()[ymin:ymax, xmin:xmax].copy()
        else:
            exp_map = None
        ra_coord, dec_coord = self.get_coordinates(xmin, xmax, ymin, ymax, wcs)

        return img, head, exp_map, ra_coord, dec_coord

    def change_header(self, head, xmin, xmax, ymin, ymax):
        """
        changes the header to adjust information to the cutout image
        """
        head['CRPIX1'] -= xmin
        head['CRPIX2'] -= ymin
        head['NAXIS1'] = int(xmax-xmin)
        head['NAXIS2'] = int(ymax-ymin)
        if head.get('NAXIS1') == 0 or head.get('NAXIS2') == 0:
            raise ValueError("Map has a 0 dimension: %i,%i." % (head.get('NAXIS1'), head.get('NAXIS2')))
        return head

    def get_coordinates(self, x_min, x_max, y_min, y_max, wcs):
        """
        :param xc: center in ra
        :param yc: center in dec
        :param x_min: min pixel x-axis
        :param x_max: max pixel x-axis
        :param y_min: min pixel y-axis
        :param y_max: max pixel y-axis
        :param wcs: coordinate class initialized with the fits file of the original image
        :return: ra_coord, dec_coord in units of arc seconds centered to the cutout position
        """
        x_coords = np.linspace(x_min, x_max-1, x_max - x_min)
        y_coords = np.linspace(y_min, y_max-1, y_max - y_min)
        x_coords, y_coords = np.meshgrid(x_coords, y_coords)
        ra_coords, dec_coords = wcs.all_pix2world(x_coords, y_coords, 0)
        ra_coords -= self.ra
        dec_coords -= self.dec
        ra_coords *= 3600
        dec_coords *= 3600
        ra_coords = util.image2array(ra_coords)
        dec_coords = util.image2array(dec_coords)
        return ra_coords, dec_coords

    @property
    def get_pixel_zero_point(self):
        """

        :return: pixel coordinate (x_0, y_0) of (ra,dec) = (0,0)
        """
        head = self._header_cutout
        wcs = pywcs.WCS(head)
        x_0, y_0 = wcs.wcs_world2pix(self.ra, self.dec, 0)
        return x_0, y_0