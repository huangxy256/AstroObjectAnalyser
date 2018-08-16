__author__ = 'amaraa & sibirrer'

import astrofunc.constants as const
import astropy.units as u
import numpy as np
from astrofunc.Footprint.footprint import CheckFootprint
from astropy.coordinates import Angle


class StrongLensSystem(object):
    """
    contains all the information about one strong lens system.
    this will be inherited by other class. We should make
    sure that all the functionality that is exposed to the user
    is only defined here in the super class (especially those used to
    acess the data). Other functionality
    can be added by the subclass (hopefully in hidden function starting with _) - except
    add_image_data which should be replaced with more tailored ones.
    """

    def __init__(self, name):
        """
        initialise with very basic entries (that all the systems should/must have)
        """
        self.name = name
        self.tile_name = name
        self.available_frames = []


    def add_info_attribute(self, attrname, info_data, replace=False):
        """
        creates an attribute to store particular information.

        """
        attset = ['name','ra','dec','num_images','radius_est','z_source',
                  'sigma_z_source','z_lens','sigma_z_lens','sys_type','ra_str','dec_str','image_pos_str',
                  'image_pos_x','image_pos_y', 'tile_number', 'data_type']
        assert attrname in attset,'%s cannot be used. Supported names are: %s'%(attrname,attset)

        if attrname == 'sys_type':
            types = ['double', 'fold', 'cusp','ring','cross']
            assert info_data in types,'system type (%s) not in supported, please pick one of: %s' % (attrname, types)

        if replace or (not hasattr(self, attrname)):
            setattr(self, attrname, info_data)
        # elif not hasattr(self,'num_images'):
        #     setattr(self,attrname, info_data)
        elif not getattr(self, attrname) == info_data:
            raise TypeError("The number of images is already set to %f" %getattr(self, attrname))

        if attrname in ['ra_str','dec_str']:
            if (hasattr(self,'ra_str') and (hasattr(self,'dec_str'))):
                self.convert_angel_units()
        if attrname in ['image_pos','ra_str','dec_str']:
            if (hasattr(self,'ra') and (hasattr(self,'dec')) and (hasattr(self, 'image_pos_str'))):
                self.convert_image_pos(self.ra, self.dec)

    def is_in_survey(self, survey_name):
        """

        :param survey_name: survey name ATTENTION: only survey names accepted as implemented in astrofunc package
        :return: bool, True or False specifying whether the systems is in the survey requested.
        """
        checkFootprint = CheckFootprint()
        return checkFootprint.check_footprint(self.ra, self.dec, surveyname=survey_name)

    def add_image_data(self, imagedata, attrname):
        """
        simplest possible version of _add_image_data. Here we simply create an
        attribute for imagedata instance. *Can be replaced by subclass*
        """
        setattr(self, attrname, imagedata)
        self.available_frames.append(attrname)

    def add_image_data_init(self, attrname, local_filename=None, local_psf_filename=None, local_wht_filename=None,
                       ra=None, dec=None, ra_cutout_cent=None, dec_cutout_cent=None,
                       cutout_scale=None, data_type='HST', sci_extension=0, wht_extension=1):
        """
        adds an image_data class to the astro data object.
        Tries to take over all possible configurations of the object class
        """
        from astroObjectAnalyser.image_data import StrongLensImageData
        if ra is None:
            ra = self.ra
        if dec is None:
            dec = self.dec
        if ra_cutout_cent is None and hasattr(self, "ra_cutout_cent"):
            ra_cutout_cent = self.ra_cutout_cent
        if dec_cutout_cent is None and hasattr(self, "dec_cutout_cent"):
            dec_cutout_cent = self.dec_cutout_cent
        imagedata = StrongLensImageData(local_filename=local_filename, local_psf_filename=local_psf_filename,
                                        local_wht_filename=local_wht_filename, ra=ra, dec=dec,
                                        ra_cutout_cent=ra_cutout_cent, dec_cutout_cent=dec_cutout_cent,
                                        cutout_scale=cutout_scale, data_type=data_type, sci_extension=sci_extension,
                                        wht_extension=wht_extension)

        setattr(self, attrname, imagedata)
        self.available_frames.append(attrname)

    def merge(self, strong_lens_system, replace=False):
        """
        merges with the another instance of the same system. For example if you are working with an instance from
        the central store (i.e. on the SAN) and you want to update it with data in another instance (for the same system)
        that you have created elsewhere, the new instance (which contains new data) can be passed to the other one,
        which will update itself.

        :param strong_lens_system: another instance of strong_lens_system (with the same name attribute as self)
        :param replace: if set to true that attributes that are in both, the one from the input will be used to replace self
        """

        assert self.name == strong_lens_system.name, 'can only merge data for systems with the same name'

        keys = strong_lens_system.__dict__.keys()
        keys.remove('name')

        for attr in keys:
            if not hasattr(self, attr):
                setattr(self, attr, getattr(strong_lens_system, attr))
            elif getattr(self, attr) == getattr(strong_lens_system, attr):
                pass
            elif replace is True:
                setattr(self, attr, getattr(strong_lens_system, attr))
            else:
                raise(UserWarning,'duplicated attribute %s has not been changed. To force replace use replace = True'%attr)
        pass

    def set_data_type(self, attrname, data_type=None):
        """
        a data type is needed to query information in a .fits file. FITS files may have different structure.
        The datatype variable is meant to deal with it.
        :param attrname: data name
        :param data_type: string
        :return: data_type changed in attrname
        """
        if data_type is None:
            data_type = getattr(self, 'data_type')
        image_data_obj = getattr(self, attrname)
        image_data_obj.data_type = data_type

    def set_extension(self, attrname, ext_image=0, ext_wht=1):
        """

        :param ext_image:
        :param ext_wht:
        :return:
        """
        image_data_obj = getattr(self, attrname)
        image_data_obj.set_extension(ext_image, ext_wht)

    def get_full_image(self, attrname):
        """

        :param attrname: image name
        :return: full image
        """
        image_data_obj = getattr(self, attrname)
        return image_data_obj.image_full()

    def show_path2file(self, attrname):
        """

        :param attrname: image file name
        :return: path to image data and weight map
        """
        image_data_obj = getattr(self, attrname)
        image_data_obj.load_data()
        image_path = image_data_obj.local_filename
        try:
            wht_path = image_data_obj.local_wht_filename
        except:
            print("no wht file found")
            wht_path = None
        try:
            psf_path = image_data_obj.local_psf_filename
        except:
            print('no PSF file found')
            psf_path = None
        return image_path, wht_path, psf_path

    def get_cutout_image(self, attrname, cutout_scale, force=False):
        """

        :param attrname: image file name
        :param cutout_scale: number of pixels of cutout
        :param force: bool, if cutout already exists, only force=True will create another cutout
        :return: cutout image
        """
        image_data_obj = getattr(self, attrname)
        if force is True:
            image_data_obj.del_cutout()
        image_data_obj.cutout_scale = int(cutout_scale)
        return image_data_obj.data_cutout

    def del_cutout_image(self, attrname):
        """

        :param attrname: image file name
        :return: deletes cutout of attrname image
        """
        image_data_obj = getattr(self, attrname)
        image_data_obj.del_cutout()

    def get_cutout_header(self, attrname):
        """

        :param attrname: image file name
        :return: header of cutout image of attrname
        """
        image_data_obj = getattr(self, attrname)
        return image_data_obj.header_cutout

    def get_header(self, attrname):
        """

        :param attrname: image file name
        :return: header of full image attrname
        """
        image_data_obj = getattr(self, attrname)
        return image_data_obj.header

    def get_header_primary(self, attrname):
        """

        :param attrname: image file name
        :return: primary header (if present)
        """
        image_data_obj = getattr(self, attrname)
        return image_data_obj.header_primary

    def get_pixel_scale(self, attrname):
        """
        returns pixel scale
        :param attrname: image file name
        :return: returns pixel scale (unrotated) (if in arcseconds, devide by 3600!)
        """
        image_data_obj = getattr(self, attrname)
        return image_data_obj.cd1, image_data_obj.cd2

    def get_pixel_number(self, attrname):
        """

        :param attrname: image file name
        :return: number of pixels in axis1, axis2
        """
        image_data_obj = getattr(self, attrname)
        return image_data_obj.naxis1, image_data_obj.naxis2

    def get_coordinate_grid_absolute(self, attrname):
        """
        pixel coordinates in arc sec relative to the center (lens system coordinates)
        ATTENTION RA is scaled. Use get_coordinate_grid_relative to have relative arc sec coordinates.
        :param attrname:
        :return:
        """
        image_data_obj = getattr(self, attrname)
        ra_coords, dec_coords = image_data_obj.get_cutout_coords
        return ra_coords, dec_coords

    def get_coordinate_grid_linear(self, attrname):
        """
        returns the linear transformation coordinate grid that is fully consistent with the linear operation from pixel to coordinates
        :param attrname:
        :return:
        """
        image_data_obj = getattr(self, attrname)
        numPix = image_data_obj.numPix
        a = np.arange(numPix)
        matrix = np.dstack(np.meshgrid(a, a)).reshape(-1, 2)

        x_idex = matrix[:, 0]
        y_idex = matrix[:, 1]
        ra_coords, dec_coords = self.pix2coord(attrname, x_idex, y_idex)
        return ra_coords, dec_coords

    def get_coordinate_grid_relative(self, attrname):
        """
        relative RA, DEC coordinates in arcsec to the center (lens system coordinates)
        :param attrname:
        :return:
        """
        image_data_obj = getattr(self, attrname)
        ra_coords, dec_coords = image_data_obj.get_cutout_coords
        cos_dec = np.cos(self.dec/360*2*np.pi)
        return ra_coords*cos_dec, dec_coords

    def get_coordinate_subgrid(self, attrname, subgrid_res=2):
        """

        :param attrname: image file name
        :param subgrid_res: subgrid resolution
        :return: relative arc secs from center of the lens system
        """
        image_data_obj = getattr(self, attrname)
        ra_coords, dec_coords = image_data_obj.get_subgrid(subgrid_res)
        return ra_coords, dec_coords

    def get_image_position(self):
        """

        :return: coordinates relative of image positions, in relative arc sec (e.g. lensed quasars)
        """
        if hasattr(self, 'image_pos_ra') and hasattr(self, 'image_pos_dec'):
            return self.image_pos_ra, self.image_pos_dec
        else:
            raise ValueError('Strong Lens system does not provide image positions')

    def pix2coord(self, attrname, x_pix, y_pix):
        """

        :param attrname:
        :param x_pix:
        :param y_pix:
        :return: relative ra, dec for pixel coordinate (arc seconds)
        """
        image_data_obj = getattr(self, attrname)
        ra, dec = image_data_obj.map_pix2coord(x_pix, y_pix)
        return ra, dec

    def coord2pix(self, attrname, ra, dec):
        """

        :param attrname:
        :param ra:
        :param dec:
        :return: pixel coordinate of realtive ra, dec coordinate (arc seconds)
        """
        image_data_obj = getattr(self, attrname)
        x_pix, y_pix = image_data_obj.map_coord2pix(ra, dec)
        return x_pix, y_pix

    def get_exposure_time(self, attrname):
        """

        :param attrname: image file name
        :return: exposure time of image
        """
        image_data_obj = getattr(self, attrname)
        return image_data_obj.exposure_time

    def get_exposure_map(self, attrname):
        """

        :param attrname: image file name
        :return: effective exposure time for each pixel (wheight map)
        """
        image_data_obj = getattr(self, attrname)
        return image_data_obj.exposure_map

    def get_CCD_gain(self, attrname):
        """

        :param attrname: image file name
        :return: CCD gain
        """
        image_data_obj = getattr(self,attrname)
        return image_data_obj.CCD_gain

    def get_psf_kernel_image(self, attrname, kernel_size):
        """

        :param attrname: image file name
        :return: point spread function as a 2x2 kernel (as a fit from nearby stars in the same image)
        """
        image_data_obj = getattr(self, attrname)
        return image_data_obj.psf_kernel(kernel_size)

    def get_psf_fit(self, attrname, psf_type):
        """
        fits a parameterized model to stars nearby
        :param attrname:
        :param psf_type: fitting function (Moffat or Gaussian implemented)
        :return: fitting parameters
        """
        image_data_obj = getattr(self, attrname)
        return image_data_obj.psf_fit(psf_type)

    def get_psf_from_file(self, attrname, kernelsize):
        """

        :param attrname: image file name
        :return: returs PSF from data file (.fits file) 2x2 kernel
        """
        image_data_obj = getattr(self, attrname)
        return image_data_obj.get_psf_from_file(kernelsize)

    def get_background(self, attrname):
        """

        :param attrname: image file name
        :return: background sigma estimate from image
        """
        image_data_obj = getattr(self, attrname)
        return image_data_obj.background

    def get_cat(self, attrname):
        """

        :param attrname: image file name
        :return: SourceExtractor catalogue of image
        """
        image_data_obj = getattr(self, attrname)
        return image_data_obj.get_cat

    def shift_cutout_center(self, attrname, delta_ra, delta_dec):
        """

        :param attrname: image file name
        :param delta_ra: shift in RA [arcsec]
        :param delta_dec: shift in DEC [arcsec]
        :return: relocates center coordinate frame of lens system
        """
        image_data_obj = getattr(self, attrname)
        cos_dec = np.cos(self.dec / 360 * 2 * np.pi)
        image_data_obj.ra_cutout_cent = image_data_obj.ra + delta_ra/3600./cos_dec
        image_data_obj.dec_cutout_cent = image_data_obj.dec + delta_dec/3600.
        print("WARNING: This command may have unexpected consequences for the relative coordinate system!")

    def shift_system_center(self, delta_ra, delta_dec):
        """

        :param delta_ra:
        :param delta_dec:
        :return:
        """
        assert hasattr(self, 'ra')
        assert hasattr(self, 'dec')
        cos_dec = np.cos(self.dec / 360 * 2 * np.pi)
        self.ra += delta_ra/3600./cos_dec
        self.dec += delta_dec/3600.
        print("WARNING: This command may have unexpected consequences for the relative coordinate system!")

    def pixel_at_angle_0(self, attrname):
        """
        get pixel coordinate of lens system coordinate
        :param attrname: image file name
        :return: x_pix, y_pix
        """
        image_data_obj = getattr(self, attrname)
        x_at_radec_0, y_at_radec_0 = image_data_obj.pixel_at_angle_0
        return x_at_radec_0, y_at_radec_0

    def coord_at_pixel_0(self, attrname):
        """
        get relative ra,dec coordinate of pixel image 0,0 of cutout
        :param attrname: image file name
        :return: ra, dec
        """
        image_data_obj = getattr(self, attrname)
        ra_at_xy_0, dec_at_xy_0 = image_data_obj.coord_at_pixel_0
        return ra_at_xy_0, dec_at_xy_0

    def get_transform_matrix_angle2pix(self, attrname):
        """
        linear transformation matrix of relative angular [arcsec] to pixel coordinates
        :param attrname: image file name
        :return: 2x2 matrix
        """
        image_data_obj = getattr(self, attrname)
        pix2coord_transform, coord2pix_transform = image_data_obj.transforms
        return coord2pix_transform

    def get_transform_matrix_pix2angle(self, attrname):
        """
        linear transformation matrix of relative angular [arcsec] to pixel coordinates (inverse
        :param attrname: image file name
        :return: 2x2 matrix
        """
        image_data_obj = getattr(self, attrname)
        pix2coord_transform, coord2pix_transform = image_data_obj.transforms
        return pix2coord_transform

    def get_transform_matrix_undistorted(self, attrname):
        """
        linear transformation matrix of relative angular [arcsec] to pixel coordinates
        :param attrname: image file name
        :return: 2x2 matrix
        """
        image_data_obj = getattr(self, attrname)
        pix2coord_transform, coord2pix_transform = image_data_obj.transforms_undistorted
        return pix2coord_transform, coord2pix_transform

    def get_HDUFile(self, attrname, force=False):
        """

        :param attrname: image file name
        :param force: force a new creation
        :return: SourceExtractor raw output HDUFile
        """
        image_data_obj = getattr(self, attrname)
        HDUFile, image_no_border = image_data_obj.get_HDUFile(force)
        return HDUFile, image_no_border

    def convert_angel_units(self):
        """
        checks and converts angel units of the class lens_system. Standard is ICRS in degrees.
        """
        assert hasattr(self, 'ra_str')
        assert hasattr(self, 'dec_str')

        # pass
        # if not (hasattr(self,ra_str) and hasattr(self,dec_str)):
        #     raise Warning
        #     return None

        if (':' in self.ra_str) and (':' in self.dec_str): #assumes e.g 12:21:32.5 convention as 12h21min32.5s if ':' exists in the object
            angle_ra = Angle(self.ra_str, unit=u.hour)
            angle_dec = Angle(self.dec_str, unit=u.degree)
            self.ra = angle_ra.degree
            self.dec = angle_dec.degree
        # #TODO needs to be expanded to (all) possible conventions to convert

    def convert_image_pos(self, ra_cent, dec_cent):
        """
        catalogue image positions relative to center of system in arcsec
        :return: catalogue of arcsec in x and y for the multiple image system
        """
        arcsec = const.arcsec
        assert hasattr(self, 'image_pos_str')
        assert hasattr(self, 'ra')
        assert hasattr(self, 'dec')
        imPos = getattr(self, 'image_pos_str')
        self.image_pos_ra = np.empty_like((len(imPos)))
        self.image_pos_dec = np.empty_like((len(imPos)))
        for i in range(len(imPos)):
            if (':' in imPos[i, 0]) and (':' in imPos[i,1]): #assumes e.g 12:21:32.5 convention as 12h21min32.5s if ':' exists in the object
                angle_ra = Angle(imPos[i, 0], unit=u.hour)
                angle_dec = Angle(imPos[i, 1], unit=u.degree)
                self.image_pos_ra[i] = (ra_cent - angle_ra.degree)/arcsec*np.cos(angle_dec)
                self.image_pos_dec[i] = (dec_cent - angle_dec.degree)/arcsec

