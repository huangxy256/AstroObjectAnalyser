__author__ = 'amaraa'

from astropy.coordinates import ICRS, FK5, Angle

import astropy.units as u
import numpy as np
from astrofunc.Footprint.footprint import CheckFootprint



# from StrongLensRep.strong_lens_data.strong_lens_image_data import StrongLensImageData
# from lenstronomy.Footprint.footprint import CheckFootprint

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
    #TODO add a method that lets the instance merge with another stronglenssystem instance

    def __init__(self, name):
        """
        initialise with very basic entries (that all the systems should/must have)
        """
        self.name = name
        self.tile_name = name
        # self.ra = ra
        # self.dec = dec
        self.available_frames = []
        # self.data_manager = DataManager()
        # self.convert_angel_units()

        # self.local_cache = tempfile.mkdtemp()


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
            assert info_data in types,'system type (%s) not in supported, please pick one of: %s' %(attrname,types)

        if (replace) or (not hasattr(self, attrname)):
            setattr(self, attrname, info_data)
        # elif not hasattr(self,'num_images'):
        #     setattr(self,attrname, info_data)
        elif not getattr(self,attrname) == info_data:
            raise TypeError("The number of images is already set to %f" %getattr(self, attrname))

        if attrname in ['ra_str','dec_str']:
            if (hasattr(self,'ra_str') and (hasattr(self,'dec_str'))):
                self.convert_angel_units()
        if attrname in ['image_pos','ra_str','dec_str']:
            if (hasattr(self,'ra') and (hasattr(self,'dec')) and (hasattr(self, 'image_pos_str'))):
                self.convert_image_pos(self.ra, self.dec)

    def is_in_survey(self, survey_name):
        """
        returns True of False specifying whether the systems is in the survey requested.
        """
        #surveys = ['DES','DES_SV','VLT']
        surveys = ['DES']

        assert survey_name in surveys,'Survey requested (%s) not in database, please pick one of: %s' %(survey_name,surveys)

        checkFootprint = CheckFootprint()
        return checkFootprint.check_footprint(self.ra,self.dec,surveyname = survey_name)

    def add_image_data(self, imagedata, attrname):
        """
        simplest possible version of _add_image_data. Here we simply create an
        attribute for imagedata instance. *Can be replaced by subclass*
        """
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


    #TODO might want to think about creating a superclass of image_data that support and define these features (that can be replaced in subclass)


    def set_data_type(self, attrname, data_type=None):
        if data_type is None:
            data_type = getattr(self, 'data_type')
        image_data_obj = getattr(self, attrname)
        image_data_obj.data_type = data_type

    def get_full_image(self, attrname):
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
            wht_path = None
        return image_path, wht_path

    def get_cutout_image(self, attrname, cutout_scale, force=False):
        image_data_obj = getattr(self, attrname)
        if force is True:
            image_data_obj.del_cutout()
        image_data_obj.cutout_scale = int(cutout_scale/2)
        return image_data_obj.data_cutout

    def del_cutout_image(self, attrname):
        image_data_obj = getattr(self, attrname)
        image_data_obj.del_cutout()

    def get_cutout_header(self, attrname):
        image_data_obj = getattr(self, attrname)
        return image_data_obj.header_cutout

    def get_header(self, attrname):
        image_data_obj = getattr(self, attrname)
        return image_data_obj.header

    def get_header_primary(self, attrname):
        image_data_obj = getattr(self, attrname)
        return image_data_obj.header_primary

    def get_pixel_scale(self, attrname):
        image_data_obj = getattr(self, attrname)
        return image_data_obj.cd1, image_data_obj.cd2

    def get_pixel_number(self, attrname):
        image_data_obj = getattr(self, attrname)
        return image_data_obj.naxis1, image_data_obj.naxis2

    def get_coordinate_grid(self, attrname):
        image_data_obj = getattr(self, attrname)
        ra_coords, dec_coords = image_data_obj.get_cutout_coords
        return ra_coords, dec_coords

    def get_coordinate_grid_relative(self, attrname):
        image_data_obj = getattr(self, attrname)
        ra_coords, dec_coords = image_data_obj.get_cutout_coords
        cos_dec = np.cos(self.dec/360*2*np.pi)
        return ra_coords*cos_dec, dec_coords

    def get_coordinate_subgrid(self, attrname, subgrid_res=2):
        image_data_obj = getattr(self, attrname)
        ra_coords, dec_coords = image_data_obj.get_subgrid(subgrid_res)
        return ra_coords, dec_coords

    def get_image_position(self):
        if hasattr(self, 'image_pos_ra') and hasattr(self, 'image_pos_dec'):
            return self.image_pos_ra, self.image_pos_dec
        else:
            raise ValueError('Strong Lens system does not provide image positions')

    def get_exposure_time(self, attrname):
        image_data_obj = getattr(self, attrname)
        return image_data_obj.exposure_time

    def get_exposure_map(self, attrname):
        image_data_obj = getattr(self,attrname)
        return image_data_obj.exposure_map

    def get_CCD_gain(self, attrname):
        image_data_obj = getattr(self,attrname)
        return image_data_obj.CCD_gain

    def get_psf_kernel(self, attrname, filter_object=None):
        image_data_obj = getattr(self,attrname)
        return image_data_obj.psf_kernel

    def get_psf_fit(self, attrname, psf_type):
        image_data_obj = getattr(self,attrname)
        return image_data_obj.psf_fit(psf_type)

    def get_background(self, attrname):
        image_data_obj = getattr(self,attrname)
        return image_data_obj.background

    def get_psf_data(self, attrname):
        image_data_obj = getattr(self,attrname)
        return image_data_obj.get_psf_data

    def get_cat(self, attrname):
        image_data_obj = getattr(self,attrname)
        return image_data_obj.get_cat

    def shift_cutout_center(self, attrname, delta_ra, delta_dec):
        image_data_obj = getattr(self, attrname)
        image_data_obj.ra_cutout_cent = image_data_obj.ra + delta_ra/3600.
        image_data_obj.dec_cutout_cent = image_data_obj.dec + delta_dec/3600.

    def get_pixel_zero_point(self, attrname):
        image_data_obj = getattr(self, attrname)
        x0, y0 = image_data_obj.get_pixel_zero_point
        return x0, y0

    def get_transform_matrix_angle2pix(self, attrname):
        image_data_obj = getattr(self, attrname)
        pix2coord_transform, coord2pix_transform = image_data_obj.transforms
        return coord2pix_transform

    def get_HDUFile(self, attrname, force=False):
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
        arcsec = 4*np.pi/360/3600
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