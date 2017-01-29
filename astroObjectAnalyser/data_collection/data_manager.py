import datetime
import time as time_mod


__author__ = 'amaraa, sibirrer'

import configparser
import os
import glob
import astrofunc.util as util
from collections import namedtuple
from darkskysync.DarkSkySync import DarkSkySync
import astropy.io.fits as pyfits


class DataManager(object):
    """
    manages the data for strong lensing systems. Acts as a layer between the
    strong_lens_system_factory and the way the data is stored in the background.

    """

    def __init__(self, image_dir=None):

        self.lens_system_data = namedtuple('strong_lens_system', 'data_cat')
        # if central_data_dir is None:
        self.dssync = DarkSkySync()
        self.dark_sky_data_dir = 'Strong_lensing/Lenstronomy_data/'
        self.dark_sky_cosmos_dir = 'acs_cosmos/v2/'
        if image_dir is None:
            self.dark_sky_image_dir = os.path.join(self.dark_sky_data_dir, 'image_data')
        else:
            self.dark_sky_image_dir = image_dir
        self.dark_sky_fits_file = os.path.join(self.dark_sky_data_dir, 'strong_lens_systems.fits')

        self.central_data_dir = os.path.join('/Volumes/astro/refreg/data', self.dark_sky_data_dir)
        self.central_image_dir = os.path.join(self.central_data_dir, 'image_data')
        self.central_archive = os.path.join(self.central_data_dir,'archive')

    def get_data(self, inputname, datatype='sysdata_file'):

        datatype_list = ['sysdata_file', 'fits']

        assert datatype in datatype_list,'%s cannot be used. Supported names are: %s' % (datatype_list, datatype)
        if datatype == 'sysdata_file':
            data = self._from_sysdata_files(inputname)
        elif datatype == 'fits':
            data = self._from_fits(inputname)
        else:
            raise ValueError("datatype %s not valid!" % datatype)
        return data

    def get_data_central(self):
        """
        access central strong_lens_system.fits file and reads it out
        :return:
        """

        filename = self.dssync.load(self.dark_sky_fits_file, force=True)[0]
        data = self._from_fits(filename)
        #TODO test this
        return data

    def load_central_image_data(self, folder, fits_image_name, force=False, data_type="HST"):
        if data_type == "cosmos":
            data_path = os.path.join(self.dark_sky_cosmos_dir, fits_image_name)
            print(data_path, 'data_path image')
        else:
            data_path = os.path.join(self.dark_sky_image_dir, folder, fits_image_name)
        filename = self.dssync.load(data_path, force=force)[0]
        return filename

    def load_central_psf_data(self, folder, fits_image_name, force=False):

        data_path = os.path.join(self.dark_sky_image_dir, folder, 'PSF_TinyTim', fits_image_name, 'result00.fits')
        print(data_path, 'data_path psf')
        filename = self.dssync.load(data_path, force=force)[0]
        return filename

    def _from_fits(self, filename):
        """
        retrieves information about strong lensing systems and returns the data in the form of a
        namedtuple
        """

        data_list = []
        tbdata = pyfits.getdata(filename)
        cat_nametuple = namedtuple('catalog_info', tbdata.names)

        for i in range(0, tbdata.size):
            data_cat = cat_nametuple(*tbdata[:][i])
            data_unit = self.lens_system_data(data_cat=data_cat)
            data_list.append(data_unit)
        return data_list

    def _check_central_dir_access(self):
        if not os.path.isdir(self.central_data_dir):
            print('The central (SAN) directory %s is not accessible'%self.central_data_dir)
            return False
        return True

    def _time_string(self):
        today = datetime.date.today()
        time_now = time_mod.localtime()
        datestring = str(today.year).zfill(4)+str(today.month).zfill(2)+str(today.day).zfill(2)
        timestring = str(time_now[3]).zfill(2) + str(time_now[4]).zfill(2) + str(time_now[5]).zfill(2)
        output_string = datestring+'_'+timestring
        return output_string

    def _from_sysdata_files(self, inputname):

        data_list = []
        file_list = []
        if not type(inputname) == type(['']):
            inputname = [inputname]

        if len(inputname) == 1:
            if os.path.isdir(inputname[0]):
                file_ending = '.sysdata'
                search_string = os.path.join(inputname, '*'+file_ending)
                file_list = glob.glob(search_string)
            else:
                file_list.extend(inputname)
        else:
            file_list.extend(inputname)

        for filename in file_list:

            assert os.path.isfile(filename), 'file %s could not be found' %filename

            systemdata = configparser.ConfigParser()

            systemdata.read(filename)
            data_dict = systemdata._sections
            print(data_dict, "data_dict")
            data_cat = util.dictionary_to_namedtuple(data_dict['catalog_data'])
            data_unit = self.lens_system_data(data_cat=data_cat)
            data_list.append(data_unit)

        return data_list

