import datetime
import time as time_mod


__author__ = 'amaraa'

import configparser
import os
import glob
import astrofunc.util as util
from collections import namedtuple
import pickle
from darkskysync.DarkSkySync import DarkSkySync
import shutil
import astropy.io.fits as pyfits
import numpy as np


class DataManager(object):
    """
    manages the data for strong lensing systems. Acts as a layer between the
    strong_lens_system_factory and the way the data is stored in the background.

    """

    def __init__(self, image_dir=None):

        self.lens_system_data = namedtuple('strong_lens_system', 'data_cat data_image')
        # if central_data_dir is None:
        self.dssync = DarkSkySync()
        self.dark_sky_data_dir = 'Strong_lensing/Lenstronomy_data/'
        self.dark_sky_cosmos_dir = 'acs_cosmos/v2/'
        if image_dir is None:
            self.dark_sky_image_dir = os.path.join(self.dark_sky_data_dir, 'image_data')
        else:
            self.dark_sky_image_dir = image_dir
        self.dark_sky_pickle_file = os.path.join(self.dark_sky_data_dir, 'central_pickle.pickle')
        self.dark_sky_fits_file = os.path.join(self.dark_sky_data_dir, 'strong_lens_systems.fits')

        self.central_data_dir = os.path.join('/Volumes/astro/refreg/data', self.dark_sky_data_dir)
        self.central_image_dir = os.path.join(self.central_data_dir, 'image_data')
        self.central_archive = os.path.join(self.central_data_dir,'archive')
        # self.central_pickle_file= os.path.join(self.central_data_dir,'central_pickle.pickle')
        # print(self.central_data_dir,'central_data_dir')
        # else:
        #     self.central_data_dir = central_data_dir

        # self.central_filename = 'strong_lens_systems.pickle'

    def get_data(self,inputname,datatype='sysdata_file'):

        datatype_list = ['sysdata_file','fits']

        assert datatype in datatype_list,'%s cannot be used. Supported names are: %s'%(datatype_list,datatype)
        if datatype == 'sysdata_file':
            data = self._from_sysdata_files(inputname)
        elif datatype == 'fits':
            data = self._from_fits(inputname)
        else:
            raise ValueError("datatype %s not valid!" % datatype)
        return data

    def get_data_central(self):

        filename = self.dssync.load(self.dark_sky_fits_file, force=True)[0]
        data = self._from_fits(filename)
        for i in range(0, len(data)):
            temp_filelist = self.dssync.avail(os.path.join(self.dark_sky_image_dir, data[i].data_cat.name))
            # print(temp_filelist)
            # print(data[i].data_cat)
            # print(data[i].data_image)
            # print(data.data_image)
            # data[i].data_image = 'temp'
            data[i].data_image.extend(temp_filelist)
        #TODO test this
        return data

    def get_pickle(self):
        """
        loads pickle file with lens system classes
        :return:
        """
        print(self.dark_sky_pickle_file)
        return pickle.load(self.dark_sky_pickle_file)

    def update_pickle(self, lens_system):
        """
        updates central pickle file
        :param lens_system:
        :return:
        """
        system_name = lens_system.name
        data_list = self.get_pickle()
        for i in range(0,len(data_list)):
            if data_list[i].name == system_name:
                data_list[i] = lens_system
                pickle.dump(data_list, self.dark_sky_pickle_file)
                return 0
        data_list.append(lens_system)
        pickle.dump(data_list, self.dark_sky_pickle_file)
        raise(Warning, 'new lens system added to  %s '%self.dark_sky_pickle_file)

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

    def add_new_systems_to_central(self,tbdata):
        """
        should append a fits table with new lensing systems to the central fits table and mkdir
        where the fits imaging data for each system can be stored.

        :param tbdata: data entry from a fits table, e.g. coming from pyfits.get_data(file)
        """

        if not self._check_central_dir_access():
            print('Warning: The SAN is not mounted')
            return

        #TODO make a directory for each entry
        #TODO make sure the columns match
        #TODO make a new fits table so that I can append the new entries
        #TODO write the new fits table (and archive the previous one)

        pass

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
            data_image = []
            data_unit = self.lens_system_data(data_cat=data_cat, data_image=data_image)
            data_list.append(data_unit)
        return data_list

    def _check_central_dir_access(self):
        if not os.path.isdir(self.central_data_dir):
            print('The central (SAN) directory %s is not accessible'%self.central_data_dir)
            return False
        return True

    def _rd_pickle_file(self,pickle_filename):
        fileObject = open(pickle_filename,'r')
        lens_sys_list = pickle.load(fileObject)
        fileObject.close()
        return lens_sys_list

    def _write_pickle_file(self,pickle_filename,lens_sys_list,safe=True):
        #TODO needs to be tested (Not done yet)
        if safe is True:
            if os.path.isfile(pickle_filename):
                archive = os.path.join(os.path.pardir(pickle_filename),'pickle_archive')
                if not os.path.isfile(archive):
                    os.mkdir(archive)
                shutil.copy(pickle_filename,os.path.join(archive,'pickle_archived_'.join(self._time_string())))
        fileObject = open(pickle_filename,'w')
        pickle.dump(lens_sys_list,fileObject)
        fileObject.close()

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
                search_string = os.path.join(inputname,'*'+file_ending)
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
            data_image = util.dictionary_to_namedtuple(data_dict['image_data'])
            data_unit = self.lens_system_data(data_cat=data_cat, data_image=data_image)
            data_list.append(data_unit)

        return data_list

    def _from_clerk_table(self):
        pass



    # def manage_fitsdata(self,fits_file):
    #     """
    #     tries to make a copy of a fits file on the SAN then makes the darksky version the working one.
    #     """
    #     if self._check_central_dir_access():
    #         shutil.copy(fits_file,self.central_image_dir)
    #         fits_avail = self.dssync.avail(self.dark_sky_image_dir)
    #         base_file=os.path.basename(fits_file)
    #         assert base_file in fits_avail,'Something has gone wrong, data should have been moved to darksky but ' \
    #                                        'does not seem to be available!'
    #         local_filename = self.dssync.load(os.path.join(self.dark_sky_image_dir,base_file))
    #         return local_filename
    #     else:
    #         raise Warning,'It looks like the SAN is not mounted'
    #         return fits_file



    # def add_system_central(self,strong_lens_system,replace=False):
    #     """
    #     function for adding an instance of the strong_lens_system class to a centrally (i.e. on the SAN)
    #     list stored with pickle.
    #     """
    #     print('Adam test3')
    #
    #     if not self._check_central_dir_access():
    #         return None
    #     print('Adam test2')
    #
    #     lens_sys_list = self.get_central_pickle()
    #     print('Adam test4')
    #
    #     system_names_exist = []
    #     for i in range(0,len(lens_sys_list)):
    #         system_names_exist.append(lens_sys_list[i].name)
    #
    #     for i in range (0,len(strong_lens_system)):
    #         name = strong_lens_system[i].name
    #         if not name in system_names_exist:
    #             lens_sys_list.append(strong_lens_system[i])
    #         else:
    #             ind = system_names_exist.index(name)
    #             lens_sys_list[ind].merg(strong_lens_system[i],replace=replace)
    #             # raise Warning,'The system %s already has an entry (#%i)'%(name,ind)
    #
    #     self._write_pickle_file(self.central_pickle_file,lens_sys_list)
    #
    # def get_central_pickle(self):
    #     # if self._check_central_dir_access():
    #     #     return None
    #     print('Adam test5')
    #     pickle_filename = self.dssync.load(self.dark_sky_pickle_file,force=True)[0]
    #     print('Adam test 6')
    #     print(pickle_filename)
    #     print(os.path.exists(pickle_filename))
    #
    #     lens_sys_list = []
    #
    #     # pickle_filename = os.path.join(self.central_data_dir,self.central_filename)
    #     if os.path.isfile(pickle_filename):
    #         lens_sys_list = self._rd_pickle_file(pickle_filename)
    #
    #     return lens_sys_list#,pickle_filename

