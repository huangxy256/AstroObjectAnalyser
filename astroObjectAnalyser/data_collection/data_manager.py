import datetime
import time as time_mod


__author__ = 'amaraa, sibirrer'

import configparser
import os
import glob
import astrofunc.util as util
from collections import namedtuple
import astropy.io.fits as pyfits
import shutil


class DataManager(object):
    """
    manages the data for strong lensing systems. Acts as a layer between the
    strong_lens_system_factory and the way the data is stored in the background.

    """

    def __init__(self, server_path=None, scratch_path=None, directory_path=None):
        """

        :param server_path: path to mounted server
        :param scratch_path: path to local scratch (can be deleted, use hidden folders)
        :param directory_path: path on server to system database
        """

        self.lens_system_data = namedtuple('strong_lens_system', 'data_cat')
        self.server_path = server_path
        self.scratch_path = scratch_path
        self.directory_path = directory_path
        self._dir_local = os.path.join(self.scratch_path, self.directory_path)
        self._dir_server = os.path.join(self.server_path, self.directory_path)

        self.fits_file_path = 'strong_lens_systems.fits'
        self._max_size = 10**9  # 1Gbyte as max file size

    def _check_central_dir_access(self):
        if not os.path.isdir(self._dir_server):
            print('The central (SAN) directory %s is not accessible'%self._dir_server)
            return False
        return True

    def _sync_load(self, relative_path, force=False):
        """
        checks the presentce of a file locally, if not, makes local copy of it from server and returns local path
        :param relatve_path: path and name of file relative to mounted server
        :return: local path to file (which is existing)
        """
        local_path = os.path.join(self.scratch_path, self.directory_path, relative_path)
        if not os.path.isfile(local_path) or force is True:
            self._copy2local(relative_path)
        return local_path

    def _copy2local(self, relative_path):
        """
        copies a file from server to the local scratch folder
        :param relative_path: relative path on server (including file name)
        :return:
        """
        self._check_central_dir_access()
        remote_path = os.path.join(self.server_path, self.directory_path, relative_path)
        if not os.path.isfile(remote_path):
            raise ValueError("%s does not exist or is not a file!" % remote_path)
        size = os.path.getsize(remote_path)
        if size > self._max_size:
            raise ValueError("The file %s  has %s bytes, which is more than the limit of %s" %
                             (remote_path, size, self._max_size))
        local_path = os.path.join(self.scratch_path, self.directory_path, relative_path)
        filepath = local_path
        try:
            with open(filepath) as f:
                pass
        except IOError as e:
            splitlocaldir = filepath.split(os.sep)
            splitlocaldir.remove(splitlocaldir[-1:][0])
            localdir = ""
            for item in splitlocaldir:
                localdir += item + os.sep
            if not os.path.exists(localdir):
                os.makedirs(localdir)
            #shutil.copyfile(sourcefile, filepath)
        shutil.copy2(remote_path, local_path)

    def get_data(self, inputname, datatype='sysdata_file'):

        datatype_list = ['sysdata_file', 'fits']
        assert datatype in datatype_list,'%s cannot be used. Supported names are: %s' % (datatype_list, datatype)
        filename = self._sync_load(inputname, force=False)
        if datatype == 'sysdata_file':
            data = self._from_sysdata_files(filename)
        elif datatype == 'fits':
            data = self._from_fits(filename)
        else:
            raise ValueError("datatype %s not valid!" % datatype)
        return data

    def get_data_central(self, force=False):
        """
        access central strong_lens_system.fits file and reads it out
        :return:
        """

        filename = self._sync_load(self.fits_file_path, force=force)
        data = self._from_fits(filename)
        return data

    def load_central_image_data(self, folder, fits_image_name, force=False):
        """

        :param folder: folder relative to directory_path (normally named after the system)
        :param fits_image_name: name of the fits file
        :param force: force a new sync
        :return:
        """
        data_path = os.path.join(folder, fits_image_name)
        filename = self._sync_load(data_path, force=force)
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
            data_cat = util.dictionary_to_namedtuple(data_dict['catalog_data'])
            data_unit = self.lens_system_data(data_cat=data_cat)
            data_list.append(data_unit)

        return data_list

