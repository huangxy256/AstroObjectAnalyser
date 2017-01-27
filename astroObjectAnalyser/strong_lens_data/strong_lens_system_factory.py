
__author__ = 'sibirrer'

# from astropy.coordinates import ICRS
# import os
# import configparser
from strong_lens_system import StrongLensSystemData
from strong_lens_image_data import StrongLensImageData
from astroObjectAnalyser.strong_lens_data.data_manager import DataManager
import os
import pickle
# import glob


class StrongLensSystemFactory(object):
    """
    able to create instances of StrongLensSystem in many ways (+ book keeping)
    """

    def __init__(self, image_dir=None):
        """
        filname of a file in a format that is readable by configparser
        """
        self.data_manager = DataManager(image_dir=image_dir)
        pass

    def create_from_sysdata(self,input_data,datatype='sysdata_file'):
        """
        creates a list of instance of StrongLensCatalogData from simple text files using configparser
        """

        data_list = self.data_manager.get_data(input_data, datatype='sysdata_file')
        lens_system_list = self.create_from_namedtuple(data_list)
        return  lens_system_list

    def create_from_central(self):
        """
        Uses the common data stored on the SAN for make instances of all the strong lens systems
        we have available.

        :return:
        """
        data_list = self.data_manager.get_data_central()
        return self.create_from_namedtuple(data_list)

    def create_from_pickle(self):
        """
        loads the lens system list from pickle file
        :return:
        """
        lens_system_list = self.data_manager.get_pickle()
        return lens_system_list

    def update_pickle(self, lens_system):
        """
        updates the pickle file
        :return:
        """
        # TODO make sure the file gets also updated locally
        self.data_manager.update_pickle(lens_system)
        print('lens system list updated')


    def create_from_namedtuple(self, data_list):

        lens_system_list = []

        for i in range(0, len(data_list)):

            lens_system = StrongLensSystemData('')
            data_cat = data_list[i].data_cat
            data_image = data_list[i].data_image
            # print(len(data_cat))

            for j in range(0, len(data_cat)):
                # print(data_cat)
                lens_system.add_info_attribute(data_cat._fields[j], data_cat[j], replace=True)
            lens_system.convert_angel_units()

            for j in range(0, len(data_image)):
                filename = data_image[j]
                name = os.path.splitext(filename)[0]
                #print(name)
                #print(data_image)
                attrname = name
                ra = lens_system.ra
                dec = lens_system.dec

                lens_system.add_image_data(filename, attrname, ra=ra, dec=dec, ra_cutout_cent=ra, dec_cutout_cent=dec,
                                           data_manager=self.data_manager,
                                           cutout_filename=None, cutout_scale=None)

            # print('Adam test: %s' %lens_system.name)

            lens_system_list.append(lens_system)
            #TODO don't forget to write unit tests

        return lens_system_list


    def find_from_central(self, system_name, pickle=False):
        """

        :param system_name: name of the strong lens system
        :return: class StrongLensSystem of system_name
        """
        if pickle==True:
            data_list = self.create_from_pickle()
        else:
            data_list = self.create_from_central()
        for i in range(0,len(data_list)):
            if data_list[i].name == system_name:
                return data_list[i]
        raise ValueError('\s is not in the list of systems on the central repository' %(system_name))

