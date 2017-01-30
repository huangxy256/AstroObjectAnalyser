
__author__ = 'sibirrer'

from astroObjectAnalyser.astro_object_superclass import StrongLensSystem
from astroObjectAnalyser.data_collection.data_manager import DataManager


class StrongLensSystemFactory(object):
    """
    able to create instances of StrongLensSystem in many ways (+ book keeping)
    """

    def __init__(self, server_path=None, scratch_path=None, directory_path=None):
        """
        filname of a file in a format that is readable by configparser
        """
        self.data_manager = DataManager(server_path=server_path, scratch_path=scratch_path,
                                        directory_path=directory_path)

    def create_from_sysdata(self, input_data, datatype='sysdata_file'):
        """
        creates a list of instance of StrongLensCatalogData from simple text files using configparser
        """

        data_list = self.data_manager.get_data(input_data, datatype=datatype)
        lens_system_list = self.create_from_namedtuple(data_list)
        return lens_system_list

    def create_from_central(self):
        """
        Uses the common data stored on the SAN for make instances of all the strong lens systems
        we have available.

        :return:
        """
        data_list = self.data_manager.get_data_central()
        return self.create_from_namedtuple(data_list)

    def create_from_namedtuple(self, data_list):
        """

        :param data_list: list of information about all the lens systems in the database
        :return: list of StrongLensSystem classes with all the informaiton of data_list incorporated
        """

        lens_system_list = []

        for i in range(0, len(data_list)):

            lens_system = StrongLensSystem('')
            data_cat = data_list[i].data_cat

            for j in range(0, len(data_cat)):
                # print(data_cat)
                lens_system.add_info_attribute(data_cat._fields[j], data_cat[j], replace=True)
            lens_system.convert_angel_units()
            lens_system_list.append(lens_system)
        return lens_system_list

    def find_from_central(self, system_name):
        """

        :param system_name: name of the strong lens system
        :return: class StrongLensSystem of system_name
        """
        data_list = self.create_from_central()
        for i in range(0, len(data_list)):
            if data_list[i].name == system_name:
                return data_list[i]
        raise ValueError('%s is not in the list of systems on the central repository' %(system_name))

