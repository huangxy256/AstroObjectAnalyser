
"""
Tests for `lenstronomy` module.
"""
from __future__ import print_function, division, absolute_import

import pytest
import os
# from lenstronomy.data_collection.data_manager import DataManager
from astroObjectAnalyser.data_collection.strong_lens_system_factory import StrongLensSystemFactory as SLsys_fac
from mock import patch


STRONG_LENS_SYSTEMS_PATH = os.path.join(os.path.dirname(__file__), 'Test_data','strong_lens_systems.fits')

class TestLenstronomy(object):

    def setup(self):
        
        #prepare unit test. Load data etc
        print("setting up " + __name__)
        test_dir = os.path.join(os.path.dirname(__file__))
        self.server_path = os.path.join(test_dir, 'Test_data')
        self.scratch_path = os.path.join(test_dir, 'Scratch_test')
        filename = os.path.join(test_dir,'test_data', 'RXJ1131_1231.sysdata')
        self.files = [filename, filename]
        
        self.lens_system_fac = SLsys_fac(server_path=self.server_path, scratch_path=self.scratch_path,
                                         directory_path="")

    def test_from_sysdata_files(self):

        # lens_systems = self.lens_system_fac.create_from_sysdata(self.files,datatype='sysdata_file')
        lens_systems = self.lens_system_fac.create_from_central()
        assert len(lens_systems) >= 7
        RXJ1131 = lens_systems[0]
        # assert len(RXJ1131.data_cat) == 6
        # assert len(RXJ1131.data_image) == 1
        #
        # assert RXJ1131.data_cat.name == 'RXJ1131-1231'
        # assert RXJ1131.data_cat.ra == '11:31:51.6'
        # assert RXJ1131.data_cat.dec == '-12:31:57'
        # assert RXJ1131.data_cat.radius_est == '1.7'
        # assert RXJ1131.data_cat.z_source == '0.658'
        # assert RXJ1131.data_cat.z_lens == '0.295'
        # assert RXJ1131.data_image.r_band == 'name_of_fits_file.fits'


    def teardown(self):
        os.remove(os.path.join(self.scratch_path, "strong_lens_systems.fits"))
        print("tearing down " + __name__)
        pass

if __name__ == '__main__':
    pytest.main("-k TestLenstronomy")

