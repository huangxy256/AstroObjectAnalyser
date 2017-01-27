
"""
Tests for `lenstronomy` module.
"""
from __future__ import print_function, division, absolute_import

import pytest
import os
# from lenstronomy.strong_lens_data.data_manager import DataManager
from astroObjectAnalyser.strong_lens_data.strong_lens_system_factory import StrongLensSystemFactory as SLsys_fac
from mock import patch
from darkskysync.DataSourceFactory import DataSourceFactory
from darkskysync.DarkSkySync import DarkSkySync

STRONG_LENS_SYSTEMS_PATH = os.path.join(os.path.dirname(__file__), 'Test_data','strong_lens_systems.fits')

class TestLenstronomy(object):

    @patch.object(DataSourceFactory, "fromConfig", autospec=False)
    @patch.object(DarkSkySync, 'load', autospec=True, side_effect="files")
    def setup(self, dss_mock, dsf_object):
        
        #prepare unit test. Load data etc
        print("setting up " + __name__)
        test_dir = os.path.join(os.path.dirname(__file__ ))
        filename = os.path.join(test_dir,'test_data','RXJ1131_1231.sysdata')
        self.files = [filename,filename]
        
        dss_mock.side_effect = [[self.files]]
        
        self.lens_system_fac = SLsys_fac()

    @patch.object(DarkSkySync, 'load', autospec=True, side_effect=[[STRONG_LENS_SYSTEMS_PATH]])
    @patch.object(DarkSkySync, 'avail', autospec=True)
    def test_from_sysdata_files(self, avail_mock, load_mock):

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
        #tidy up
        print("tearing down " + __name__)
        pass

if __name__ == '__main__':
    pytest.main("-k TestLenstronomy")

