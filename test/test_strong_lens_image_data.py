__author__ = 'sibirrer'

"""
Tests for `StrongLensImageData` module.
"""



from astroObjectAnalyser.strong_lens_data.strong_lens_image_data import StrongLensImageData
import astrofunc.util as util

import os
import pytest
import numpy as np
import numpy.testing as npt
from mock import patch
from test import IMAGE_DATA_PATH

class TestStrongLensImageData(object):

    @patch("darkskysync.DarkSkySync", autospec=False)
    def setup(self, dss_mock):
        
        print("setting up " + __name__)
        dss_mock.load_central_image_data.side_effect = [IMAGE_DATA_PATH]

        self.fits_filename = 'RXJ1131_1231_test.fits'
        self.name = 'RXJ1131-1231'
        self.ra_cutout_cent = 172.96421#274.40179 #in degree ICRS coordinates
        self.dec_cutout_cent = -12.533066#45.881676 #in degree ICRS coordinates
        self.system1 = StrongLensImageData(self.fits_filename, self.name, self.ra_cutout_cent, self.dec_cutout_cent,
                                self.ra_cutout_cent, self.dec_cutout_cent
                              ,data_manager=dss_mock, cutout_filename=None, cutout_scale=None)
        self.system1.data_type = 'HST'
        #prepare unit test. Load data etc
        # filename = 'Test_data/RXJ1131_1231_test.fits'
        # path = os.path.dirname(__file__)
        # self.fits_filename = os.path.join(path, filename)
        # cutout_filename = None
        # cutout_scale = None
        # self.system1 = StrongLensImageData(self.fits_filename,self.ra_cutout_cent,self.dec_cutout_cent, cutout_filename, cutout_scale)


        pass

    def test_init(self):
        """
        tests initialization without a given cutout image
        """

        assert self.system1.header_primary['TELESCOP'] == 'HST'
        assert self.system1.header['WCSNAME'] == 'DRZWCS'
        assert self.system1.cd1 == 1.3887809127949737e-05
        assert self.system1.cd2 == 1.3887808930277022e-05
        assert self.system1.naxis1 == 400
        assert self.system1.naxis2 == 400

    def test_cutout(self):
        """
        tests the cutout routine and whether the header is also changed accordingly
        """
        xw = 50  # number of pixel to be cutout in x-axis
        yw = 50  # number of pixel to be cutout in y-axis
        img, head, exp_map, ra_coords, dec_coords = self.system1._cutout(self.fits_filename, self.ra_cutout_cent, self.dec_cutout_cent, xw, yw)
        ra_coords = util.array2image(ra_coords)
        dec_coords = util.array2image(dec_coords)
        assert head['NAXIS1'] == 2*xw
        assert head['NAXIS2'] == 2*yw
        assert len(img) == 2*yw
        assert len(img[0]) == 2*xw
        assert len(ra_coords) == yw*2
        assert len(ra_coords[0]) == xw*2
        npt.assert_almost_equal(ra_coords[yw, xw], 0, decimal=8)
        npt.assert_almost_equal(dec_coords[yw, xw], 0, decimal=8)

    def test_image_full(self):
        """
        no testing of this routine
        """

    def test_make_sugrid(self):
        ra_sub, dec_sub = self.system1.get_subgrid(subgrid_res=2)
        ra_sub = util.array2image(ra_sub)
        ra, dec = self.system1.get_cutout_coords
        cos_dec = np.cos(self.system1.dec / 360 * 2 * np.pi)
        ra = util.array2image(ra*cos_dec)
        ra_resized = util.averaging(ra_sub, len(ra_sub), len(ra_sub)/2)
        assert ra[0][0] == ra_resized[0][0]
        assert ra[9][7] == ra_resized[9][7]


    def test_image_cutout(self):
        """
        test routine when providing the cutout image already
        """
        filename = 'Test_data/RXJ1131_1231_74010_cutout.fits'
        path = os.path.dirname(__file__)
        cutout_filename = os.path.join(path, filename)
        self.system1.image_cutout(self.ra_cutout_cent, self.dec_cutout_cent, cutout_scale=50, cutout_filename=cutout_filename)
        assert self.system1.header_cutout['NAXIS1'] == 348
        assert self.system1.header_cutout['NAXIS2'] == 349
        assert abs(self.system1.data_cutout[153, 173] - 1.1558878) < 0.0001

    def teardown(self):
        #tidy up
        print("tearing down " + __name__)
        pass

if __name__ == '__main__':
    pytest.main()