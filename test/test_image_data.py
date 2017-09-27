__author__ = 'sibirrer'

"""
Tests for `StrongLensImageData` module.
"""

import os

import astrofunc.util as util
import numpy as np
import numpy.testing as npt
import pytest

from astroObjectAnalyser.image_data import StrongLensImageData


class TestStrongLensImageData(object):

    def setup(self):
        local_path = os.getcwd()
        fits_path = 'test/Test_data/RXJ1131_1231_test.fits'
        self.fits_filename = os.path.join(local_path, fits_path)
        self.name = 'RXJ1131-1231'
        self.ra = 172.96421  # in degree ICRS coordinates
        self.dec = -12.533066  # in degree ICRS coordinates
        self.system1 = StrongLensImageData(local_filename=self.fits_filename, ra=self.ra, dec=self.dec,
                                           ra_cutout_cent=self.ra, dec_cutout_cent=self.dec,  cutout_scale=None)
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
        img, head, exp_map, ra_coords, dec_coords = self.system1._cutout(self.fits_filename, self.ra, self.dec, xw, yw)
        ra_coords = util.array2image(ra_coords)
        dec_coords = util.array2image(dec_coords)
        assert head['NAXIS1'] == 2*xw
        assert head['NAXIS2'] == 2*yw
        assert len(img) == 2*yw
        assert len(img[0]) == 2*xw
        assert len(ra_coords) == yw*2
        assert len(ra_coords[0]) == xw*2
        cos_dec = np.cos(self.system1.dec / 360 * 2 * np.pi)
        print(cos_dec)
        print(np.sqrt((ra_coords[yw, xw]*cos_dec)**2 + dec_coords[yw, xw]**2))
        ra, dec = self.system1.map_pix2coord(0.5, 0.5)
        npt.assert_almost_equal(ra_coords[yw, xw], 0, decimal=8)
        npt.assert_almost_equal(dec_coords[yw, xw], 0, decimal=8)

    def test_image_full(self):
        """
        no testing of this routine
        """
        pass

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
        self.system1.image_cutout(self.ra, self.dec, cutout_scale=50, cutout_filename=cutout_filename)
        assert self.system1.header_cutout['NAXIS1'] == 348
        assert self.system1.header_cutout['NAXIS2'] == 349
        assert abs(self.system1.data_cutout[153, 173] - 1.1558878) < 0.0001

    def test_cutout_coordinates(self):
        """
        test routine when providing the cutout image already
        """
        filename = 'Test_data/RXJ1131_1231_74010_cutout.fits'
        path = os.path.dirname(__file__)
        cutout_filename = os.path.join(path, filename)
        self.system1.image_cutout(self.ra, self.dec, cutout_scale=50, cutout_filename=cutout_filename)

        pix2coord, coord2pix = self.system1.transforms
        from numpy.linalg import inv
        pix2coord_new = inv(coord2pix)
        npt.assert_almost_equal(pix2coord_new[0, 0], pix2coord[0, 0], decimal=8)
        npt.assert_almost_equal(pix2coord_new[0, 1], pix2coord[0, 1], decimal=8)
        npt.assert_almost_equal(pix2coord_new[1, 0], pix2coord[1, 0], decimal=8)
        npt.assert_almost_equal(pix2coord_new[1, 1], pix2coord[1, 1], decimal=8)

        ra, dec = self.system1.map_pix2coord(0, 0)
        ra_0, dec_0 = self.system1.coord_at_pixel_0
        assert ra == ra_0
        assert dec_0 == dec

        x_0, y_0 = self.system1.pixel_at_angle_0
        x, y = self.system1.map_coord2pix(0, 0)
        assert x_0 == x
        assert y_0 == y
        x, y = 50., 50.
        ra, dec = self.system1.map_pix2coord(x, y)
        x_, y_ = self.system1.map_coord2pix(ra, dec)
        npt.assert_almost_equal(x, x_, decimal=8)
        npt.assert_almost_equal(y, y_, decimal=8)

    def teardown(self):
        #tidy up
        print("tearing down " + __name__)
        pass

if __name__ == '__main__':
    pytest.main()