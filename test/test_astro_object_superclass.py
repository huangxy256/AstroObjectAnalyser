
# Copyright (C) 2014 ETH Zurich, Institute for Astronomy

"""
Tests for `StrongLensRep` module.
"""
from __future__ import print_function, division, absolute_import, unicode_literals

import pytest
import os
from astroObjectAnalyser.astro_object_superclass import StrongLensSystem


class TestStrongLensSystem(object):

    def setup(self):
        #prepare unit test. Load data etc
        print("setting up " + __name__)
        self.system_base = StrongLensSystem('RXJ1131-1231')#, '11:31:51.6', '-12:31:57')
        self.system_base.add_info_attribute('ra_str','11:31:51.6')#, , '-12:31:57')
        self.system_base.add_info_attribute('dec_str','-12:31:57')#, , '-12:31:57')

        test_dir = os.path.join(os.path.dirname(__file__ ))
        self.fits_filename = os.path.join(test_dir,'Test_data','RXJ1131_1231_test.fits')

        pass

    def test_init(self):
        name = 'RXJ1131-1231'
        ra = '11:31:51.6'
        dec = '-12:31:57'
        system1 = StrongLensSystem(name)#, ra, dec)
        assert system1.name == name
        # assert system1.ra == 11.531000000000006
        # assert system1.dec == -12.5325
        # assert os.path.isdir(system1.local_cache)


    def test_add_info_attribute(self):
        system = self.system_base
        system.add_info_attribute('num_images',5,replace=False)
        assert system.num_images == 5
        system.add_info_attribute('num_images',3,replace=True)
        assert system.num_images == 3
        system.add_info_attribute('num_images',3)
        assert system.num_images == 3
        #TODO fix the assert error test below
        try:
            system.add_info_attribute('num_images',5,replace=False)
            pytest.fail('Error message !!!')
        except TypeError as e:
            assert True

        system.add_info_attribute('ra_str','11:31:51.6',replace=True)
        # assert not hasattr(system,'ra')
        assert hasattr(system,'ra_str')
        assert system.ra_str == '11:31:51.6'

        system.add_info_attribute('dec_str','-12:31:57',replace=True)
        assert hasattr(system,'ra')
        assert hasattr(system,'dec_str')
        assert system.dec_str == '-12:31:57'
        assert system.ra == 172.96499999999997
        assert system.dec == -12.5325

    def teardown(self):
        #tidy up
        print("tearing down " + __name__)
        pass

if __name__ == '__main__':
    pytest.main()

