
__author__ = 'sibirrer'

from astroObjectAnalyser.strong_lens_data.strong_lens_image_data import StrongLensImageData
from astroObjectAnalyser.strong_lens_system_superclass import StrongLensSystem
from astroObjectAnalyser.strong_lens_data.data_manager import DataManager


class StrongLensSystemData(StrongLensSystem):
    """
    extends the features of StrongLensSystem. We should try to keep all the extensions as hidden functions
    (except for add_image_data)
    """

    def add_image_data(self, fits_filename, attrname, folder=None, image_dir=None, ra=None, dec=None, ra_cutout_cent=None, dec_cutout_cent=None,
                     data_manager=None, cutout_filename=None, cutout_scale=None, data_type='HST'):
        """
        """
        if not hasattr(self, 'available_frames'):
            self.available_frames = []
        if data_manager is None:
            self.data_manager = DataManager(image_dir=image_dir)
        else:
            self.data_manager = data_manager
        if folder is None:
            folder = self.name
        imagedata = StrongLensImageData(fits_filename, folder, ra=ra, dec=dec, ra_cutout_cent=ra_cutout_cent, dec_cutout_cent=dec_cutout_cent,
                 data_manager=self.data_manager,
                 cutout_filename=None, cutout_scale=None, data_type=data_type)

        setattr(self, attrname, imagedata)
        self.available_frames.append(attrname)


