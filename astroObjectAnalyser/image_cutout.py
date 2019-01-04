
#  from astroObjectAnalyser.image_data_new import ImageData
import astropy.units as u
from astropy.coordinates import Angle
import numpy as np


class Frame(object):
    """
    manages cutouts of a image defined in the ImageData class with a joint coordinate system
    """
    def __init__(self, ra_0, dec_0):
        """

        :param imageData: ImageData instance
        :param ra_0: center of relative coordinate system
        :param dec_0: center of relative coordinate system
        """
        self._ra_0, self._dec_0 = self.convert_angle_units(ra_0, dec_0)
        self._cos_dec = np.cos(self._dec_0 / 360 * 2 * np.pi)

    def cutout(self, imageData, ra_c, dec_c, width, relative_coords=False):
        """

        :param ra_c: center of cutout
        :param dec_c: center of cutout
        :param width: width of cutout in arcseconds
        :param relative_coords: bool, if True, uses relatve coordinates (in arc seconds) from the defined center of the
        coordinate frame, else uses absolute coordinates
        :return:
        """
        ra_c, dec_c = self.convert_angle_units(ra_c, dec_c)
        xmin, xmax, ymin, ymax = imageData.cutout_range(ra_c, dec_c, width, width, units='arcseconds')
        img = imageData.image_full[int(ymin):int(ymax), int(xmin):int(xmax)].copy()
        wht_map = imageData.exposure_full[int(ymin):int(ymax), int(xmin):int(xmax)].copy()
        ra_at_xy_0, dec_at_xy_0 = imageData.pix2coord(int(xmin), int(ymin))
        ra_at_xy_0 -= self._ra_0
        ra_at_xy_0 *= self._cos_dec * 3600
        dec_at_xy_0 -= self._dec_0
        dec_at_xy_0 *= 3600
        xc = (xmax + xmin)/2
        yc = (ymax + ymin)/2
        transform_pix2angle, transform_angle2pix = imageData.transform(xc, yc)
        return img, wht_map, ra_at_xy_0, dec_at_xy_0, transform_pix2angle

    def convert_angle_units(self, ra, dec):
        """
        convert hexadecimal angle units into degrees

        :param ra:
        :param dec:
        :return:
        """
        if type(ra) is str:
            if (':' in ra) and (':' in dec): #assumes e.g 12:21:32.5 convention as 12h21min32.5s if ':' exists in the object
                angle_ra = Angle(ra, unit=u.hour)
                angle_dec = Angle(dec, unit=u.degree)
                ra = angle_ra.degree
                dec = angle_dec.degree
        return ra, dec
