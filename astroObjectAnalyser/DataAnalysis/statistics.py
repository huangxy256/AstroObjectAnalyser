__author__ = 'sibirrer'


import numpy as np
import matplotlib.pyplot as plt # plotting package
import matplotlib.cm as cm # colormaps
import pyfits

import lenstronomy.util as util


class Statistics(object):
    """
    this class contains routines and methods to estimate the background noise of a given image
    """

    def plot_pixel_histogram(self,image):
        imgh = util.image2array(image)

        # print some statistics about the image
        print 'Image minimum = ', min(imgh)
        print 'Image maximum = ', max(imgh)
        print 'Image mean = ', np.mean(imgh)
        print 'Image standard deviation = ', np.std(imgh)

        # now plot a histogram of the image values
        plt.figure(2)
        plt.hist(imgh, bins=100, histtype='stepfilled')

        plt.show() # display the plots

        plow = -1.
        phi = 1.
        q = np.where((imgh >= plow) & (imgh <= phi))
        imghcut = imgh[q]

        print 'Image minimum = ', min(imghcut)
        print 'Image maximum = ', max(imghcut)
        print 'Image mean = ', np.mean(imghcut)
        print 'Image standard deviation = ', np.std(imghcut)

        return 0

